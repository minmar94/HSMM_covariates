###### GOMPERTZ hazard (discrete) ######
gomp <- function(i, alpha, betai) exp(alpha + betai*i)

###### Transition probabilities (from hazards) between t and t+1 ######
p_k <- function(i, alpha, betai, beta_x, x_t) 1 - exp(-(gomp(i + 0.5, alpha, betai)*exp(as.numeric(x_t%*%beta_x))))

###### Bivariate density ######
density <- function(theta1, theta2, tau = c(0, 0, 0.1, 0.1, 0)){
  
  mu1 <- tau[1]
  mu2 <- tau[2]
  k1 <- tau[3]
  k2 <- tau[4]
  rho <- tau[5]
  
  const <- (1-rho^2)*(1-k1^2)*(1-k2^2)/(4*pi^2)
  c0 <- (1+rho^2)*(1+k1^2)*(1+k2^2)-8*abs(rho)*k1*k2
  c1 <- 2*(1+rho^2)*k1*(1+k2^2)-4*abs(rho)*(1+k1^2)*k2
  c2 <- 2*(1+rho^2)*k2*(1+k1^2)-4*abs(rho)*(1+k2^2)*k1
  c3 <- -4*(1+rho^2)*k1*k2+2*abs(rho)*(1+k1^2)*(1+k2^2)
  c4 <- 2*rho*(1-k1^2)*(1-k2^2)
  
  f <- const/(c0-c1*cos(theta1-mu1)-c2*cos(theta2-mu2)-c3*cos(theta1-mu1)*cos(theta2-mu2)-c4*(sin(theta1-mu1)*sin(theta2-mu2)))
  
  return(f)
  
}

###### Conditional density ######
copula.2given1 <- function(theta1, tau = c(0, 0, 0.1, 0.1, 0)){
  
  mu1=tau[1]
  mu2=tau[2]
  k1=tau[3]
  k2=tau[4]
  rho=tau[5]
  
  A1 <- matrix(c(-1,k1,-k1,1),2,2)
  A2 <- matrix(c(abs(rho),-k2,-abs(rho)*k2,1),2,2)
  A <- A1%*%A2
  tA <- t(A)
  q <- rho/abs(rho)
  
  zz1 <- cos(theta1)+1i*sin(theta1)
  mmu1 <- cos(mu1)+1i*sin(mu1)
  mmu2 <- cos(mu2)+1i*sin(mu2)
  zzz <- (zz1*Conj(mmu1))^q
  complex.reg <- -mmu2*(tA[1,1]*zzz+tA[1,2])/(tA[2,1]*zzz+tA[2,2])
  
  return(complex.reg)
}

###### Simulate density ######
copula.sim <- function(n, tau){
  
  theta1 <- circular::rwrappedcauchy(n, mu=circular::circular(tau[1]), rho=tau[3])
  theta2 <- c()
  for(i in 1:n){
    theta2[i] <- circular::rwrappedcauchy(1,
                                          mu=circular::circular(Arg(copula.2given1(theta1[i],tau=tau))),
                                          rho=Mod(copula.2given1(theta1[i],tau=tau)))
    theta1 <- ifelse(theta1>pi,theta1-2*pi,theta1)
    theta2 <- ifelse(theta2>pi,theta2-2*pi,theta2)
  }
  return(cbind(theta1,theta2))
}

############################################
## Wrapper for simulation
############################################
simulate_hsmm_torus <- function(n, K, delta, alpha, betai, beta_x, pars, xmat, omega, seed = 7777){
  
  if(is.null(seed)) seed <- sample(1111:111111,1)
  set.seed(seed)
  # Initialize multinomial matrix for the markov chain
  u <- matrix(0, ncol = K, nrow = n)
  # Simulate initial state
  u[1, ] <- rmultinom(n = 1, size = 1, prob = delta)
  state_lab <- which(u[1,]==1) # get the label
  
  # Random generation of omega
  # omega <- matrix(runif(K*K), K, K)
  # diag(omega) <- 0
  # for(k in 1:K) omega[k,-k] <- omega[k,-k]/sum(omega[k,-k]) 
  
  # Dati
  ymat <- matrix(NA, nrow = n, ncol = 2)
  ymat[1,] <- copula.sim(n = 1, tau = pars[state_lab[1],])
  
  i <- 1 # age of the state: assumption is that at the beginning we have observed a transition
  Jt <- c(1,rep(0,n-1)) # jump == 1 if there is a transition
  times <- 1:n # time points
  
  for(t in 2:n){
    
    # compute transition probabilities
    pki <- p_k(times[i], alpha = alpha[state_lab[t-1]], betai = betai[state_lab[t-1]], beta_x = beta_x[state_lab[t-1],], x_t = xmat[t-1,])
    pvec <- c(1 - pki, pki*omega[state_lab[t-1],-state_lab[t-1]])
    ru <- rmultinom(1, 1, pvec) # New state
    
    if(ru[1,1] == 1){ # If no transition, than updata the age of the state
      
      u[t, ] <- u[t-1, ] # new state
      state_lab[t] <- which(u[t,]==1) # new state label
      ymat[t, ] <- copula.sim(n = 1, tau = pars[state_lab[t],])
      i <- i + 1 # age of the state
      
    }else{
      
      Jt[t] <- 1 # record the transition
      lab_app <- c(state_lab[t-1], (1:K)[-which(u[t-1,]==1)]) 
      u[t, lab_app[ru == 1]] <- 1 # new state
      state_lab[t] <- which(u[t,]==1) # new state label
      ymat[t, ] <- copula.sim(n = 1, tau = pars[state_lab[t],])
      i <- 1 # age of the state
      
    }
  }
  # Tn is the transition times
  Tn <- which(Jt == 1)
  # N of transitions
  Nt <- cumsum(Jt)
  
  # Dwell times
  runs <- rle(state_lab)
  M <- max(runs$lengths) # Max dwell
  
  return(list(ymat = ymat, state = state_lab, U = u, Omega = omega, Jt = Jt, Tn = Tn, Nt = Nt, runs = runs, M = M))
  
}


#############################
# hazards: qui si costruisce p.array
haz <- function(xmat, beta.array, K, M){
  
  ld <- rep(M, K)
  x <- as.matrix(xmat)
  n <- nrow(xmat)
  d <- rep(1:K, ld)
  p.array <- array(0, dim = c(n-1, sum(ld)))
  
  for(t in 1:(n-1)){
    for(k in 1:K){
      pos <- which(d == k)
      p.array[t, pos] <- p_k(1:ld[k], alpha = beta.array[k,1], betai = beta.array[k,2], beta_x = beta.array[k,-c(1,2)], x_t = xmat[t,])
    }
  }
  return(p.array)
}

#############################
# Augmented transition matrix for the HMM representation of the HSMM
Gamma.f <- function(Omega, p.array, K, M){
  
  n <- nrow(p.array)
  ld <- rep(M, K)
  d <- rep(1:K,ld)
  p <- p.array
  Gamma <- array(0, dim = c(n, sum(ld), sum(ld)))    
  for(t in 1:n){
    for(k in 1:K){
      pos <- which(d==k)
      if(length(pos) >= 2){
        if(length(pos)==2){
          Gamma[t,min(pos):(max(pos)-1), (min(pos)+1):max(pos)] <- 1-p[t,pos][-length(p[t,pos])]
        } else{
          diag(Gamma[t,min(pos):(max(pos)-1), (min(pos)+1):max(pos)]) <- 1-p[t,pos][-length(p[t,pos])]
        }
      }
      Gamma[t,max(pos),max(pos)] <- 1-p[t,pos][length(p[t,pos])]
      
      for(i in (1:K)[-k]){
        Gamma[t,min(pos):max(pos),c(1,which(diff(d)==1)+1)[i]] <- p[t,pos]*Omega[k,i]
      }
    }
    
  }
  return(Gamma)
}

# EM initialization
EM.start <- function(y, x, K, M, dim.theta, help_mu = T, theta_in = NULL){
  
  n <- nrow(as.matrix(y)) # n. obs
  ld <- rep(M, K)        # n. states
  x <- as.matrix(x)
  q <- ncol(x) # n. covariates
  
  # Initialize beta array
  beta.array <- array(0, dim=c(K, q+2)) # n. states times n. covariates (+2 because of the intercept and the time-effect) 
  # 
  beta.array[,1] <- runif(K, -5, 0)
  beta.array[,2] <- runif(K, 0.01, 0.4)
  beta.array[,3:(q+2)] <- runif(K*q, -.5, 0.5)
  
  # Initialize theta array
  theta.array <- array(dim = c(K, dim.theta))
  
  # Initialize design matrix
  d <- rep(1:K, ld)
  d_matrix <- array(0, dim = c(sum(ld)*(n-1), q+2))
  d_matrix[,1] <- rep(d, n-1)
  d_matrix[,2] <- rep(unlist(lapply(1:K, function(k) 1:ld[k])), n-1)
  for(j in 1:q) d_matrix[,2+j] <- rep(x[1:(n-1),j], each = sum(ld))
  d_matrix <- as.data.frame(d_matrix)
  colnames(d_matrix) <- c("state", "time", paste0("V", 1:ncol(x)))
  
  # HMM matrix
  Omega <- matrix(runif(K*K), K, K)     
  diag(Omega) <- 0
  # Normalize omega by row
  for(k in 1:K) Omega[k,-k]<- Omega[k,-k]/sum(Omega[k,-k]) 
  
  # HSMM matrix
  p.array <- haz(x, beta.array, K, M)
  Gamma_mat <- Gamma.f(Omega, p.array, K, M)
  
  # Initialize probabilities
  post.pi <- LaplacesDemon::rdirichlet(n, alpha = rep(1/K, K))
  for(k in 1:K){
    trig.mom <- colSums(post.pi[,k]*(cos(y)+1i*sin(y)))/sum(post.pi[,k])
    theta.array[k,] <- c(Arg(trig.mom), Mod(trig.mom), 0)
  }
  
  # Help by fixing the random start of the circular averages means
  if(help_mu) theta.array[,1:2] <- theta_in
  
  # initialization time 1 probabilities
  Pi <- runif(K) 
  Pi <- Pi/sum(Pi)
  
  return(list(d_matrix = d_matrix, 
              ld = ld, 
              theta.array = theta.array, 
              beta.array = beta.array, 
              Pi = Pi, 
              Omega = Omega, 
              Gamma = Gamma_mat, 
              p.array = p.array))
  
}


# weighted loglikelihood
w.loglik <- function(theta, obs = cbind(0, 0), weights = 1){
  
  theta[3] <- (tanh(theta[3])+1)/2
  theta[4] <- (tanh(theta[4])+1)/2
  theta[5] <- tanh(theta[5])
  f <- sum(weights*log(density(obs[,1], obs[,2], theta)))
  return(-f)
  
}


# backward.forward
backward.forward <- function(fit, Gamma, Pi, K, M){
  
  ld <- rep(M, K)
  n <- nrow(fit)
  
  psi <- alpha <- Beta <- array(dim = c(n, sum(ld)))
  C	<- array(dim = n)
  
  post.pi.aug <- array(dim = c(n, sum(ld)))
  post.bi.pi.aug <- array(dim = c(n-1, sum(ld), sum(ld)))
  
  psi[1,] <- Pi
  psi_per_fit <- psi[1,]*fit[1,]
  C[1] <- sum(psi_per_fit)
  alpha[1,] <- (psi_per_fit)/C[1] 
  
  for(i in 2:n){
    psi[i, ] <- alpha[i-1,]%*%Gamma[i-1,,]
    psi_per_fit <- psi[i,]*fit[i,]
    C[i] <- sum(psi_per_fit)
    alpha[i,] <- (psi_per_fit)/C[i] 
    
  }
  
  l <- sum(log(C))
  
  Beta[n,] <- 1/C[n]
  
  for(i in (n-1):1){
    Beta[i, ] <- c(Gamma[i,,]%*%(fit[i+1,]*Beta[i+1,]))/C[i]
  }
  
  alphaBeta <- alpha*Beta
  post.pi.aug = alphaBeta/rowSums(alphaBeta)
  
  for(h in 1:sum(ld)){
    for(k in 1:sum(ld)){
      post.bi.pi.aug[,h,k] <- alpha[1:(n-1),h]*Gamma[,h,k]*fit[2:n,k]*Beta[2:n,k] 
    }
  }
  
  return(list(l = l, post.pi.aug = post.pi.aug, post.bi.pi.aug = post.bi.pi.aug))
}

##############################
# EM algorithm               #
# for parametric dwell time  #
##############################
cloglog.hsmm.EM <- function(y, x = NULL, init = NULL, max.iter = 50, tol = 10^-4, interact = F,
                            verbose = TRUE, dim.theta = NULL,
                            link = "cloglog"){
  
  if(is.null(init)) init <- EM.start(y, x, ld, dim.theta)
  # initialization
  ld <- init$ld
  K <- length(ld)
  M <- max(ld)
  d <- rep(1:K,ld)
  y <- as.matrix(y)
  x <- as.matrix(x)
  q <- ncol(x)
  Gamma_mat <- init$Gamma
  Omega <- init$Omega
  d_matrix <- init$d_matrix
  theta.array <- init$theta.array
  theta.array.step <- beta.array.step <- list()
  beta.array <- init$beta.array
  Pi <- init$Pi
  n <- nrow(y)
  fit <- array(dim = c(n, sum(ld)))
  
  llk <- c()
  old.llk <- -Inf
  dif <- Inf
  step <- euclid_beta <- 0
  while((dif > tol | euclid_beta > 0.01) & step<max.iter){
    step <- step+1
    
    for(k in 1:K){
      pos <- which(d==k)
      fit[,pos] <- density(y[,1], y[,2], theta.array[k,])
    }
    
    bw <- backward.forward(fit = fit, Gamma = Gamma_mat, Pi = Pi, K, M)
    llk[step] <- bw$l
    post.pi.aug <- bw$post.pi.aug
    post.bi.pi.aug <- bw$post.bi.pi.aug
    # Normalization of post.bi.pi.aug for numerical issues
    for(i in 1:(n-1)) post.bi.pi.aug[i,,] <- post.bi.pi.aug[i,,]/sum(post.bi.pi.aug[i,,])
    post.pi <- t(sapply(1:n, function(i){sapply(split(post.pi.aug[i,], d), sum)}))
    
    ###########
    # M -step #
    ###########
    
    # initial (augmented) probabilities
    Pi <- post.pi[1,]
    
    # update Omega (not necessary if K =2)
    for(h in 1:K){
      posh <- which(d==h)
      s <- sum(post.bi.pi.aug[,posh,-posh])
      for(k in (1:K)[-h]){
        posk <- which(d==k)
        Omega[h,k] <- sum(post.bi.pi.aug[,posh,posk])/s
      }
    }
    
    cases_fin <- noncases_fin <- c()
    for(t in 1:(n-1)){
      cases_int <- noncases_int <- c()
      for(h in 1:K){
        posh <- which(d==h)
        cases_int <- c(cases_int, rowSums(post.bi.pi.aug[t, posh,-posh]))
        noncases_int <- c(noncases_int, rowSums(post.bi.pi.aug[t, posh,posh]))
      }
      cases_fin <- c(cases_fin, cases_int)
      noncases_fin <- c(noncases_fin, noncases_int)
    }
    
    if(interact){
      frml_pt1 <- paste("d_matrix", colnames(d_matrix)[-c(1,2)], sep = "$")
      frml <- as.formula(paste0("cbind(cases_fin, noncases_fin) ~ factor(d_matrix$state)*d_matrix$time +",
                                paste0("factor(d_matrix$state)*", frml_pt1, collapse =  " + ")))
      mod <- suppressWarnings(glm(frml, family = binomial(link = link)))
      coefs <- coef(mod)
      ids_app <- c(1:(K+1), (K+q+2):(2*K+q))
      coefs_covar <- matrix(coefs[setdiff(1:length(coefs), ids_app)], nrow = q)
      for(k in 2:K) coefs[k] <- coefs[1] + coefs[k]
      for(j in (K+q+2):(2*K+q)) coefs[j] <- coefs[j] + coefs[K+1]
      for(p in 1:q){
        for(k in 2:K)
          coefs_covar[p, k] <- coefs_covar[p,1] + coefs_covar[p,k]
      }
      
      coefs <- c(coefs[ids_app], c(t(coefs_covar)))
    }else{
      frml_pt1 <- paste("d_matrix", colnames(d_matrix)[-c(1,2)], sep = "$")
      frml <- as.formula(paste0("cbind(cases_fin,noncases_fin) ~ factor(d_matrix$state)*d_matrix$time +",
                                paste0(frml_pt1, collapse =  " + ")))
      mod <- suppressWarnings(glm(formula = frml, family = binomial(link = link)))
      coefs <- coef(mod)
      for(k in 2:K) coefs[k] <- coefs[1] + coefs[k]
      for(j in (K+q+2):(2*K+q)) coefs[j] <- coefs[j] + coefs[K+1]
      coefs_covar <- rep(coefs[(K+2):((K+1+q))], each = K)
      ids_app <- c(1:(K+1), (K+q+2):(2*K+q))
      coefs <- c(coefs[ids_app], coefs_covar)
    }
    beta.array[,1:ncol(beta.array)] <- coefs
    p.array <- haz(x, beta.array, K, M)
    Gamma_mat <- Gamma.f(Omega, p.array, K, M)
    # update parameters toroidal density
    for(k in 1:K){
      opt <- optim(theta.array[k,], w.loglik, obs = y, weights = post.pi[, k])
      theta.array[k,] <- c(opt$par[1], opt$par[2], (tanh(opt$par[3])+1)/2, (tanh(opt$par[4])+1)/2, tanh(opt$par[5]))
    }
    theta.array.step[[step]] <- theta.array
    beta.array.step[[step]] <- beta.array
    if(step > 1) euclid_beta <- sqrt(sum((beta.array.step[[step]]-beta.array.step[[step-1]])^2))
    if(step > 1) dif <- abs((llk[step]-old.llk)/old.llk)
    if(is.nan(dif)) dif <- 0
    old.llk <- llk[step]
    if(verbose==TRUE) cat("iteration ", step,"; loglik = ", llk[step], "\n")
  }
  return(list(K = K, theta.array = theta.array, Pi = post.pi[1,], 
              Gamma = Gamma_mat, Omega = Omega,
              p.array = p.array, beta.array = beta.array, theta.array.step = theta.array.step, beta.array.step = beta.array.step,
              loglik = llk, iteration = step, post.pi = post.pi, 
              post.pi.aug = post.pi.aug, post.bi.pi.aug = post.bi.pi.aug))
}


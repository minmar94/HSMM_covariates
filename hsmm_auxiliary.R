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
simulate_hsmm_torus <- function(n, K, delta, alpha, betai, beta_x, pars, xmat, omega, seed = NULL){
  
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
  
  return(list(ymat = ymat, state = state_lab, U = u, Omega = omega, Jt = Jt, Tn = Tn, Nt = Nt, runs = runs, M = M, seed = seed))
  
}


#############################
# hazards: qui si costruisce p.array
haz <- function(xmat, beta.array, K, M){
  
  ld <- rep(M, K)
  x <- as.matrix(xmat)
  n <- nrow(x)
  d <- rep(1:K, ld)
  p.array <- array(0, dim = c(n-1, sum(ld))) # modifico con n perchè x è con le differenze prime, altrimenti va fino ad n-1
  
  for(t in 1:(n-1)){ # modifico con n perchè x è con le differenze prime, altrimenti va fino ad n-1
    for(k in 1:K){
      pos <- which(d == k)
      p.array[t, pos] <- p_k(1:ld[k], alpha = beta.array[k,1], betai = beta.array[k,2], beta_x = beta.array[k,-c(1,2)], x_t = x[t,])
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

# Augmented transition matrix for the HMM representation of the HSMM with time-varying Omega
Gamma.f_inhomo <- function(Omega, p.array, K, M){
  
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
        Gamma[t,min(pos):max(pos),c(1,which(diff(d)==1)+1)[i]] <- p[t,pos]*Omega[t,k,i]
      }
    }
    
  }
  return(Gamma)
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




# EM initialization -----------------------------------------------------
EM.start <- function(y, x = NULL, K, M, model_type, dim.theta, help_mu = T, theta_in = NULL, seed = sample(1:1000, 1)){
  
  set.seed(seed)
  n <- nrow(as.matrix(y)) # n. obs
  ld <- rep(M, K)        # n. states
  if(is.null(x)){
    x <- as.matrix(rep(1, n))
    q <- 0
  }else{
    x <- as.matrix(x) # covariate
    q <- ncol(x) # number of covariates
  }
  
  # Initialize beta array
  beta.array <- array(0, dim=c(K, q+2)) # n. states times n. covariates (+2 because of the intercept and the time-effect) 
  # 
  beta.array[,1] <- runif(K, -5, 0)
  beta.array[,2] <- runif(K, 0.01, 0.4)
  if(q>1) beta.array[,3:(q+2)] <- runif(K*q, -.5, 0.5)
  
  # Initialize theta array
  theta.array <- array(dim = c(K, dim.theta))
  
  
  # set the function to compute the Gamma matrix
  if(K == 2){
    compute_gamma <- Gamma.f
  }else{
    if(model_type %in% c("simple", "pars")){
      compute_gamma <- Gamma.f
    }else{
      compute_gamma <- Gamma.f_inhomo
    }
  }
  
  # Initialize design matrix for the estimation of the dwell-time
  d <- rep(1:K, ld)
  d_matrix <- array(0, dim = c(sum(ld)*(n-1), q+2))
  d_matrix[,1] <- rep(d, n-1)
  d_matrix[,2] <- rep(unlist(lapply(1:K, function(k) (1:ld[k])+0.5)), n-1)
  for(j in 1:q) d_matrix[,2+j] <- rep(x[1:(n-1),j], each = sum(ld)) # covariata osservata, tranne l'ultimo termine
  # for(j in 1:q) d_matrix[,2+j] <- rep(c(diff(x[,j])[1], diff(x[,j])[-length(diff(x[,j]))]), each = sum(ld)) # differenze prime della covariata
  d_matrix <- as.data.frame(d_matrix)
  colnames(d_matrix) <- c("state", "time", paste0("V", 1:ncol(x)))
  
  if(model_type %in% c("simple", "pars") | K == 2){
    # HMM matrix
    Omega <- matrix(runif(K*K), K, K)     
    diag(Omega) <- 0
    # Normalize omega by row
    for(k in 1:K) Omega[k,-k]<- Omega[k,-k]/sum(Omega[k,-k]) 
  }else{
    # HMM matrix
    Omega <- array(0, dim = c(n-1, K, K))     
    # Normalize omega by row
    for(t in 1:(n-1)) {
      Omega[t,,] <- runif(K*K)
      diag(Omega[t,,]) <- 0
      for(k in 1:K) Omega[t, k,-k] <- Omega[t, k,-k]/sum(Omega[t, k,-k]) 
    } 
  }
  # HSMM matrix
  p.array <- haz(x, beta.array, K, M)
  # p.array <- haz(diff(x), beta.array, K, M)
  Gamma_mat <- compute_gamma(Omega, p.array, K, M)
  
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
              p.array = p.array,
              seed = seed))
  
}






fit_hsmm_torus <- function(y, x = NULL, init = NULL, 
                           link = "cloglog", model_type = "full",
                           max.iter = 50, tol = 10^-4, verbose = TRUE){
  
  # Checks
  if(is.null(init)) stop("Initial values for the EM algorithm must be provided!")
  
  
  # Get initial values for the EM
  # Data
  y <- as.matrix(y) # outcome
  n <- nrow(y)
  if(is.null(x)){
    x <- as.matrix(rep(1, n))
    q <- 0
  }else{
    x <- as.matrix(x) # covariate
    q <- ncol(x) # number of covariates
  }
  
  # Initialization
  # Model
  ld <- init$ld # maximum lenght of each dwell-time (k = 1,...,K)
  K <- length(ld) # number of latent classes
  M <- max(ld) # largest dwell time overall
  if(M > 1 & model_type %in% c("simple", "inhomo")) stop("for the selected model_type, M must be equal to 1")
  d <- rep(1:K,ld) # auxiliary 
  fit <- array(dim = c(n, sum(ld))) # density
  llk <- c() # likelihood
  old.llk <- -Inf # initial value of the loglik
  dif <- Inf # inizial diff
  step <- 1 # initial step
  
  # EM
  Gamma_mat <- init$Gamma # transition probability matrix: moving from h to k at time t after having dwelled m
  Omega <- init$Omega # conditional transition probability matrix of the hsmm
  d_matrix <- init$d_matrix # design matrix for the dwell times
  theta.array <- init$theta.array # initial guesses of the density
  beta.array <- init$beta.array # initial guesses of the density
  Pi <- init$Pi
  
  # set the function to compute the Gamma matrix
  if(K == 2){
    compute_gamma <- Gamma.f
  }else{
    if(model_type %in% c("simple", "pars")){
      compute_gamma <- Gamma.f
    }else{
      compute_gamma <- Gamma.f_inhomo
    }
  }
  
  
  #### EM ####
  while(dif > tol & step <= max.iter){
    
    ################
    #### E-STEP ####
    ################
    
    # Compute the density with the actual parameters
    for(k in 1:K){
      pos <- which(d==k)
      fit[,pos] <- density(y[,1], y[,2], theta.array[k,]) # here the density is a bivariate wrapped Cachy
    }
    
    # Backward forward algorithm
    bw <- backward.forward(fit = fit, Gamma = Gamma_mat, Pi = Pi, K, M)
    llk[step] <- bw$l # save the current value of the likelihood
    post.pi.aug <- bw$post.pi.aug # univariate posterior probabilities
    post.bi.pi.aug <- bw$post.bi.pi.aug # bivariate posterior probabilities
    # Normalization of post.bi.pi.aug for numerical issues
    # for(i in 1:(n-1)) post.bi.pi.aug[i,,] <- post.bi.pi.aug[i,,]/sum(post.bi.pi.aug[i,,])
    
    # Aggregate over m
    post.pi <- t(sapply(1:n, function(i){sapply(split(post.pi.aug[i,], d), sum)}))
    
    
    ################
    #### M-STEP ####
    ################
    
    # initial (augmented) probabilities
    Pi <- post.pi[1,]
    
    
    # Here the code separates accounting for the different model specifications. According to the number of latent classes (K) the code simplifies.
    # Possible model specifications are:
    # - full:   omega_hkt (time-varying) and gompertz-type dwell time estimated parametrically allowing for covariates (HSMM)
    # - pars:   omega_hk (homogeneous) and gompertz-type dwell time estimated parametrically allowing for covariates (HSMM)
    # - inhomo: omega_hkt (time-varying) and geometric dwell time (M = 1) estimated parametrically allowing for covariates
    # - simple: omega_hk (homogeneous) and geometric dwell time (M = 1) estimated parametrically allowing for covariates
    
    # !!!In the current version, the same set of covariates is used to estimate the time-varying omegas!!!
    
    # First: maximization of the Omegas
    alpha_coefs <- list() # initialize a list to save the coefficients of the multinomial model. if K = 2, the list is empty.ci sta 
    if(K == 2){
      # Trivial solution
      Omega <- matrix(c(0,1,1,0), nrow = K, byrow = T)
    }else if(K > 2){
      # Create the data to estimate a multinomial model
      pi_kht <- array(NA, dim = c(K, n-1, K-1))
      cnt <- 1
      for(k in 1:K){
        for(t in 1:(n-1)){
          cnt2 <- 1
          for(h in 1:K){
            if(k != h){
              posk <- which(d==k)
              posh <- which(d==h)
              pi_kht[cnt, t, cnt2] <- sum(post.bi.pi.aug[t, posk, posh])
              cnt2 <- cnt2+1
            }
          }
        }
        
        # Create the design matrix to estimate the omegas
        desing_mat <- as.data.frame(pi_kht[k,,])
        colnames(desing_mat) <- c(1:K)[-k]
        
        # formula multinomial
        formula_omega <- as.formula("pi ~ 1") # ok for model_type %in% c("pars", "simple")
        
        if(model_type %in% c("full", "inhomo")){
          desing_mat <- as.data.frame(cbind(pi_kht[k,,], x[1:(n-1)]))
          # for(j in 1:q) desing_mat[,ncol(desing_mat)+j] <- c(diff(x[,j])[1], diff(x[,j])[-length(diff(x[,j]))])#, each = sum(ld))
          colnames(desing_mat) <- c(c(1:K)[-k], "x")
          
          formula_omega <- as.formula("pi ~ .")
        } 
        desing_mat <- desing_mat %>% pivot_longer(1:(K-1), names_to = "pi", values_to = "w") %>% mutate(pi = factor(pi)) # long format
        
        if(K == 3){
          # Fit the logistic model
          fit_omega <- suppressWarnings( glm(formula_omega, weights = w, data = desing_mat, family = "binomial") )
          # Get the predictions of the omegas
          preds <- predict(fit_omega, type = "response") # need the probs!
          preds <- cbind(preds, 1-preds) # add the baseline probabilty
        }else{
          # Fit the multinomial model
          fit_omega <- suppressWarnings( nnet::multinom(formula_omega, weights = w, data = desing_mat, trace = F) )
          # Get the predictions of the omegas
          preds <- predict(fit_omega, type = "probs") # need the probs!
        }
        
        # Save alpha coefficients of the multinomial model
        alpha_coefs[[k]] <- coef(fit_omega)
        
        # Assing estimated probabilities to Omega
        idextr <- seq(1, nrow(preds), by = K-1) # index to save: probs are repeated due to the long data format
        appids <- 1:K
        appids[k] <- 1
        appids[-k] <- (2:K)
        if(model_type %in% c("full", "inhomo")) Omega[,k,] <- cbind(0, preds[idextr, ])[, appids]
        if(model_type %in% c("pars", "simple")) Omega[k,] <- cbind(0, preds[1, , drop = F])[, appids]
        cnt <- cnt+1
      }
    }
    
    # Second: parametric estimation of the dwell-times
    # Prepare the data
    cases_fin <- noncases_fin <- c()
    for(t in 1:(n-1)){
      cases_int <- noncases_int <- c()
      for(h in 1:K){
        posh <- which(d==h)
        if(model_type %in% c("inhomo", "simple")){
          postbipi_1 <- t(as.matrix(post.bi.pi.aug[t, posh,-posh]))
          postbipi_2 <- t(as.matrix(post.bi.pi.aug[t, posh, posh]))
        }else{
          postbipi_1 <- post.bi.pi.aug[t, posh,-posh]
          postbipi_2 <- post.bi.pi.aug[t, posh, posh]
        }
        cases_int <- c(cases_int, rowSums(postbipi_1))
        noncases_int <- c(noncases_int, rowSums(postbipi_2))
      }
      cases_fin <- c(cases_fin, cases_int)
      noncases_fin <- c(noncases_fin, noncases_int)
    }
    
    # !!!Here we only allow for the interaction, namely each parameter of the regression is state-specific!!!
    # formula for the dwell-time estimation
    formula_dwell <- paste("d_matrix", colnames(d_matrix)[-c(1,2)], sep = "$")
    formula_dwell <- as.formula(paste0("cbind(cases_fin, noncases_fin) ~ factor(d_matrix$state)*d_matrix$time +",
                                       paste0("factor(d_matrix$state)*", formula_dwell, collapse =  " + ")))
    # Fit the binomial model with the predetermined link function (c-loglog as it is for now the default)
    fit_beta <- suppressWarnings( glm(formula_dwell, family = binomial(link = link)) )
    
    # Get the regression coefficients
    coefs <- coef(fit_beta)
    ids_app <- c(1:(K+1), (K+q+2):(2*K+q))
    coefs_covar <- matrix(coefs[setdiff(1:length(coefs), ids_app)], nrow = q) # output beta matrix
    for(k in 2:K) coefs[k] <- coefs[1] + coefs[k] # first the intercepts
    for(j in (K+q+2):(2*K+q)) coefs[j] <- coefs[j] + coefs[K+1] # then the dwell-time coefficients
    for(p in 1:q){ # then the covariates effect
      for(k in 2:K)
        coefs_covar[p, k] <- coefs_covar[p,1] + coefs_covar[p,k]
    }
    # Put togher
    coefs <- c(coefs[ids_app], c(t(coefs_covar)))
    beta.array[,1:ncol(beta.array)] <- coefs
    
    if(model_type %in% c("inhomo", "simple")) beta.array[,2] <- 0 # if the model is HMM, then no dwell-time parameter
    
    # Update the hazards
    p.array <- haz(x, beta.array, K, M)
    # p.array <- haz(c(diff(x)[1], diff(x)[-length(diff(x))]), beta.array, K, M)
    Gamma_mat <- compute_gamma(Omega, p.array, K, M)
    
    # Update the density (toroidal) parameters
    for(k in 1:K){
      opt <- optim(theta.array[k,], w.loglik, obs = y, weights = post.pi[, k])
      theta.array[k,] <- c(opt$par[1], opt$par[2], (tanh(opt$par[3])+1)/2, (tanh(opt$par[4])+1)/2, tanh(opt$par[5]))
    }
    
    # Check convergence
    if(step > 5) dif <- abs(llk[step]-old.llk)
    old.llk <- llk[step]
    step <- step + 1
    if(verbose==TRUE) cat("iteration ", step-1,"; loglik = ", llk[step-1], "\n")
    
  }
  
  return(list(K = K, theta.array = theta.array, Pi = post.pi[1,], 
              Gamma = Gamma_mat, Omega = Omega, alpha_coefs = alpha_coefs,
              p.array = p.array, beta.array = beta.array,
              loglik = llk, iteration = step, post.pi = post.pi, 
              post.pi.aug = post.pi.aug, post.bi.pi.aug = post.bi.pi.aug))
  
}



compute_icl <- function(data, model_obj, model_type){
  
  post_pi <- model_obj$post.pi
  tau <- model_obj$theta.array
  K <- model_obj$K
  iterconv <- model_obj$iteration
  ll <- model_obj$loglik[length(model_obj$loglik)]
  M <- ncol(model_obj$post.pi.aug)/K
  
  if(model_type == "full") numpar <-   ncol(tau)*K + (K-1) + K*(K-2)*2 +3*K
  if(model_type == "pars") numpar <-   ncol(tau)*K + (K-1) + K*(K-1)-K +3*K
  if(model_type == "inhomo") numpar <- ncol(tau)*K + (K-1) + K*(K-2)*2 +2*K
  if(model_type == "simple") numpar <- ncol(tau)*K + (K-1) + K*(K-1)-K +2*K
  
  size <- 0.5*numpar*log(nrow(data))
  
  ICLapp <- c()
  for(k in 1:K){ ICLapp[k] <- w.loglik(theta = tau[k,], obs = data, weights = post_pi[,k]) }

  comp.like <- sum(ICLapp)
  entropy <- -sum(post_pi*log(post_pi+10e-20))
  entropy_aug <- -sum(model_obj$post.pi.aug*log(model_obj$post.pi.aug+10e-20))
  if(model_type %in% c("simple", "inhomo")){
    M <- 75
    post.pi.aug <- array(rep(model_obj$post.pi/M, M),dim=c(nrow(post_pi), K, M))
    entropy_aug <- -sum(post.pi.aug*log(post.pi.aug+10e-20))
  }
  entropy_old <- -sum(model_obj$post.bi.pi.aug*log(model_obj$Gamma+10e-20))
  
  BIC <- -2*ll + 2*size
  # ICL <- comp.like + size + entropy_aug
  ICL <- BIC + entropy_aug
  
  
  return(list(loglik = ll, numpar = numpar, size = size, icl = ICL, bic = BIC, complik = comp.like, 
              entropy = entropy, entropy_aug = entropy_aug, 
              iter = iterconv, K = K, model_type = model_type))
  
}



#############################
p.to.d <- function(p, max.range){
  dlist <- c()
  #for(k in 1:length(p)){
  dlist[1] <- p[1]
  for(i in 2:length(p)){
    dlist[i] <- p[i]*prod(1-p[1:(i-1)])
  }
  for(i in (length(p)+1):max.range){
    dlist[i] <- dlist[length(p)]*(1-p[length(p)])^(i-length(p))    
  }
  #}
  return(dlist)
}

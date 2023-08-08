#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Packages
require(tidyverse)
require(magrittr)
require(foreach)
require(doParallel)

set.seed(7777)
# Auxiliary ---------------------------------------------------------------
source("auxiliary_hsmm.R")

# Simulation study --------------------------------------------------------
K <- as.numeric(args[1]) # n. of latent regimes
ntimes <- as.numeric(args[2]) # n. of obs.
nsimul <- 250 # n. of simulated artificial datasets

# Initial probabilities
delta <- rep(1/K, K)

# Covariates
ncov <- 1
xmat <- matrix(rnorm(ntimes*ncov, sd = 3), ncol = ncov)

# Hazard true parameters
betat_0K <- list(c(-8, -3), c(-8, -5, -3), c(-8, -6, -4, -2))[[K]]
betat_1K <- list(c(0.35, 0.075), c(0.4, 0.15, 0.05), c(0.4, 0.3, 0.05, 0.15))[[K]]
betax_K <- list(c(-.5, .5), c(-.5, .2, .7), c(-.5, .2, .7, -.1))[[K]]
beta.array <- cbind(betat_0K, betat_1K, betax_K) 
# Density true parameters
taus <- list(c(0.5, 0.5, 0.2, 0.3, 0.6, 2, 2, 0.2, 0.8, 0.1),
             c(0.5, 0.5, 0.2, 0.3, 0.6, 2, 2, 0.2, 0.8, 0.1, 2, -2, 0.5, 0.5, -0.6),
             c(0.5, 0.5, 0.2, 0.3, 0.6, 2, 2, 0.2, 0.8, 0.1, -2, -2, 0.7, 0.9, -0.3, 2, -2, 0.5, 0.5, -0.6))

beta.array <- cbind(betat_0K, betat_1K, betax_K)

tau0 <- matrix(taus[[K-1]], nrow = K, byrow = T)

# Simulation --------------------------------------------------------------
numcores <- detectCores() - 2
cl <- makeCluster(numcores)
registerDoParallel(cl)

varying_M <- as.logical(as.numeric(args[3])) # set M for the estimation
if(varying_M){
  Mperc <- as.numeric(args[4]) 
}else{
  Mperc <- 1
}

#### RANDOM STARTS ####
rnd <- 100
seeds <- sample(123:100000, nsimul, replace = F)

tinit <- Sys.time()
out <- foreach(nsim = 1:nsimul, .export = ls(), .inorder = T) %dopar% {
  
  # Simulation
  sims <- simulate_hsmm_torus(n = ntimes, K = K, delta = delta, alpha = betat_0K, betai = betat_1K, 
                              beta_x = as.matrix(betax_K), pars = tau0, xmat = xmat, seed = seeds[nsim])
  # set M
  M <- round(sims$M*Mperc)
  
  # Random initialization
  loglik_rnd <- c()
  betas_rnd <- theta_rnd <- gamma_rnd <- list()
  for(j in 1:rnd){
    skip_to_next <- FALSE
    init <- EM.start(y = sims$ymat, x = xmat, K = K, M = M, help_mu = T, theta_in = tau0[,1:2], dim.theta = 5)
    tryCatch(
      hsmm2 <- cloglog.hsmm.EM(y = sims$ymat, x = xmat, interact = T, init = init, max.iter = 15, tol = 10^-4, verbose = F, link = "cloglog"),
      error = function(e) { skip_to_next <<- TRUE}
    )
    
    if(skip_to_next){ 
      next 
    }else{
      loglik_rnd[j] <- hsmm2$loglik[length(hsmm2$loglik)]
      theta_rnd[[j]] <- hsmm2$theta.array
      betas_rnd[[j]] <- hsmm2$beta.array
      gamma_rnd[[j]] <- hsmm2$Gamma
    }     
  }
  
  # Get the best solution
  id_opt <- which.max(loglik_rnd)[1]
  init <- EM.start(y = sims$ymat, x = xmat, K = K, M = M, help_mu = T, theta_in = tau0[,1:2], dim.theta = 5)
  init$beta.array <- betas_rnd[[id_opt]]
  init$theta.array <- theta_rnd[[id_opt]]
  init$Gamma <- gamma_rnd[[id_opt]]
  
  # Run the model until convergence
  tryCatch(
    hsmm2new <- cloglog.hsmm.EM(y = sims$ymat, x = xmat, interact = T, init = init, max.iter = 150, tol = 10^-6, verbose = F, link = "cloglog"),
    error = function(e) {hsmm2new = "Error"; return(hsmm2new)}
  )
  
  if(!is.character(hsmm2new)) hsmm2new <- hsmm2new[-c(4, 5, 12, 13)]
  return(hsmm2new)
}
tend <- Sys.time() - tinit
stopCluster(cl)
gc()

# Save output
outfile <- paste0("WS/SimExtCov_", K, "_", ntimes, "M_", percname, ".RData")
save.image(file = outfile)

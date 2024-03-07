#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Packages
require(tidyverse)
require(magrittr)
require(foreach)
require(doParallel)

# Auxiliary ---------------------------------------------------------------
source("hsmm_auxiliary.R")

# Data --------------------------------------------------------------------
load("dat.RData")

model_type <- args[1] # type of model
M <- as.numeric(args[2]) # max. dwell time approximation
K <- as.numeric(args[3]) # number of latent states

load(paste0("run_hsmm_", K, "_", model_type, ".RData"))

# Preparation --------------------------------------------------------------
B <- as.numeric(args[4]) # N. of bootstrap sample

# Get MLE
delta <- fit_hsmm$Pi
betat_0K <- fit_hsmm$beta.array[,1,drop=T]
betat_1K <- fit_hsmm$beta.array[,2,drop=T]
betax_K <- fit_hsmm$beta.array[,3,drop=T]
tau0 <- fit_hsmm$theta.array
X <- as.matrix(dat[, c("Wind.Vv")])
omega_sim <- fit_hsmm$Omega

# Bootstrap --------------------------------------------------------------
cl <- makeCluster(20)
registerDoParallel(cl)

pars_boot <- foreach(b = 1:B, .combine = "rbind", .inorder = F, .export = as.vector(lsf.str()), .packages = c("tidyverse", "magrittr")) %dopar% {
  
  # Simulate
  sims <- simulate_hsmm_torus(n = nrow(dat), K = K, delta = delta, alpha = betat_0K, betai = betat_1K,
                              beta_x = as.matrix(betax_K), pars = tau0, xmat = X, omega = omega_sim, seed = NULL)
  
  Y <- sims$ymat
  init <- EM.start(y = Y, x = X, K = K, M = M, help_mu = F, theta_in = NULL, dim.theta = 5, model_type = model_type)
  
  init$beta.array <- fit_hsmm$beta.array
  init$theta.array <- fit_hsmm$theta.array
  if(M > 1){
    init$Gamma <- fit_hsmm$Gamma
    init$Omega <- fit_hsmm$Omega
  }
  
  hsmm2 = "Empty"
  tryCatch(
    hsmm2 <- fit_hsmm_torus(y = Y, x = X, init = init, max.iter = 100, tol = 10^-4, verbose = T, link = "cloglog", model_type = model_type),
    error = function(e) {hsmm2 = "Error"}
  )
  if(!is.character(hsmm2)){
    return(cbind(hsmm2$beta.array, hsmm2$theta.array, hsmm2$Omega, b, 1:K, sims$seed, hsmm2$loglik[length(hsmm2$loglik)]))
  }else{
    return(rep(NA, 16))
  }
  
}
stopCluster(cl)

save(pars_boot, file = paste0("Bootstrap", args[5], "_", K, "_", model_type, ".RData"))

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Packages
require(tidyverse)
require(magrittr)
require(foreach)
require(doParallel)

# Auxiliary ---------------------------------------------------------------
source("auxiliary_hsmm.R")

# Data --------------------------------------------------------------------
load("dat.RData")
load("auxiliary_obj_boot.RData")

Y <- dat[, c("Wind.Dir", "Wave.Dir")]
# X <- scale(dat[, c("Wind.Vv", "Wave.Hs")])
X <- dat[, c("Wind.Vv")]

# Estimation --------------------------------------------------------------
K <- as.numeric(args[1])
M <- as.numeric(args[2])
B <- as.numeric(args[3]) # N. of bootstrap sample

cl <- makeCluster(100)
registerDoParallel(cl)

pars_boot <- foreach(b = 1:B, .combine = "rbind", .inorder = F, .export = as.vector(lsf.str()), .packages = c("tidyverse", "magrittr")) %dopar% {
                       
                       init <- EM.start(y = Y, x = X, K = K, M = M, help_mu = F, theta_in = NULL, dim.theta = 5)
                       
                       init$beta.array <- beta_init[[1]][[ids_opt[1]]]
                       init$theta.array <- theta_init[[1]][[ids_opt[1]]]
                      
                       hsmm2 = "Empty"
                       tryCatch(
                         hsmm2 <- cloglog.hsmm.EM(y = Y, x = X, interact = T, init = init, max.iter = 150, tol = 10^-6, verbose = F, link = "cloglog"),
                         error = function(e) {hsmm2 = "Error"}
                       )
                       if(!is.character(hsmm2)){
                         return(cbind(hsmm2$beta.array, hsmm2$theta.array, b, 1:K))
                       }else{
                         return(rep(NA, 8))
                       }
                       
                     }
stopCluster(cl)

save(pars_boot, file = paste0("WS/Bootstrap", args[4], "_K", K, ".RData"))

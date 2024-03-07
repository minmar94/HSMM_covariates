#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Packages
require(tidyverse)
require(magrittr)


# Auxiliary ---------------------------------------------------------------
source("hsmm_auxiliary.R")

# Data --------------------------------------------------------------------
load("dat.RData")

Y <- dat[, c("Wind.Dir", "Wave.Dir")]
X <- dat[, c("Wind.Vv")]

dat %>% 
  ggplot(aes(Wind.Dir, Wave.Dir, color = Wind.Vv)) +
  geom_point(size = 2) +
  labs(x = "Wind direction [rad]", y = "Wave direction [rad]", color = "Wind speed (m/s)") +
  scale_x_continuous(breaks = seq(-3,3, by = 1.5), labels = c("S", "W", "N", "E", "S")) +
  scale_y_continuous(breaks = seq(-3,3, by = 1.5), labels = c("S", "W", "N", "E", "S")) +
  scale_color_distiller() +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 16))



# Estimation --------------------------------------------------------------
# Data preparation
model_type <- args[1] # type of model
M <- as.numeric(args[2]) # max. dwell time approximation
K <- as.numeric(args[3]) # number of latent states

rnd_init <- 100
theta_in <- matrix(c(-2, 2, 0, 2, 0, 0, -2, 0, -1, 1, 2, -1, -1, 0), byrow = T, ncol = 2)

# EM initialization
liks_init <- -Inf
set.seed(1234)
# EM short-run strategy
for(i in 1:rnd_init){
  
  # Initialization
  init <- EM.start(y = Y, x = X, K = K, M = M, model_type = model_type, help_mu = T, theta_in = theta_in[1:K,], dim.theta = 5)
  fit_init <- fit_hsmm_torus(Y, X, init, model_type = model_type, max.iter = 10, link = "cloglog")
  best_lik <- fit_init$loglik[length(fit_init$loglik)]
  
  # Save the best short-run
  if(best_lik > liks_init){
    liks_init <- best_lik
    init_opt <- init
  }
  
}

# Estimate till convergence
fit_hsmm <- fit_hsmm_torus(Y, X, init_opt, model_type = model_type, link = "cloglog", max.iter = 150, tol = 10e-04, verbose = T)

model_sel <- compute_icl(Y, fit_hsmm, model_type)


save(fit_hsmm, model_sel, file = paste0("JRSSC_Review/GithubRepo/run_hsmm_", K, "_", model_type, ".RData"))



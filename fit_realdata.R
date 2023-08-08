
# Packages
require(tidyverse)
require(magrittr)
require(foreach)
require(doParallel)

# Auxiliary ---------------------------------------------------------------
source("auxiliary_hsmm.R")

# Data --------------------------------------------------------------------
load("dat.RData")

Y <- dat[, c("Wind.Dir", "Wave.Dir")]
# X <- scale(dat[, c("Wind.Vv")])
X <- dat[, c("Wind.Vv")]

dat %>% 
  ggplot(aes(Wind.Dir, Wave.Dir, color = Wind.Vv)) +
  geom_point(size = 2) +
  labs(x = "Wind direction [rad]", y = "Wave direction [rad]", color = "Wind speed (m/s)") +
  scale_x_continuous(breaks = seq(-3,3, by = 1.5), labels = c("S", "W", "N", "E", "S")) +
  scale_y_continuous(breaks = seq(-3,3, by = 1.5), labels = c("S", "W", "N", "E", "S")) +
  theme_bw() +
  theme(legend.position = "top", text = element_text(size = 16))


# Estimation --------------------------------------------------------------
# Model with no mixture
# optim(abs(rnorm(5)), w.loglik, obs = Y)$value

K <- 2:4
theta_in <- matrix(c(-2, 2, 0, 2, 0, 0, -2, 0, -1, 1, 2, -1), byrow = T, ncol = 2)
rndstart <- 100

beta_init <- theta_init <- list()
llmat <- matrix(NA, nrow = length(K), ncol = rndstart)
hsmm2 <- NULL

# Optimal initialization with multiple random starts
for(k in 1:length(K)){
  
  beta_init[[k]] <- theta_init[[k]] <- list()
  
  for(i in 1:rndstart){
    
    init <- EM.start(y = Y, x = X, K = K[k], M = 10, help_mu = T, theta_in = theta_in[1:K[k],], dim.theta = 5)
    
    tryCatch(
      hsmm2 <- cloglog.hsmm.EM(y = Y, x = X, interact = T, init = init, max.iter = 15, tol = 10^-4, verbose = TRUE, link = "cloglog"),
      error = function(e){print("Error")}
    )
    
    if(!is.null(hsmm2)){
      llmat[k, i] <- hsmm2$loglik[length(hsmm2$loglik)]
      beta_init[[k]][[i]] <- hsmm2$beta.array
      theta_init[[k]][[i]] <- hsmm2$theta.array
    }
    hsmm2 = NULL
  }
  # print(k)
}

# Get best solution
ids_opt <- apply(llmat, 1, which.max)
hsmm_all <- list()
Ms <- c(25, 50, 75)

par_mk <- expand_grid(K = K,M =  Ms) %>% as.data.frame() %>% arrange(K)

for(k in 1:nrow(par_mk)){
  
  init <- EM.start(y = Y, x = X, K = par_mk$K[k], M = par_mk$M[k], help_mu = F, theta_in = NULL, dim.theta = 5)
  
  init$beta.array <- beta_init[[par_mk$K[k]-1]][[ids_opt[par_mk$K[k]-1]]]
  init$theta.array <- theta_init[[par_mk$K[k]-1]][[ids_opt[par_mk$K[k]-1]]]
  
  tryCatch(
    hsmm_all[[k]] <- cloglog.hsmm.EM(y = Y, x = X, interact = T, init = init, max.iter = 150, tol = 10^-6, verbose = TRUE, link = "cloglog")
  )
  # print(k)
  gc()
}

saveRDS(hsmm_all, file = "run_realdata.rds")

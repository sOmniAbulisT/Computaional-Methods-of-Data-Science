rm(list = ls())

#--- parameters ---#
pA <- 0.8; pB <- 0.4
alpha <- 0.1; beta <- 0.08

flitering <- function(piA, piB, y, pA, pB, alpha, beta){
  # predict
  predA <- piA*(1 - alpha) + piB*beta
  predB <- piA*alpha + piB*(1 - beta)
  
  likA <- if (y == 0) pA else (1 - pA)
  likB <- if (y == 0) pB else (1 - pB)
  # update 
  uA <- predA * likA; uB <- predB * likB
  Z <- uA + uB
  c(piA = uA/Z, piB = uB/Z)
}

filter_step(piA, piB, 0, pA, pB, alpha, beta)
filter_step(piA, piB, 1, pA, pB, alpha, beta)

###----------------------------------------------###
###       Rare-event Important Sampling          ###
###  Estimate P(X >= 200) for X ~ Poission(0.95) ###
###----------------------------------------------###
lambda1 <- c(100, 150, 200, 250)

result <- data.frame(lambda1 = numeric(),
                     Estimate = numeric(),
                     SE = numeric(),
                     CI_L = numeric(),
                     CI_U = numeric(),
                     ESS_ratio = numeric())

for(i in lambda1){
  #--- Generating samples from q(x) ---#
  x <- rpois(100000, i)
  
  #--- Computing w(x) ---#
  w <- exp(dpois(x, 0.95, log = TRUE)-dpois(x, i, log = TRUE))
  
  #--- Indicater function ---#
  I <- as.numeric(x >= 200)
  
  #--- Important sampling estimator ---#
  theta_hat <- mean(w*I)
  
  #--- MC SE and 95% CI ---#
  SE <- sd(w*I)/sqrt(100000)
  CI <- theta_hat+c(-1, 1)*qnorm(0.975)*SE
  
  #--- Effective sample size ---#
  ESS <- (sum(w)^2)/sum(w^2)
  ratio <- ESS/100000
  
  result <- rbind(result, 
                  data.frame(lambda1 = i, 
                             Estimate=theta_hat, 
                             probaility=log10(theta_hat), 
                             SE=SE, 
                             CI_L=CI[1], 
                             CI_U=CI[2], 
                             ESS_ratio=ratio))
}

###############################################
####----- ICA (Kurtosis Maximization) -----####
###############################################

ICA <- function(Z, ncomp = nrow(Z), tol = 1e-12, iter = 1000){
  
  m <- nrow(Z); n <- ncol(Z)
  k <- min(ncomp, m)
  W <- matrix(0, m, k)
  iters <- integer(k)
  
  #--- iteration ---#
  for(i in seq_len(k)){
    #--- initial
    w <- rnorm(m); w <- w / sqrt(sum(w^2))
    
    for(t in seq_len(iter)){
      u  <- as.numeric(t(w) %*% Z)      
      g  <- u^3
      gp <- 3 * u^2
      
      wk <- (Z %*% (g / n)) - mean(gp) * w
      if(i > 1){
        Proj <- W[, 1:(i - 1), drop = FALSE]
        wk <- wk - Proj %*% (t(Proj) %*% wk)
      }
      wk <- wk / sqrt(sum(wk^2))
      if(abs(1 - abs(drop(t(wk) %*% w))) < tol){
        w <- wk
        iters[i] <- t
        break
      }
      
      w <- wk
      iters[i] <- t
    }
    
    W[, i] <- w   
  }
  
  S <- t(W) %*% Z      
  
  return(list(S = S))
}

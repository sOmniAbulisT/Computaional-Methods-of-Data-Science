##################################################################
####----- This function is EigenValues and EigenProblems -----####
##################################################################

# PowerIter <- function(X, x0 = NULL, iter = 2000, tol = 1e-12){
#   n <- nrow(X)
#   if(is.null(x0)) x0 <- rep(1, n)
#   x <- x0 / max(abs(x0))
#   
#   lambda_prev <- NA
#   k <- 0
#   repeat{
#     k <- k + 1
#     y <- X %*% x
#     xk <- as.numeric(y) / max(abs(y))
#     lambda <- as.numeric(t(xk) %*% X %*% xk) / as.numeric(t(xk) %*% xk)
#     
#     if(!is.na(lambda_prev) && abs(lambda - lambda_prev) <= tol) break
#     if(max(abs(xk - x)) <= tol) break
#     if(k >= iter) break
#     
#     x <- xk
#     lambda_prev <- lambda
#   }
#   
#   return(list(EigenValue = lambda, EigenVector = xk))
# }

###################################
####----- Power Iteration -----####
###################################

PowerIter <- function(X, x0 = NULL, iter = 2000, tol = 1e-12, k = nrow(X)){
  n <- nrow(X)
  values <- numeric(k); vectors <- matrix(NA_real_, n, k)
  
  for(j in seq_len(k)){
    x <- if(is.null(x0)) rep(1, n) else x0
    if(j > 1){
      for(i in seq_len(j-1)){
        coeff <- as.numeric(t(x)%*%vectors[, i])
        x <- x - coeff * vectors[, i]
      } 
    }
    if (max(abs(x)) == 0) x <- runif(n)
    x <- x / max(abs(x))
    
    ### Iteration
    lam_prev <- NA_real_
    for(m in seq_len(iter)){
      y <- X %*% x
      normy <- max(abs(y))
      xk <- y / normy
      lambda <- as.numeric(t(xk) %*% X %*% xk) / as.numeric(t(xk) %*% xk)
      
      if(!is.na(lam_prev) && abs(lambda - lam_prev) <= tol){x <- xk; break}
      if(max(abs(xk - x)) <= tol){x <- xk; break}
      x <- xk; lam_prev <- lambda
    }
    
    v <- x / sqrt(sum(x^2))
    lambda <- as.numeric(t(v)%*%X%*%v)
    values[j] <- lambda
    vectors[, j] <- v
    
    X <- X - lambda * (v %*% t(v))
  }
  
  return(list(EigenValue = values, EigenVector = vectors))
} |> compiler::cmpfun()

#######################################
####----- Othogonal Iteration -----####
#######################################

OthoIter <- function(X, iter = 2000, tol = 1e-12){
  Ak <- X; n <- nrow(X); Qacc <- diag(n); off_hist <- numeric(0)
  
  for (k in seq_len(iter)) {
    S <- QR(Ak); Qk <- S$Q; Rk <- S$R
    
    Ak <- Rk %*% Qk
    Qacc <- Qacc %*% Qk
    
    off <- if(n > 1) max(abs(Ak[lower.tri(Ak)])) else 0
    off_hist[k] <- off
    
    if(off < tol) break
  }
  
  eigenvalues <- c(); i <- 1
  while(i <= n){
    if(i == n || abs(Ak[i + 1, i]) < tol){
      eigenvalues <- c(eigenvalues, Ak[i, i])
      i <- i + 1
    } else {
      a <- Ak[i,i]; b <- Ak[i,i+1]; c <- Ak[i+1,i]; d <- Ak[i+1,i+1]
      tr <- a + d; det <- a*d - b * c
      disc <- tr^2 - 4 * det
      rtd <- sqrt(as.complex(disc))
      eigenvalues <- c(eigenvalues, (tr + rtd)/2, (tr - rtd)/2)
      i <- i+2
    }
  } 
  return(list(EigenValues = eigenvalues, EigenVectors = Qacc))
} |> compiler::cmpfun()

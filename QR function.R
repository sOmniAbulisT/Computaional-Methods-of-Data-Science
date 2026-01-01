#####################################################
####----- This function is QR decomposition -----####
#####################################################
# Using (Modified) Gram-Schmidt Algorithm
# X: Import an m x n arbitrary matrix
# Q: An m x n orthogonal matrix 
# R: An n x n upper triangular matrix

QR <- function(X){
  stopifnot(is.matrix(X))
  
  B <- X; m <- nrow(X); n <- ncol(X); rmax <- min(m, n)
  Q <- matrix(0, nrow = m, ncol = rmax)
  R <- matrix(0, nrow = rmax, ncol = n)
  r <- 0
  
  for(k in seq_len(n)){
    v <- B[, k]
    
    if(r > 0){
      coeff <- t(Q[, seq_len(r)]) %*% v
      v <- v - Q[, seq_len(r)] %*% coeff
      R[seq_len(r), k] <- R[seq_len(r), k] + coeff[, 1]
    }
    nk <- sqrt(sum(v^2))
    
    if(nk <= 1e-12) next
    
    r <- r + 1
    Q[, r] <- v / nk
    R[r, k] <- nk
    if (r == rmax) break
  }
  
  if (r == 0) {
    Q <- matrix(0, m, 0)
    R <- R[0, ]
  } else {
    Q <- Q[, seq_len(r)]
    R <- R[seq_len(r), ]
  }

  return(list(Q = Q, R = R))
} |> compiler::cmpfun()

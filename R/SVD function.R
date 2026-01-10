#######################
####----- SVD -----####
#######################

SVD <- function(X){ 
  m <- nrow(X); n <- ncol(X) 
  xTx <- crossprod(X) 
  eig <- PowerIter(xTx) 
  
  S <- matrix(0, m, n) 
  for(i in seq_len(m)){ 
    S[i, i] <- sqrt(pmax(eig$EigenValue[i], 0))
  } 
  
  V <- eig$EigenVector 
  
  U <- matrix(0, m, m) 
  for(j in seq_len(m)){ 
    U[, j] <- (X %*% V[, j]) / S[j, j] 
  } 
  
  return(list(U = U, S = S, V = V)) 
}


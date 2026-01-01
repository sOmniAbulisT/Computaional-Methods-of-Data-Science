

PCA <- function(X, method = c("power", "SVD")){
  
  ### Center Data
  n <- nrow(X); p <- ncol(X)
  mu <- colSums(X) / n
  Xc <- sweep(X, MARGIN = 2, STATS = mu, FUN = "-")
  
  ### Compute the Covariance Matrix
  COV <- (t(Xc) %*% Xc) / (n - 1)
  
  ### Solve Eigenvalue Problem
  Eigen <- EigenValue(COV)
  lambda <- Eigen$eigenvalues
  v <- Eigen$Q
  
  
  
} |> compiler::cmpfun()
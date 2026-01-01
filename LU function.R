#####################################################
####----- This function is LU decomposition -----####
#####################################################

# X: Import an n x n square matrix
# U: An n x n upper triangular matrix
# L: An n x n lower triangular matrix

LU <- function(X){
  #--- Basic Check ---#
  if(!is.matrix(X)) X <- as.matrix(X)
  if(nrow(X) != ncol(X)) stop("X must be square. ")
  n <- nrow(X)
  L <- diag(1, n)
  U <- matrix(0, n, n)
  
  ### Row (U)
  for (k in seq_len(n)) {
    for(j in k:n){
      K <- if(k == 1) 0 else sum(L[k, seq_len(k-1)] * U[seq_len(k-1), j])
      U[k, j] <- X[k, j] - K
    }
    
    if(k < n){
      ### Column (L)
      for(i in (k + 1):n){
        K <- if(k == 1) 0 else sum(L[i, seq_len(k-1)] * U[seq_len(k-1), k])
        L[i, k] <- (X[i, k] - K) / U[k, k]
      }
    }
  }
  
  return(list(UpperTriangular = U, LowerTriangular = L))
} |> compiler::cmpfun()

#' Lower - Upper triangular matrix Decomposition 
#' INPUT: A - a squared matrix having dimension N*N
#'       TOL - small tolerance number to detect failure when the matrix is near degenerate
#' OUTPUT: The matrix A is changed, it contains a copy of both matrices L and U as such that
#'         PA=LU.
#' The permutation matrix is not sorted as a matrix, but in an integer vector P of size 
#' N+1 containing column indexes 
#' where the permutation matrix has "l". The last element 

LUDecompose <- function(A, TOL = 1e-9){
  #--- Basic Check ---#
  n <- nrow(A)
  if(ncol(A) != n){
    stop("A matrix must be squared.")
  }
  
  P <- 0:n
  #--- LU with partial pivoting ---#
  for(i in seq_len(n)){
    
    maxA <- 0; imax <- i
    for(k in i:n){
      absA <- abs(A[k, i])
      if(absA > maxA){
        maxA <- absA
        imax <- k
      }
    }
    
    if(maxA < TOL){
      stop("Matrix is degenerate")
    }
    
    if(imax != i){
      #--- pivoting P --- #
      j <- P[i]
      P[i] <- P[imax]
      P[imax] <- j
      
      #--- pivoting rows of A ---#
      ptr <- A[i, ]
      A[i, ] <- A[imax, ]
      A[imax, ] <- ptr
      
      #--- counting pivots starting from N (for determinant) ---#
      P[n+1] <- P[n+1]+1
    }
    
      #--- Gaussian elimination ---#
      if(i+1 <= n){
        for(j in (i+1):n){
          A[j, i] <- A[j, i]/A[i, i] # L[j, i]
          
          if(i+1 <= n){
            for(k in (i+1):n){
              A[j, i] <- A[j, k]-A[j, i]*A[i, k] # U[j,i]
            }
          }
        }
      }
  }
  
   L <- diag(1, n); U <- matrix(0, n, n)
   
   for(i in seq_len(n)){
     for(j in seq_len(n)){
       if(i>j){
         L[i, j] <- A[i, j] # lower triangular matrix
       }else{
         U[i, j] <- A[i, j] # upper triangular matrix
       }
     }
   }
  #--- LU decomposition done ---#
  return(list(L=L, U=U, P=P))
}

#' Solving the inverse matrix via LU decomposition 
#' INPUT: A and P matrices filled in LUDecompose with N*N dimension
#' OUTPUT: The inverse of the initial matrix
#' 

LUInverse <- function()
#' Lower - Upper triangular matrix Decomposition 
#' INPUT: A - a squared matrix having dimension N*N
#'       TOL - small tolerance number to detect failure when the matrix is near degenerate
#' OUTPUT: The matrix A is changed, it contains a copy of both matrices L and U as such that
#'         PA=LU.
#' The permutation matrix is not sorted as a matrix, but in an integer vector P of size 
#' N+1 containing column indexes 
#' where the permutation matrix has "l". The last element 

LUDecompose <- function(A, TOL = 1e-9){
  
}

#' Solving the inverse matrix via LU decomposition 
#' INPUT: A and P matrices filled in LUDecompose with N*N dimension
#' OUTPUT: The inverse of the initial matrix
#' 

LUInverse <- function()
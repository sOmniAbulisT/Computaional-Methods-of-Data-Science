#include <Rcpp.h>
#include <cmath>
using namespace Rcpp; 

// [[Rcpp::export]]
List LUDecompose(NumericMatrix mat){
  
  int n = mat.nrow();
  if(mat.ncol() != n) stop("Matrix must be square"); 
  
  NumericMatrix A = clone(mat); 
  
  //Initial permutation matrix
  NumericMatrix P_mat(n, n); 
  for(int i = 0; i < n; i++){
    P_mat(i, i)=1;
  }
  
  for(int i = 0; i < n; i++){
    
    // pivoting
    double maxVal = 0.0;
    int imax = i; 
    
    for(int k = i; k < n; k++){
      if(std::abs(A(k, i)) > maxVal){
        maxVal = std::abs(A(k, i)); 
        imax = k; 
      }
    }
    
    if(maxVal < 1e-12) stop("Matrix is singular"); 
    
    // row swapping
    if(imax != i){
      for(int col = 0; col < n; col++){
        double temp = A(i, col); 
        A(i, col) = A(imax, col); 
        A(imax, col) = temp; 
      }
      
      for(int col = 0; col < n; col++){
        double temp = P_mat(i, col); 
        P_mat(i, col) = P_mat(imax, col);
        P_mat(imax, col) = temp; 
      }
    }
    
    // Doolittle Algorithm
    for(int j = i + 1; j < n; j++){
      // L matrix
      A(j, i) /= A(i, i); 
      
      for(int k = i + 1; k < n; k++){
        // U matrix
        A(j, k) -= A(j, i)*A(i, k);
      }
    }
  }
  
  NumericMatrix L(n, n); 
  NumericMatrix U(n, n); 
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i > j){
        // lower triangular
        L(i, j) = A(i, j); 
        U(i, j) = 0; 
      } else if(i == j) {
        // diagonal
        L(i, j) = 1; 
        U(i, j) = A(i, j); 
      } else {
        // upper triangular
        L(i, j) = 0; 
        U(i, j) = A(i, j); 
      }
    }
  }
  
  return List::create(
    Named("L") = L, 
    Named("U") = U, 
    Named("P") = P_mat
  ); 
}


#include <Rcpp.h>
#include <cmath>
using namespace Rcpp; 

// [[Rcpp::export]]

List CholeskyDecompose(NumericMatrix mat){
  int n = mat.nrow(); 
  if(mat.ncol() != n) stop("Matrix must be squared. "); 
  
  NumericMatrix L(n, n); 
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j <= i; j++){
      double sum = 0.0; 
      
      for(int k = 0; k < j; k++){
        sum += L(i, k)*L(j, k); 
      }
      
      if(i == j){
        double val = mat(i, i) - sum; 
        
        if(val <= 0) stop("Matrix is not positive definite. "); 
        L(i, j) = std::sqrt(val); 
      } else {
        L(i, j) = (mat(i, j) - sum) / L(j, j); 
      }
    }
  }
  
  return List::create(
    Named("L") = L
  ); 
}
  
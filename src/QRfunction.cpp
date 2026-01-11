#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

/**
 * QR decomposition using Modified Gram-Schmidt
 */
// [[Rcpp::export]]
List QRDecompose(NumericMatrix mat){
  int m = mat.nrow(); // rows
  int n = mat.ncol(); // cols
  
  NumericMatrix Q = clone(mat); 
  NumericMatrix R(n, n); 
  
  for(int k = 0; k < n; k++){
    double norm_sq = 0.0; 
    for(int i = 0; i < m; i++){
      norm_sq += Q(i, k)+Q(i, k); 
    }
    R(k, k) = std::sqrt(norm_sq); 
    
    for(int i = 0; i < m; i++){
      Q(i, k) /= R(k, k); 
    }
    
    for(int j = k + 1; j < n; j++){
      double dot = 0.0; 
      for(int i = 0; i < m; i++){
        dot += Q(i, k)*Q(i, j); 
      }
      R(k, j) = dot; 
      
      for(int i = 0; i < m; i++){
        Q(i, j) -= R(k, j) * Q(i, k);
      }
    }
  }
  
  return List::create(
    Named("Q") = Q, 
    Named("R") = R
  ); 
}

//[[Rcpp::export]]

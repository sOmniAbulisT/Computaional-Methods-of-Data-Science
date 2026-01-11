#include <Rcpp.h>
#include <cmath>
using namespace Rcpp; 

// [[Rcpp::export]]

List CholDecompose(NumericMatrix mat){
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

// [[Rcpp::export]]
NumericVector CholSolve(List chol_result, NumericVector b){
  NumericMatrix L = chol_result["L"]; 
  int n = L.nrow(); 
  if(b.size() != n) stop("Vector b dimension does not match matrix. "); 
  
  // ----- Forward Substitution -----
  NumericVector y(n); 
  for(int i = 0; i < n; i++){
    double sum = 0.0; 
    for(int k = 0; k < i; k++){
      sum += L(i, k)*y[k]; 
    }
    y[i] = (b[i] - sum)/L(i, i); 
  }
  
  // ----- Backward Substitution -----
  NumericVector x(n);
  for(int i = n - 1; i >= 0; i--){
    double sum = 0.0; 
    for(int k = i + 1; k < n; k++){
      sum += L(k, i)*x[k]; 
    }
    x[i] = (y[i]-sum)/L(i, i); 
  }
  
  return x; 
}

/**
 * Calculate Inverse Matrix using Cholesky Decomposition
 */
// [[Rcpp::export]]
NumericMatrix CholInvert(List chol_result){
  NumericMatrix L = chol_result["L"]; 
  int n = L.nrow(); 
  NumericMatrix IA(n, n); 
  
  for(int j = 0; j < n; j++){
    NumericVector b(n); 
    b[j] = 1.0; 
    
    NumericVector y(n);
    for(int i = 0; i < n; i++){
      double sum = 0.0; 
      for(int k = 0; k < i; k++){
        sum += L(i, k)*y[k]; 
      }
      y[i] = (b[i] - sum) / L(i, i); 
    }
    
    for(int i = n - 1; i >= 0; i--){
      double sum = 0.0; 
      for(int k = i + 1; k < n; k++){
        sum += L(k, i) * IA(k, j); 
      }
      IA(i, j) = (y[i] - sum) / L(i, i); 
    }
  }
  
  return IA; 
}

/**
 * Calculate Determinant using Cholesky Decomposition
 * * det(A) = (prod(diag(L)))^2
 * * @param chol_result List containing "L". 
 * @return Determinant of A.
 */
//[[Rcpp::export]]
double CholDeterminant(List chol_result){
  NumericMatrix L = chol_result["L"]; 
  int n = L.nrow(); 
  
  double det_L = 1.0; 
  for(int i = 0; i < n; i++){
    det_L *= L(i, i); 
  }
  
  // det(A)=det(L*L^T)=det(L)*det(L^T)=det(L)^2
  return det_L*det_L; 
}
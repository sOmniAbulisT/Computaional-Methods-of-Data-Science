#include <Rcpp.h>
#include <cmath>
using namespace Rcpp; 

// [[Rcpp::export]]
List LUDecompose(NumericMatrix mat){
  
  int n = mat.nrow();
  if(mat.ncol() != n) stop("Matrix must be square"); 
  
  NumericMatrix A = clone(mat); 
  
  // permutation vector
  IntegerVector P(n); 
  for(int i = 0; i < n; i++){
    double maxVal = 0.0; 
    int imax = i; 
    
    for(int k = i; k < n; k++){
      if(std::abs(A(k, i)) > maxVal){
        maxVal = std::abs(A(k, i)); 
      }
    }
    
    if(maxVal < 1e-12) stop("Matrix is singular"); 
    
    // row swapping
    if(imax != i){
      int tempP = P[i]; 
      P[i] = P[imax]; 
      P[imax] = tempP; 
      
      for(int col = 0; col < n; col++){
        double tempVal = A(i, col); 
        A(i, col) = A(imax, col);
        A(imax, col) = tempVal; 
      }
    }
    
    // Doolittle 
    for(int j = i+1; j < n; j++){
      // L matrix; A[j][i] = A[j][i]/A[i][i]
      A(j, i) /= A(i, i); 
      
      // U matrix; A[j][k]=A[j][k]-L[j][i]*U[i][k]
      for(int k = i+1; k<n; k++){
        A(j, k)-=A(j, i)*A(i, k)
      }
    }
  }
  
  return List::creat(Named("LU_Matrix") = A, Named("Permutation") = P);
}


#include <Rcpp.h>
#include <cmath>
using namespace Rcpp; 

/**
 * LU Decomposition with Partial Pivoting (Doolittle's Method)
 * * This function decomposes a square matrix A into P * A = L * U, where:
 * - L is a unit lower triangular matrix (diagonal elements are 1)
 * - U is an upper triangular matrix
 * - P is a permutation matrix recording row swaps
 *
 * @param mat A square NumericMatrix (NxN)
 * @return A List containing:
 * - "L": The lower triangular matrix (NxN)
 * - "U": The upper triangular matrix (NxN)
 * - "P": The permutation matrix (NxN)
 * - "swaps": Integer count of row swaps performance 
 */

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
  int swaps = 0; 
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
      swaps++; 
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
    Named("P") = P_mat, 
    Named("swaps") = swaps
  ); 
}

/**
 * Solve a System of Linear Equations (Ax = b) using LU Decomposition
 *
 * This function solves for x in the equation Ax = b.
 * Since A = P^-1 * L * U, the equation becomes L * U * x = P * b.
 * The solution is found in two steps:
 * 1. Forward Substitution: Solve L * y = P * b for y.
 * 2. Backward Substitution: Solve U * x = y for x.
 *
 * @param lu_result A List returned by LUDecompose containing L, U, and P.
 * @param b A NumericVector representing the right-hand side vector.
 * @return A NumericVector x representing the solution.
 */

// [[Rcpp::export]]
NumericVector LUSolve(List lu_result, NumericVector b){
  NumericMatrix L = lu_result["L"]; 
  NumericMatrix U = lu_result["U"]; 
  NumericMatrix P = lu_result["P"]; 
  
  int n = L.nrow(); 
  if(b.size() != n) stop("Vector b dimension does not match matrix."); 
  
  NumericVector b_perm(n); 
  
  // b_perm = P*b
  for(int i = 0; i < n; i++){
    double sum = 0; 
    for(int j = 0; j < n; j++){
      sum += P(i, j)*b[j];
    }
    b_perm[i] = sum;
  }
  
  // Forward Substitution
  NumericVector y(n);
  for(int i = 0; i < n; i++){
    double sum = 0; 
    for(int j = 0; j < i; j++){
      sum += L(i, j)*y[j]; 
    }
    y[i] = b_perm[i]-sum; 
  }
  
  // Backward Substitution 
  NumericVector x(n);
  for(int i = n - 1; i >= 0; i--){
    double sum = 0; 
    for(int j = i + 1; j < n; j++){
      sum += U(i, j)*x[j];
    }
    x[i] = (y[i]-sum)/U(i, i); 
  }
  
  return x; 
}

/**
 * Calculate the Inverse Matrix using LU Decomposition
 *
 * This function computes A^-1 by solving Ax = e_j for each column j,
 * where e_j is the j-th column of the identity matrix.
 * It essentially runs LUSolve n times.
 *
 * @param lu_result A List returned by LUDecompose containing L, U, and P.
 * @return A NumericMatrix representing the inverse of A.
 */

// [[Rcpp::export]]
NumericMatrix LUInvert(List lu_result){
  NumericMatrix L = lu_result["L"]; 
  NumericMatrix U = lu_result["U"]; 
  NumericMatrix P = lu_result["P"];
  
  int n = L.nrow(); 
  NumericMatrix IA(n, n); 
  
  for(int j = 0; j < n; j++){
    NumericVector b(n); 
    for(int i = 0; i < n; i++){
      b[i] = P(i, j);     
    }
    
    NumericVector y(n); 
    for(int i = 0; i < n; i++){
      double sum = 0; 
      for(int k = 0; k < i; k++){
        sum += L(i, k)*y[k]; 
      }
      y[i] = b[i] - sum; 
    }
    
    for(int i = n-1; i >= 0; i--){
      double sum = 0; 
      for(int k = i + 1; k < n; k++){
        sum += U(i, k)*IA(k, j);
      }
      IA(i, j) = (y[i] - sum) / U(i, i); 
    }
  }
  return IA; 
}

/**
 * Calculate the Determinant of a Matrix using LU Decomposition
 *
 * Computes det(A) using the property: det(A) = det(P)^-1 * det(L) * det(U)
 * - det(L) = 1 (unit lower triangular)
 * - det(U) = product of diagonal elements
 * - det(P) = (-1)^swaps
 *
 * @param lu_result A List returned by LUDecompose containing "U" and "swaps".
 * @return A double representing the determinant of A.
 */

// [[Rcpp::export]]
double LUDeterminant(List lu_result){
  NumericMatrix U = lu_result["U"]; 
  int swaps = lu_result["swaps"]; 
  
  // det(U)
  double det = 1.0; 
  int n = U.nrow(); 
  
  for(int i = 0; i < n; i++){
    det *= U(i, i); 
  }
  
  if(swaps % 2 != 0){
    det = -det; 
  }
  
  return det; 
}
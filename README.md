# Computational Methods of Data Science

![R](https://img.shields.io/badge/Language-R-blue)
![Status](https://img.shields.io/badge/Status-Active-green)

## Overview
This repository contains my implementations of various algorithms for **Computational Methods of Data Science**. 

The primary goal of this project is to **implement statistical algorithms from scratch** using R (and Rcpp where applicable), 
without relying on built-in high-level functions. This approach focuses on understanding the underlying mathematical principles, 
numerical linear algebra, and optimization techniques.

## File structure 
- **'R/'**: Source code of core algorithms and functions
- **'data/'**: Datasets used for testing and validation (e.g., 'CAmaxTemp.txt').

## Algorithms Implemented
The following methods have been implemented:

### 1. Matrix Decomposition & Dimensionality Reduction
Located in the 'R/' directory:
- **PCA (Principal Component Analysis):** Implementation of Eigen-decomposition and SVD-based approaches ('PCA function.R').
- **SVD (Singular Value Decomposition):** Computation of singular values and vectors ('SVD function.R').
- **ICA (Independent Component Analysis):** Implementation for blind source separation ('ICA function.R').
- **QR Decomposition:** Gram-Schmidt process or Householder reflections ('QR function.R').
- **LU Decomposition:** Matrix factorization for solving linear systems ('LU function.R').
- **Eigenvalues:** Numerical methods for computing eigenvalues ('Eigenvalue function.R').

### 2. Optimization Algorithm (Work in Progress)
- Genetic Algorithms (GA)
- Simulated Annealing (SA)

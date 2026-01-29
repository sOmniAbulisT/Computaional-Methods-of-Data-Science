#' Calculated Coefficient of Determination for Genomic Selection 
#' 
#' @description
#' This function computes the CD value to evaluate the reliability of genomic prediction.
#' It supports:
#' 1. Additive (A) and Additive + Dominance (A + D) models.
#' 2. Mulit-Environment Trails (MGE model) with GxE interaction. 
#' 3. Optimizated matrix oprations for fast computation 
#'
#' @param kinshipA Matrix. Additive kinship matrix.
#' @param kinshipD Matrix (option). Dominance kinship matrix.
#' @param train List. A list containing indices of training individuals for each environment.
#' @param varA Numeric. Additive variance component (sigma^2_A).
#' @param covAxE Numeric vector. Covariance for AxE interaction. Default is rep(10, env).
#' @param varD Numeric. Dominance variance component (sigma^2_D).
#' @param covDxE Numeric vector. Covariance for DxE interaction. Default is rep(0, env).
#' @param varE Numeric vector. Residual variance (sigma^2_E). Default is NULL.
#' @param methods Character. "CDmean(v2)" (default) or "CDmean_MET" (for multi-environment mean).
#'
#' @return Numeric. The calculated CD value.
#'

CD <- function(kinshipA, kinshipD = NULL, train, varA = 10, covAxE = 0, 
               varD = 0, covDxE = 0, varE = NULL, methods = "CDmean(v2)"){
  env <- length(train); Nc <- nrow(kinshipA)
  
  if(is.null(covAxE)) covAxE <- rep(10, env)
  if(is.null(covDxE)) covDxE <- rep(0, env)
  
  #--- Variance components of Additive MGE model ---#
  OmegaA <- matrix(varA, nrow = env, ncol = env)
  diag(OmegaA) <- diag(OmegaA) + covAxE
  
  #--- Variance components of Dominance MGE model (option) ---#
  if(!is.null(kinshipD)){
    OmegaD <- matrix(varD, nrow = env, ncol = env)
    diag(OmegaD) <- diag(OmegaD) + covDxE
  }else{
    OmegaD <- matrix(varD, nrow = env, ncol = env)
    diag(OmegaD) <- diag(OmegaD) + covDxE
  }
  
  #--- Variance components of Omega_E ---#
  if(is.null(varE)) varE <- diag(OmegaA)
  
  n_train <- sapply(train, length)
  
  #--- Block Diagonal Matrix (Mt) ---#
  
  M_block <- lapply(seq_len(env), function(i){
    ni <- n_train[i]
    
    if(ni > 0){
      return((diag(1, ni)-matrix(1/ni, nrow = ni, ncol = ni))*(1/varE[i]))
    }else{
      return(matrix(0, ni, ni))
    }
  }) 
  
  Mt <- Matrix::bdiag(M_block) |> as.matrix()
  
  #--- Initialize Gt, Gct and Mt ---#
  Gt <- matrix(0, nrow = sum(n_train), ncol = sum(n_train))
  Gct <- matrix(0, nrow = Nc*env, ncol = sum(n_train))
  
  #--- Fill in the elements of the matrices Gct and Gt ---#
  end_idx <- cumsum(n_train); start_idx <- c(1, head(end_idx, -1) + 1) 
  
  for(i in seq_len(env)){
    idx_i <- start_idx[i]:end_idx[i]
    train_i <- train[[i]]
    
    for(j in seq_len(env)){
      idx_j <- start_idx[j]:end_idx[j]
      train_j <- train[[j]]
      
      # Additive
      cov_A <- OmegaA[i, j]
      Diag_Value <- cov_A * kinshipA[train_i, train_j] # train * train
      Off_Diag <- cov_A * kinshipA[, train_j] # candidate * train
      
      # Dominance
      if(!is.null(kinshipD)){
        cov_D <- OmegaD[i, j]
        Diag_Value <- Diag_Value + (cov_D * kinshipD[train_i, train_j])
        Off_Diag <- Off_Diag + (cov_D * kinshipD[, train_j])
      }
      
      # fill elements in Gt matrix
      Gt[idx_i, idx_j] <- Diag_Value
      
      # fill elements in Gct matrix
      row_start <- (i - 1)*Nc + 1; row_end <- i*Nc
      Gct[row_start:row_end, idx_j] <- Off_Diag
    }
  }
  
  #--- A matrix ---#
  MGt <- Mt%*%Gt
  diag(MGt) <- diag(MGt) + 1
  Solve_MGt <- solve(MGt, Mt)
  front <- Gct%*%Solve_MGt
  A_diag <- rowSums(front*Gct)
  
  #--- B matrix ---#
  B_diag <- rep(0, Nc*env)
  kinshipA_diag <- diag(kinshipA)
  kinshipD_diag <- if(!is.null(kinshipD)) diag(kinshipD) else rep(0, Nc)
  
  for(k in seq_len(env)){
    idx_range <- ((k-1)*Nc + 1) : (k*Nc)
    B_diag[idx_range] <- (OmegaA[k,k] * kinshipA_diag) + (OmegaD[k,k] * kinshipD_diag)
  }
  
  if(methods == "CDmean(v2)"){
    CD_value <- mean(A_diag/B_diag)
  } else if(methods == "CDmean_MET"){
    warning("CDmean_MET has not implemented, return NULL. ")
    CD_value <- NULL
  } else if(methods == "CDranking"){
    warning("CDranking has not implemented, return NULL. ")
    CD_value <- NULL  
  }
  
  return(CD_value)
}

#' Calculated the Mean Squared Predictive Error for genomic selection
#'
#'

MSPE <- function(Xt, Xc, lambda = 1){
  
}

#' Calculated the r-score for genomic selection
#'
#'

Rscore <- function(Xt, Xc){
  
}
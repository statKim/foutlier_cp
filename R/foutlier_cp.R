# devtools::install_github("statKim/mrfDepth")
library(tidyverse)
library(mrfDepth)
library(fda)
library(mvtnorm)
library(progress)  # show progress bar

# Parallel computation
library(doSNOW)
library(doRNG)
library(foreach)



#' Conformal Outlier Detection for Multivariate Functional Data
#' 
#' @param X n-m-p dimensional training data (n: # of observations, m: # of timepoints, p: # of variables)
#' @param X_test n_test-m-p dimensional test data
#' @param type the option for computing nonconformity scores. 4 optons are supported. ("depth_transform" (default), "depth", "esssup", "focsvm")
#' @param type_depth depth for multivariate functional depth. See `mrfDepth::mfd()` for available depth functions.
#' @param transform a vector containing the curve transformations for `type = "depth_transform"`. Only available on "D0", "D1" and "D2". Default is c("D0","D1","D2")
#' @param train_type if it is set "mixed", initial outlier detection is performed for input training set. "clean" (default) and "mixed" are available.
#' @param alpha the given FDR level (also be coverage level for "esssup"). Default is 0.1.
#' @param ccv CCV (calibration-conditional valid) adjustments are performed for marginal conformal p-values. It performs Simes, asymptotic adjustments. Default is TRUE.
#' @param individual conformal p-values for each individual transformation are computed. (Default is FALSE)
#' @param n_cores number of cores for parallel computing in `type = "depth_transform"`. Default is 1.
#' @param seed random seed number. Default is NULL.
#' @param ... See additional options for `split_conformal_fd()`
foutlier_cp <- function(X, X_test, 
                        type = "depth_transform", 
                        type_depth = "projdepth",
                        transform = c("D0","D1","D2"),
                        alpha = 0.1, 
                        train_type = "clean",
                        ccv = TRUE, 
                        individual = FALSE,
                        n_cores = 1,
                        seed = NULL, ...) {
  # Marginal and CCV conformal p-values
  cp_obj <- split_conformal_fd(X = X, 
                               X_test = X_test, 
                               type = type, 
                               type_depth = type_depth,
                               transform = transform,
                               alpha = alpha, 
                               train_type = train_type,
                               ccv = ccv, 
                               individual = individual,
                               n_cores = n_cores,
                               seed = seed, ...)
  
  # BH procedure for marginal and CCV p-values
  idx_bh <- apply(cp_obj$conf_pvalue, 2, BH, alpha = alpha, simplify = F)
  
  out <- list(
    cp_obj = cp_obj,
    idx_out = idx_bh
  )
  
  # BH procedure for each individual transformation
  if (isTRUE(individual) & type == "depth_transform") {
    out$idx_out_indiv <- lapply(cp_obj$conf_pvalue_indiv, function(pvalue){
      apply(pvalue, 2, BH, alpha = alpha, simplify = F)
    })
  }
  
  return(out)
}


#' Split Conformal Prediction for Multivariate Functional Data
#' 
#' @param X n-m-p dimensional training data (n: # of observations, m: # of timepoints, p: # of variables)
#' @param X_test n_test-m-p dimensional test data
#' @param type the option for computing nonconformity scores. 4 optons are supported. ("depth_transform" (default), "depth", "esssup", "focsvm")
#' @param type_depth depth for multivariate functional depth. See `mrfDepth::mfd()` for available depth functions.
#' @param transform a vector containing the curve transformations for `type = "depth_transform"`. Only available on "D0", "D1" and "D2". Default is c("D0","D1","D2")
#' @param train_type if it is set "mixed", initial outlier detection is performed for input training set. "clean" (default) and "mixed" are available.
#' @param alpha coverage level (Only used for `type = "esssup"`)
#' @param rho a proportion of the proper training set for the split conformal prediction. Default is 0.5.
#' @param n_calib a number of calibration set (If it is set, `rho` is ignored.)
#' @param weight 
#' @param ccv CCV (calibration-conditional valid) adjustments are performed for marginal conformal p-values. It performs Simes, asymptotic adjustments. Default is TRUE.
#' @param delta parmeters for CCV adjustments.
#' @param k parmeters for CCV adjustments.
#' @param individual conformal p-values for each individual transformation are computed. (Default is TRUE)
#' @param n_cores number of cores for parallel computing in `type = "depth_transform"`. Default is 1.
#' @param mfd_alpha alpha for `mrfDepth::mfd()`. See `mrfDepth::mfd()` for details.
#' @param seed random seed number. Default is NULL.
#' @param ... additional options for `FOCSVM()`
split_conformal_fd <- function(X, X_test, 
                               type = "depth_transform", 
                               type_depth = "projdepth",
                               transform = c("D0","D1","D2"),
                               train_type = "clean",
                               alpha = 0.1, 
                               rho = 0.5, n_calib = NULL,
                               weight = FALSE,
                               ccv = TRUE, delta = 0.1, k = NULL,
                               individual = TRUE,
                               n_cores = 1,
                               mfd_alpha = 0,
                               seed = NULL, ...) {
  n <- nrow(X[[1]])  # number of training data
  m <- ncol(X[[1]])  # number of timepoints
  p <- length(X)   # number of functional covariates
  n_test <- nrow(X_test[[1]])   # number of test data
  
  # Check supported transformations
  if (sum(transform %in% c("D0","D1","D2")) < length(transform)) {
    stop("Not supproted for `transform`!")
  }
  
  # Under zero-assumption, remove outliers from the mixed training set
  if (train_type == "mixed") {
    obj <- get_clean_null(X,  
                          type = type, 
                          type_depth = type_depth,
                          transform = transform,
                          weight = weight)
    idx_train_null <- obj$idx_clean_null
    X <- lapply(X, function(x){ x[idx_train_null, ] })
    n <- nrow(X[[1]])
  } else {
    idx_train_null <- NULL
  }
  
  if (is.null(n_calib)) {
    n_train <- round(n * rho)   # size of proper training set
    n_calib <- n - n_train   # size of calibration set
  } else {
    # Given number of calibration set
    n_train <- n - n_calib
  }
  
  # Depth options for `mfd`
  if (p >= n_train) {
    # High-dimensional case
    depthOptions <- list(type = "Rotation")
  } else {
    depthOptions <- NULL
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Split data
  idx_proper_train <- sample(1:n, n_train)
  idx_calib <- setdiff(1:n, idx_proper_train)
  X_train <- lapply(X, function(x){ x[idx_proper_train, ] })
  X_calib <- lapply(X, function(x){ x[idx_calib, ] })
  
  # Only use raw curves for type == "depth"
  if (type == "depth") {
    type <- "depth_transform"
    transform <- "D0"
  }
  
  if (type == "esssup") {
    # Point predictor
    pred <- lapply(X_train, function(x){ colMeans(x) })
    
    # Modulation function (t-function)
    abs_resid_train <- mapply(function(X_p, pred_p) {
      apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
    }, X_train, pred, SIMPLIFY = F)
    score_H <- sapply(abs_resid_train, function(resid_p) { 
      apply(resid_p, 2, max)   # esssup_t resid
    }) %>% 
      apply(1, max)
    gamma <- sort(score_H)[ ceiling((1 - alpha) * (n_train + 1)) ]
    idx_H <- which(score_H <= gamma)   # index of H_1
    s_ftn <- lapply(abs_resid_train, function(resid_p) {
      apply(resid_p[, idx_H], 1, max)
    }) 
    
    # Non-conformity score with modulation
    nonconform_score_calib <- mapply(function(X_p, pred_p, s_ftn_p){
      apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
    }, X_calib, pred, s_ftn) %>% 
      apply(1, max)
    idx_cutoff <- ceiling((1 - alpha) * (n_calib + 1))
    k_s <- sort(nonconform_score_calib)[idx_cutoff]
    
    # # Coverage check
    # sum(nonconform_score_calib <= k_s) / n_calib
    
    # # Conformal prediction band
    # lb <- mapply(function(pred_p, s_ftn_p){
    #   pred_p - k_s*s_ftn_p
    # }, pred, s_ftn, SIMPLIFY = F)
    # ub <- mapply(function(pred_p, s_ftn_p){
    #   pred_p + k_s*s_ftn_p
    # }, pred, s_ftn, SIMPLIFY = F)
    
    
    # Conformal p-value (marginal)
    nonconform_score_test <- mapply(function(X_p, pred_p, s_ftn_p){
      apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
    }, X_test, pred, s_ftn) %>% 
      apply(1, max)
    conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
      (1 + sum(nonconform_score_calib >= s)) / (n_calib + 1)
    })
  } else if (type == "depth_transform") {
    # Transform data structure for `mrfDepth::mfd()`
    arr_train <- array(NA, c(m, n_train, p))
    arr_calib <- array(NA, c(m, n_calib, p))
    arr_test <- array(NA, c(m, n_test, p))
    for (i in 1:p) {
      arr_train[, , i] <- t(X_train[[i]])
      arr_calib[, , i] <- t(X_calib[[i]])
      arr_test[, , i] <- t(X_test[[i]])
    }
    
    # Compute functional depth with transformations
    nonconform_score_calib_indiv <- matrix(NA, n_calib, length(transform))
    nonconform_score_test_indiv <- matrix(NA, n_test, length(transform))
    
    if (n_cores > 1) {
      # Using parallel computation
      n_cores <- min(length(transform), n_cores)
      cl <- makeCluster(n_cores)
      registerDoSNOW(cl)
      pkgs <- c("mrfDepth")
      res_cv <- foreach(s = 1:length(transform), .packages = pkgs) %dopar% {
        trans_type <- transform[s]  # transform type
        
        # Transform into 1st or 2nd derivatives
        if (trans_type == "D0") {
          # Raw curves
          arr_train_trans <- arr_train
          arr_calib_trans <- arr_calib
          arr_test_trans <- arr_test
        } else if (trans_type == "D1") {
          # 1st derivatives
          arr_train_trans <- array(NA, c(m-1, n_train, p))
          arr_calib_trans <- array(NA, c(m-1, n_calib, p))
          arr_test_trans <- array(NA, c(m-1, n_test, p))
          for (i in 1:p) {
            arr_train_trans[, , i] <- apply(arr_train[, , i], 2, diff)
            arr_calib_trans[, , i] <- apply(arr_calib[, , i], 2, diff)
            arr_test_trans[, , i] <- apply(arr_test[, , i], 2, diff)
          }
        } else if (trans_type == "D2") {
          # 2nd derivatives
          arr_train_trans <- array(NA, c(m-2, n_train, p))
          arr_calib_trans <- array(NA, c(m-2, n_calib, p))
          arr_test_trans <- array(NA, c(m-2, n_test, p))
          for (i in 1:p) {
            arr_train_trans[, , i] <- apply(arr_train[, , i], 2, function(x){ diff(diff(x)) })
            arr_calib_trans[, , i] <- apply(arr_calib[, , i], 2, function(x){ diff(diff(x)) })
            arr_test_trans[, , i] <- apply(arr_test[, , i], 2, function(x){ diff(diff(x)) })
          }
        }
        
        # Non-conformity scores using multivariate functional depths for calibration and test set
        # Lower depth is outlier => we take "-" to make nonconformity score
        depth_values_calib <- mfd(arr_train_trans, arr_calib_trans, 
                                  type = type_depth, 
                                  alpha = mfd_alpha,
                                  depthOptions = depthOptions)
        
        # Multivariate functional depth for test set
        depth_values_test <- mfd(arr_train_trans, arr_test_trans, 
                                 type = type_depth,
                                 alpha = mfd_alpha,
                                 depthOptions = depthOptions)
        
        out <- list(
          calib_score = -as.numeric(depth_values_calib$MFDdepthZ),
          test_score = -as.numeric(depth_values_test$MFDdepthZ)
        )
        
        return(out)
      }
      # End parallel backend
      stopCluster(cl)
      
      for (s in 1:length(transform)) {
        nonconform_score_calib_indiv[, s] <- res_cv[[s]]$calib_score
        nonconform_score_test_indiv[, s] <- res_cv[[s]]$test_score
      }
    } else {
      # Not use prallel computation
      for (s in 1:length(transform)) {
        trans_type <- transform[s]  # transform type
        
        # Transform into 1st or 2nd derivatives
        if (trans_type == "D0") {
          # Raw curves
          arr_train_trans <- arr_train
          arr_calib_trans <- arr_calib
          arr_test_trans <- arr_test
        } else if (trans_type == "D1") {
          # 1st derivatives
          arr_train_trans <- array(NA, c(m-1, n_train, p))
          arr_calib_trans <- array(NA, c(m-1, n_calib, p))
          arr_test_trans <- array(NA, c(m-1, n_test, p))
          for (i in 1:p) {
            arr_train_trans[, , i] <- apply(arr_train[, , i], 2, diff)
            arr_calib_trans[, , i] <- apply(arr_calib[, , i], 2, diff)
            arr_test_trans[, , i] <- apply(arr_test[, , i], 2, diff)
          }
        } else if (trans_type == "D2") {
          # 2nd derivatives
          arr_train_trans <- array(NA, c(m-2, n_train, p))
          arr_calib_trans <- array(NA, c(m-2, n_calib, p))
          arr_test_trans <- array(NA, c(m-2, n_test, p))
          for (i in 1:p) {
            arr_train_trans[, , i] <- apply(arr_train[, , i], 2, function(x){ diff(diff(x)) })
            arr_calib_trans[, , i] <- apply(arr_calib[, , i], 2, function(x){ diff(diff(x)) })
            arr_test_trans[, , i] <- apply(arr_test[, , i], 2, function(x){ diff(diff(x)) })
          }
        } else {
          stop("Not supproted for `transform`!")
        }
        
        
        # Non-conformity scores using multivariate functional depths for calibration and test set
        # Lower depth is outlier => we take "-" to make nonconformity score
        depth_values <- mfd(arr_train_trans, arr_calib_trans, 
                            type = type_depth, 
                            alpha = mfd_alpha,
                            depthOptions = depthOptions)
        nonconform_score_calib_indiv[, s] <- -as.numeric(depth_values$MFDdepthZ)
        
        # Multivariate functional depth for test set
        depth_values <- mfd(arr_train_trans, arr_test_trans, 
                            type = type_depth,
                            alpha = mfd_alpha,
                            depthOptions = depthOptions)
        nonconform_score_test_indiv[, s] <- -as.numeric(depth_values$MFDdepthZ)
      }
    }
    
    
    if (length(transform) > 1) {
      # # Scaling depths for each transformed curves
      # mean_calib <- colMeans(nonconform_score_calib)
      # sd_calib <- apply(nonconform_score_calib, 2, sd)
      # nonconform_score_calib <- t( (t(nonconform_score_calib) - mean_calib) / sd_calib )
      # nonconform_score_test <- t( (t(nonconform_score_test) - mean_calib) / sd_calib )
      
      # Weights for weighted average
      if (isTRUE(weight)) {
        weight_calib <- t( apply(nonconform_score_calib_indiv, 1, function(x){ exp(x) / sum(exp(x)) }) )
        weight_test <- t( apply(nonconform_score_test_indiv, 1, function(x){ exp(x) / sum(exp(x)) }) )
      } else {
        weight_calib <- 1/length(transform)
        weight_test <- 1/length(transform)
      }
      
      # Aggregate scores from transformations
      nonconform_score_calib <- rowSums(nonconform_score_calib_indiv * weight_calib)
      nonconform_score_test <- rowSums(nonconform_score_test_indiv * weight_test)
      # nonconform_score_calib <- apply(nonconform_score_calib_indiv, 1, max)
      # nonconform_score_test <- apply(nonconform_score_test_indiv, 1, max)
    } else {
      nonconform_score_calib <- as.numeric(nonconform_score_calib_indiv)
      nonconform_score_test <- as.numeric(nonconform_score_test_indiv)
    }
    
    # Conformal p-value (marginal)
    conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
      (1 + sum(nonconform_score_calib >= s)) / (n_calib + 1)
    })
  }
  
  out <- list(
    idx_train_null = idx_train_null,
    idx_proper_train = idx_proper_train,
    idx_calib = idx_calib,
    type = type,
    type_depth = type_depth,
    transform = transform,
    nonconform_score = list(calib = nonconform_score_calib,
                            test = nonconform_score_test),
    # weight_calib = weight_calib,
    # weight_test = weight_test,
    conf_pvalue = data.frame(marginal = conf_pvalue_marg)
  )
  
  # Calibration-conditional valid (CCV) conformal p-value
  if (isTRUE(ccv)) {
    out$conf_pvalue$simes <- ccv_conf_pvalue(out$conf_pvalue$marginal, method = "simes", 
                                             delta = delta, n_calib = n_calib, k = k)
    out$conf_pvalue$asymp <- ccv_conf_pvalue(out$conf_pvalue$marginal, method = "asymp", 
                                             delta = delta, n_calib = n_calib, k = k)
  }
  
  # Conformal p-values for each transformation
  if (isTRUE(individual) & type == "depth_transform") {
    out$nonconform_score_indiv <- list(calib = nonconform_score_calib_indiv,
                                       test = nonconform_score_test_indiv)
    out$conf_pvalue_indiv <- lapply(1:length(transform), function(i){
      # marginal p-value
      df <- data.frame(
        marginal = sapply(nonconform_score_test_indiv[, i], function(s){
          (1 + sum(nonconform_score_calib_indiv[, i] >= s)) / (n_calib + 1)
        })
      )
      
      # CCV p-value
      if (isTRUE(ccv)) {
        df$simes <- ccv_conf_pvalue(df$marginal, method = "simes", 
                                    delta = delta, n_calib = n_calib, k = k)
        df$asymp <- ccv_conf_pvalue(df$marginal, method = "asymp", 
                                    delta = delta, n_calib = n_calib, k = k)
      }
      
      return(df)
    })
  }
  
  class(out) <- "split_conformal_fd"
  
  return(out)
}



### Conformal p-value (CCV)
ccv_conf_pvalue <- function(marg_conf_pvalue, method = "simes", delta = 0.1, n_calib, k = NULL) {
  n <- n_calib
  
  if (method == "simes") {
    # Simes adjustment when k = n/2
    if (is.null(k)) {
      k <- ceiling(n/2)
    }
    
    b <- rep(1, n)
    b[1] <- 1 - delta^(1/k)
    sub <- delta
    for (i in 2:(k+1)) {
      sub <- sub * (n-k+2-i)/(n+2-i)
      b[i] <- 1 - sub^(1/k)
    }
  } else if (method == "asymp") {
    # Asymptotic adjustment
    c_n <- (-log(-log(1-delta)) + 2*log(log(n)) + 1/2*log(log(log(n))) - 1/2*log(pi)) / sqrt(2*log(log(n)))
    
    b <- sapply(1:n, function(i){
      min(1, i/n + c_n*sqrt( i*(n-i)/(n^3) ))
    })
  } else if (method == "mc") {
    
  }
  
  # Adjusted function for marginal conformal p-value
  h <- function(t) { 
    idx <- ceiling((n+1)*t)
    out <- ifelse(idx == 0, 0,
                  ifelse(idx == n+1, 1, b[idx]))
    return(out)
  }
  
  # Calibration-conditional valid p-value
  conf_pvalue_ccv <- h(marg_conf_pvalue)
  
  return(conf_pvalue_ccv)
}




#' Get clean null training indices from the mixed training set
get_clean_null <- function(X, 
                           type = "depth_transform", 
                           type_depth = "projdepth",
                           transform = c("D0","D1","D2"),
                           individual = TRUE,
                           mfd_alpha = 0,
                           n_cores = 1,
                           weight = FALSE) {
  n <- nrow(X[[1]])  # number of training data
  m <- ncol(X[[1]])  # number of timepoints
  p <- length(X)   # number of functional covariates
  
  # Check supported transformations
  if (sum(transform %in% c("D0","D1","D2")) < length(transform)) {
    stop("Not supproted for `transform`!")
  }
  
  # Depth options for `mfd`
  if (p >= n) {
    # High-dimensional case
    depthOptions <- list(type = "Rotation")
  } else {
    depthOptions <- NULL
  }
  
  # Only use raw curves for type == "depth"
  if (type == "depth") {
    type <- "depth_transform"
    transform <- "D0"
  }
  
  
  # Outlier detection using non-conformity scores (Not CP based method)
  if (type == "esssup") {
    alpha <- 0.1
    
    # Point predictor
    pred <- lapply(X, function(x){ colMeans(x) })
    
    # Modulation function (t-function)
    abs_resid <- mapply(function(X_p, pred_p) {
      apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
    }, X, pred, SIMPLIFY = F)
    score_H <- sapply(abs_resid, function(resid_p) { 
      apply(resid_p, 2, max)   # esssup_t resid
    }) %>% 
      apply(1, max)
    gamma <- sort(score_H)[ ceiling((1 - alpha) * (n + 1)) ]
    idx_H <- which(score_H <= gamma)   # index of H_1
    s_ftn <- lapply(abs_resid, function(resid_p) {
      apply(resid_p[, idx_H], 1, max)
    }) 
    
    # Non-conformity score with modulation
    nonconform_score <- mapply(function(X_p, pred_p, s_ftn_p){
      apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
    }, X, pred, s_ftn) %>% 
      apply(1, max)
  } else if (type == "depth_transform") {
    # Transform data structure for `mrfDepth::mfd()`
    arr_X <- array(NA, c(m, n, p))
    for (i in 1:p) {
      arr_X[, , i] <- t(X[[i]])
    }
    
    # Compute functional depth with transformations
    nonconform_score_indiv <- matrix(NA, n, length(transform))
    
    if (n_cores > 1) {
      # Using parallel computation
      n_cores <- min(length(transform), n_cores)
      cl <- makeCluster(n_cores)
      registerDoSNOW(cl)
      pkgs <- c("mrfDepth")
      res_cv <- foreach(s = 1:length(transform), .packages = pkgs) %dopar% {
        trans_type <- transform[s]  # transform type
        
        # Transform into 1st or 2nd derivatives
        if (trans_type == "D0") {
          # Raw curves
          arr_X_trans <- arr_X
        } else if (trans_type == "D1") {
          # 1st derivatives
          arr_X_trans <- array(NA, c(m-1, n, p))
          for (i in 1:p) {
            arr_X_trans[, , i] <- apply(arr_X[, , i], 2, diff)
          }
        } else if (trans_type == "D2") {
          # 2nd derivatives
          arr_X_trans <- array(NA, c(m-2, n, p))
          for (i in 1:p) {
            arr_X_trans[, , i] <- apply(arr_X[, , i], 2, function(x){ diff(diff(x)) })
          }
        }
        
        # Non-conformity scores for given data
        # Multivariate functional depth
        # Lower depth is outlier => we take "-" to make nonconformity score
        depth_values <- mfd(arr_X_trans,
                            type = type_depth,
                            alpha = mfd_alpha,
                            depthOptions = depthOptions)
        out <- -as.numeric(depth_values$MFDdepthZ)
        
        return(out)
      }
      # End parallel backend
      stopCluster(cl)
      
      for (s in 1:length(transform)) {
        nonconform_score_indiv[, s] <- res_cv[[s]]
      }
    } else {
      for (s in 1:length(transform)) {
        trans_type <- transform[s]  # transform type
        
        # Transform into 1st or 2nd derivatives
        if (trans_type == "D0") {
          # Raw curves
          arr_X_trans <- arr_X
        } else if (trans_type == "D1") {
          # 1st derivatives
          arr_X_trans <- array(NA, c(m-1, n, p))
          for (i in 1:p) {
            arr_X_trans[, , i] <- apply(arr_X[, , i], 2, diff)
          }
        } else if (trans_type == "D2") {
          # 2nd derivatives
          arr_X_trans <- array(NA, c(m-2, n, p))
          for (i in 1:p) {
            arr_X_trans[, , i] <- apply(arr_X[, , i], 2, function(x){ diff(diff(x)) })
          }
        }
        
        # Non-conformity scores for given data
        # Multivariate functional depth
        # Lower depth is outlier => we take "-" to make nonconformity score
        depth_values <- mfd(arr_X_trans,
                            type = type_depth,
                            alpha = mfd_alpha,
                            depthOptions = depthOptions)
        nonconform_score_indiv[, s] <- -as.numeric(depth_values$MFDdepthZ)
      }
    }
    
    # Weights for weighted average
    if (isTRUE(weight)) {
      weights <- t( apply(nonconform_score_indiv, 1, function(x){ exp(x) / sum(exp(x)) }) )
    } else {
      weights <- 1/length(transform)
    }
    
    # Aggregate scores from transformations
    nonconform_score <- rowSums(nonconform_score_indiv * weights)
  }
  
  # Find training indices without outliers using boxplot
  cutoff <- max(boxplot(nonconform_score, plot = FALSE)$stats)
  idx_clean_null <- which(nonconform_score <= cutoff)
  
  out <- list(
    idx_clean_null = idx_clean_null,
    type = type,
    type_depth = type_depth,
    transform = transform,
    nonconform_score = nonconform_score,
    cutoff = cutoff
  )
  
  # Find training indices without outliers for each transformation
  if (isTRUE(individual) & type == "depth_transform") {
    out$idx_clean_null_indiv <- lapply(1:length(transform), function(i){
      
      # Find outliers using boxplot
      cutoff <- max(boxplot(nonconform_score_indiv[, i], plot = FALSE)$stats)
      idx_clean_null_indiv <- which(nonconform_score_indiv[, i] <= cutoff)
      
      df <- list(
        nonconform_score = nonconform_score_indiv[, i],
        cutoff = cutoff,
        idx_clean_null = idx_clean_null_indiv
      )
      
      return(df)
    })
  }
  
  return(out)
}



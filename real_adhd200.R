library(tidyverse)
library(fdaoutlier)
library(progress)

source("R/foutlier_cp.R")

########################################
### ADHD-200 fMRI Data
########################################
load("RData/ADHD-200_Schaefer17_400regions.RData")

# Fast Fourier Transform with smoothing splines
# Guo, X., Li, Y., & Hsing, T. (2023). Variable Selection and Minimax Prediction in High-dimensional Functional Linear Models. arXiv preprint arXiv:2310.14419.
X_fft <- X
X_fft_sm <- X
for (i in 1:p) {
  print(i)
  X_i_fft <- apply(X[, , i], 1, function(x) {
    Mod(fft(x)) * (2/m)
  })
  
  X_i_fft_sm <- apply(X[, , i], 1, function(x) {
    smooth.spline(Mod(fft(x)) * (2/m))$y
  })
  
  X_fft[, , i] <- t(X_i_fft)
  X_fft_sm[, , i] <- t(X_i_fft_sm)
}


### Simulation for X, X_fft, X_fft_sm
B <- 100  # number of simulations
alpha <- 0.2   # coverage level
n_cores <- 10   # number of cores

res <- list()
for (sim_model_idx in 1:3) {
  if (sim_model_idx == 2) {
    X <- X_fft
  } else if (sim_model_idx == 3) {
    X <- X_fft_sm
  }
  
  n <- dim(X)[1]   # number of curves
  m <- dim(X)[2]   # number of timepoints
  p <- dim(X)[3]  # number of functional variables
  
  
  ### Make additional outliers from ADHD group
  # 1st, 2nd derivatives
  X_deriv_1 <- array(NA, c(m-1, n, p))
  X_deriv_2 <- array(NA, c(m-2, n, p))
  for (i in 1:p) {
    X_deriv_1[, , i] <- apply(X[, , i], 1, diff)
    X_deriv_2[, , i] <- apply(X[, , i], 1, function(x){ diff(diff(x)) })
  }
  
  # True outlier index - c(17, 70)
  df <- data.frame(
    id = 1:nrow(X),
    y = y
  ) %>%
    filter(y == 0)
  data_control <- lapply(1:p, function(i){ X[y == 0, , i] })   # Control group
  
  # Candidates of outliers from ADHD group
  idx_adhd <- which(y == 1)
  # - Choose 20 depth based outliers
  type_depth <- "projdepth"
  depth_values <- list(
    mfd(aperm(X[idx_adhd, , ], c(2,1,3)), type = type_depth, 
        depthOptions = list(type = "Rotation"))$MFDdepthX,
    mfd(X_deriv_1[, idx_adhd, ], type = type_depth, 
        depthOptions = list(type = "Rotation"))$MFDdepthX,
    mfd(X_deriv_2[, idx_adhd, ], type = type_depth, 
        depthOptions = list(type = "Rotation"))$MFDdepthX
  )
  idx_adhd_idx <- lapply(depth_values, function(x){ order(x)[1:20] }) %>% 
    unlist() %>% 
    unique()
  idx_adhd_cand <- idx_adhd[idx_adhd_idx]
  
  
  
  ### Outlier detection with different B splits
  fdr_res <- data.frame(
    marg = rep(NA, B),
    simes = rep(NA, B),
    asymp = rep(NA, B)
  )
  fdr_bh <- list(
    T_projdepth = fdr_res,
    projdepth = fdr_res,
    esssup = fdr_res
  )
  tpr_bh <- fdr_bh
  
  fdr_comparison <- data.frame(
    # ms = rep(NA, B),
    seq = rep(NA, B)
  )
  tpr_comparison <- fdr_comparison
  
  # Progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent [Execution time :elapsedfull || Estimated time remaining: :eta]",
    total = B,   # total number of ticks to complete (default 100)
    clear = FALSE,   # whether to clear the progress bar on completion (default TRUE)
    width = 80   # width of the progress bar
  )
  progress <- function(n){
    pb$tick()
  } 
  
  
  # Repetitions
  for (b in 1:B) {
    set.seed(b)
    
    # Show the progress bar
    progress(b)
    
    # Add the 10 ADHD labels as outliers
    idx_adhd_selected <- sample(idx_adhd_cand, 10)   # 10 sampled ADHD curves
    data <- lapply(1:p, function(i){ rbind(data_control[[i]], X[idx_adhd_selected, , i]) })
    idx_outliers <- (nrow(data[[1]])-9):nrow(data[[1]])
    idx_outliers  # 10 outliers
    
    
    ### Split data into training and test set
    n <- nrow(data[[1]])
    p <- length(data)
    
    prop_train <- 0.8  # proportion of training set
    
    n_train <- round(n * prop_train)   # size of training set
    n_test <- n - n_train   # size of test set
    
    # Split training and test data
    idx_train <- sample(setdiff(1:n, idx_outliers), n_train)
    idx_test <- setdiff(1:n, idx_train)
    
    data_train <- lapply(data, function(x){ x[idx_train, ] })
    data_test <- lapply(data, function(x){ x[idx_test, ] })
    
    
    ### Conformal outlier detection
    # Transformations + projdepth
    obj_T_projdepth <- foutlier_cp(X = data_train, 
                                   X_test = data_test,
                                   type = "depth_transform", 
                                   type_depth = "projdepth",
                                   alpha = alpha,
                                   n_cores = n_cores,
                                   individual = TRUE,
                                   seed = b)
    fdr_bh$T_projdepth[b, ] <- sapply(obj_T_projdepth$idx_out, function(x){
      get_fdr(idx_test[x], idx_outliers)
    })
    tpr_bh$T_projdepth[b, ] <- sapply(obj_T_projdepth$idx_out, function(x){
      get_tpr(idx_test[x], idx_outliers)
    })
    
    # esssup
    obj_esssup <- foutlier_cp(X = data_train, 
                              X_test = data_test,
                              type = "esssup",
                              alpha = alpha,
                              seed = b)
    fdr_bh$esssup[b, ] <- sapply(obj_esssup$idx_out, function(x){
      get_fdr(idx_test[x], idx_outliers)
    })
    tpr_bh$esssup[b, ] <- sapply(obj_esssup$idx_out, function(x){
      get_tpr(idx_test[x], idx_outliers)
    })
    
    # projdepth
    fdr_bh$projdepth[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[1]], function(x){
      get_fdr(idx_test[x], idx_outliers)
    })
    tpr_bh$projdepth[b, ] <- sapply(obj_T_projdepth$idx_out_indiv[[1]], function(x){
      get_tpr(idx_test[x], idx_outliers)
    })
    
    
    ### Existing functional outlier detection (Coverage guarantee X)
    idx_comparison <- list(
      # ms = c(),
      seq = c()
    )
    arr_train <- abind::abind(data_train, along = 3)
    arr_test <- abind::abind(data_test, along = 3)
    
    # Parallel computation
    cl <- makeCluster(n_cores)
    registerDoSNOW(cl)
    pkgs <- c("fdaoutlier")
    res_cv <- foreach(i = 1:n_test, .packages = pkgs) %dopar% {
      df <- array(NA, dim = c(n_train+1, m, p))
      df[1:n_train, , ] <- arr_train
      df[n_train+1, , ] <- arr_test[i, , ]
      
      out <- list()
      
      # # MS plot
      # outlier_ms <- msplot(dts = df, plot = F, seed = b)$outliers
      # if (length(outlier_ms) > 0 & ((n_train+1) %in% outlier_ms)) {
      #   out$ms <- idx_test[i]
      # } else {
      #   out$ms <- integer(0)
      # }
      
      # Sequential transformation
      seqobj <- seq_transform(df,
                              sequence = c("O","D1","D2"),
                              depth_method = "erld",
                              erld_type = "one_sided_right",
                              seed = b)
      outlier_seq <- unique(unlist(seqobj$outliers))
      if (length(outlier_seq) > 0 & ((n_train+1) %in% outlier_seq)) {
        out$seq <- idx_test[i]
      } else {
        out$seq <- integer(0)
      }
      
      return(out)
    }
    # End parallel backend
    stopCluster(cl)
    
    # idx_comparison$ms <- unlist(sapply(res_cv, function(x){ x$ms }))
    idx_comparison$seq <- unlist(sapply(res_cv, function(x){ x$seq }))
    
    
    fdr_comparison[b, ] <- sapply(idx_comparison, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr_comparison[b, ] <- sapply(idx_comparison, function(x){
      get_tpr(x, idx_outliers)
    })
  }
  
  # Results
  res[[sim_model_idx]] <- list(
    bh = list(fdr = fdr_bh,
              tpr = tpr_bh),
    comparison = list(fdr = fdr_comparison,
                      tpr = tpr_comparison)
  )
  
}
save(res, file = "RData/adhd_200_res.RData")



# Summary the results
res2 <- list()
for (i in 1:length(res)) {
  res2[[i]] <- list(
    fdr = cbind(res[[i]]$bh$fdr,
                res[[i]]$comparison$fdr),
    tpr = cbind(res[[i]]$bh$tpr,
                res[[i]]$comparison$tpr)
  )
}

lapply(res2, function(sim){
  sub <- paste0(
    rbind(fdr = colMeans(sim$fdr),
          tpr = colMeans(sim$tpr)) %>% 
      round(3) %>% 
      format(nsmall = 3),
    " (",
    rbind(fdr = apply(sim$fdr, 2, sd),
          tpr = apply(sim$tpr, 2, sd)) %>% 
      round(3) %>% 
      format(nsmall = 3),
    ")"
  )
  dim(sub) <- c(2, ncol(sim$fdr))
  rownames(sub) <- c("FDR","TPR")
  colnames(sub) <- colnames(sim$fdr)
  sub <- data.frame(sub)
  sub
})






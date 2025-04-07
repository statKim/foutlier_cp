######################################################
### Simulation for the mixed training setting
### - Outlier detection for mixed training set
######################################################
library(tidyverse)
library(fdaoutlier)
library(progress)
source("R/foutlier_cp.R")

n <- 1000   # number of training data (proper training + calibration)
m <- 51
p <- 20

B <- 100   # number of repetitions
outlier_rate <- 0.2   # proportion of training outliers
alpha <- 0.1  # coverage level

sim_model <- 1:4  # simulation models

# Simulation
res <- list()
for (sim_model_idx in 1:length(sim_model)) {
  print(sim_model_idx)
  
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
  
  
  # Simulation for each simulation model
  fdr <- data.frame(
    T_projdepth = rep(NA, B),
    projdepth = rep(NA, B),
    esssup = rep(NA, B),
    ms = rep(NA, B),
    seq = rep(NA, B)
  )
  tpr <- fdr
  
  for (b in 1:B) {
    set.seed(b)
    
    # Show the progress bar
    progress(b)
    
    ### Data generation
    # Generate multivariate functional data with outliers
    data_obj <- foutlier_sim_mfd(n = n, m = m, p = p, outlier_rate = outlier_rate, 
                                 model = sim_model_idx)
    # data_obj$idx_outliers
    # plot(data_obj, p = 1)
    data_train <- data_obj$data
    idx_outliers <- data_obj$idx_outliers
    
    # Transformations + projdepth
    obj_T_projdepth <- get_clean_null(data_train, 
                                      type = "depth_transform",
                                      type_depth = "projdepth",
                                      n_cores = 3)
    fdr$T_projdepth[b] <- get_fdr(setdiff(1:n, obj_T_projdepth$idx_clean_null), idx_outliers)
    tpr$T_projdepth[b] <- get_tpr(setdiff(1:n, obj_T_projdepth$idx_clean_null), idx_outliers)
    
    # esssup
    obj_esssup <- get_clean_null(data_train, type = "esssup")
    fdr$esssup[b] <- get_fdr(setdiff(1:n, obj_esssup$idx_clean_null), idx_outliers)
    tpr$esssup[b] <- get_tpr(setdiff(1:n, obj_esssup$idx_clean_null), idx_outliers)
    
    # projdepth
    fdr$projdepth[b] <- get_fdr(setdiff(1:n, obj_T_projdepth$idx_clean_null_indiv[[1]]$idx_clean_null),
                                idx_outliers)
    tpr$projdepth[b] <- get_tpr(setdiff(1:n, obj_T_projdepth$idx_clean_null_indiv[[1]]$idx_clean_null),
                                idx_outliers)
    
    
    ### Existing functional outlier detection
    idx_comparison <- list()
    
    df <- abind::abind(data_train, along = 3)
    # MS plot
    idx_comparison$ms <- msplot(dts = df, plot = F)$outliers
    
    # Sequential transformation
    seqobj <- seq_transform(df, 
                            sequence = c("O","D1","D2"),
                            depth_method = "erld",
                            erld_type = "one_sided_right", 
                            save_data = F)
    idx_comparison$seq <- unique(unlist(seqobj$outliers))
    
    
    fdr[b, c("ms","seq")] <- sapply(idx_comparison, function(x){
      get_fdr(x, idx_outliers)
    })
    tpr[b, c("ms","seq")] <- sapply(idx_comparison, function(x){
      get_tpr(x, idx_outliers)
    })
  }
  
  # Results
  res[[sim_model_idx]] <- list(
    fdr = fdr,
    tpr = tpr
  )
}
save(res, file = paste0("RData/sim_p", p, "_n_", n, "_mixed_train_outlier.RData"))


# Summary the results
res2 <- list()
for (i in 1:length(res)) {
  res2[[i]] <- list(
    fdr = res[[i]]$fdr,
    tpr = res[[i]]$tpr
  )
}
res3 <- lapply(res2, function(sim){
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
res3

rbind(res3[[1]], res3[[2]], res3[[3]], res3[[4]]) %>% 
  t() %>% 
  xtable::xtable()


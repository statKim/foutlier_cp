# Get FDR (False discovery rate)
get_fdr <- function(idx_reject, idx_true) {
  if (length(idx_reject) == 0) {
    return(0)
  } else {
    sum(!(idx_reject %in% idx_true)) / length(idx_reject)
  }
}


# Get TPR (True positive rate; Power)
get_tpr <- function(idx_reject, idx_true) {
  if (length(idx_reject) == 0) {
    return(0)
  } else {
    sum(idx_true %in% idx_reject) / length(idx_true)
  }
}


# BH procedure
BH <- function(pvalue, alpha = 0.1) {
  x <- pvalue
  n_test <- length(pvalue)
  
  if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
    idx_rej <- integer(0)
  } else {
    idx_rej <- order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
  }
  
  return(idx_rej)
}



# Generate simulated multivariate functional data
foutlier_sim_mfd <- function(n, m = 51, p = 20, outlier_rate, model = 1, rho = 0.3, ...) {
  # Generate multivarate functional data using `fdaoutlier` package
  sim_ftn_list <- list(
    function(...){ simulation_model1(q = 2, ...) },
    function(...){ simulation_model2(q = 2, ...) },
    function(...){ simulation_model3(q = 1.5, ...) },
    function(...){ simulation_model5(cov_alpha2 = 0.5, ...) }
  )
  
  sim_ftn <- sim_ftn_list[[model]]   # generating function
  
  # Generate multivariate functional data
  data_list <- list()
  if (outlier_rate == 0) {
    # Generate multivariate functional data without outliers
    for (j in 1:p) {
      sim_obj <- sim_ftn(n = n, p = m, outlier_rate = 0)
      data_list[[j]] <- sim_obj$data
    }
    idx_outliers <- NULL
  } else if (outlier_rate > 0) {
    # Generate multivariate functional data with outliers
    idx_outliers <- (n - n*outlier_rate + 1):n
    data_test <- list()
    for (j in 1:p) {
      sim_obj <- sim_ftn(n = n, p = m, outlier_rate = outlier_rate)
      sim_data_p <- sim_obj$data
      
      data_list[[j]] <- matrix(0, n, m)
      # non-outliers
      data_list[[j]][-idx_outliers, ] <- sim_data_p[-sim_obj$true_outliers, ]
      # outliers
      data_list[[j]][idx_outliers, ] <- sim_data_p[sim_obj$true_outliers, ]
    }
  }
  
  out <- list(
    data = data_list,
    idx_outliers = idx_outliers
  )
  
  class(out) <- "foutlier_sim_mfd"
  
  return(out)
}




# Plot for `foutlier_sim_data` object
plot.foutlier_sim_mfd <- function(obj, p = 1, xlabel = "", ylabel = "", plot_title = NULL,
                                  title_cex = 1.5, show_legend = TRUE, legend_pos = "bottomright") {
  data <- obj$data[[p]]
  idx_outliers <- obj$idx_outliers
  
  n <- nrow(data)
  p <- ncol(data)
  gr <- seq(0, 1, length.out = p)
  
  if (length(idx_outliers) > 0) {
    data_null <- t(data[-idx_outliers, , drop = F])
    data_out <- t(data[idx_outliers, , drop = F])
    
    plot(x = gr, type = "n", ylab = ylabel, xlab = xlabel,
         ylim = range(data) + c(-.5*sd(data[, p]), .5*sd(data[, p])),
         col.lab = "gray20", axes = F)
    # grid(col = "grey75", lwd = .3)
    grid()
    matlines(data_null,
             col = "grey61",
             lty = "solid",
             lwd = .4)
    matlines(data_out,
             col = "#D55E00",
             lty = "solid",
             lwd = 1.3)
  } else {
    plot(x = gr, type = "n", ylab = ylabel, xlab = xlabel,
         ylim = range(data) + c(-.5*sd(data[, p]), .5*sd(data[, p])),
         col.lab = "gray20", axes = F)
    # grid(col = "grey75", lwd = .3)
    grid()
    matlines(t(data),
             col = "grey61",
             lty = "solid",
             lwd = .4)
  }
  
  axis(1, at = seq(1, m, length.out = 5), labels = seq(0, 1, length.out = 5), 
       col = "white", col.ticks = "grey61",
       lwd.ticks = .5, tck = -0.025,
       cex.axis = 0.9, col.axis = "gray30")
  axis(2, col = "white", col.ticks = "grey61",
       lwd.ticks = .5, tck = -0.025,
       cex.axis = 0.9, col.axis = "gray30")
  
  box(col = "grey51")
  if(show_legend){
    legend(legend_pos, legend = c("normal", "outlier"),
           lty = c("solid", "solid"),
           lwd = c(.4, 1.3),
           col = c("grey61", "#D55E00"),
           text.col = "gray40", bty = "n",
           box.lwd = .1, xjust = 0, inset = .01)
  }
  if (!is.null(plot_title)) {
    mtext(plot_title, 3, adj = 0.5, line = 1, cex = title_cex,
          col = "gray20")
  }
}

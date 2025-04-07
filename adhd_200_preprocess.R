library(tidyverse)

########################################
### ADHD-200 fMRI Data
########################################
# Phenotypic data
path <- "../../../ADHD-200/"
data_pheno <- read_tsv(paste0(path, "adhd200_preprocessed_phenotypics.tsv"))

# Load training ADHD-200 data
file_list <- list.files(paste0(path, "nilearn_preprocess/Peking_train/"))
subject_id <- rep(NA, length(file_list))
for (i in 1:length(file_list)) {
  df <- read_csv(paste0(path, "nilearn_preprocess/Peking_train/",  file_list[i]), show_col_types = F)
  
  if (i == 1) {
    region_name <- data.frame(RegionName = df$RegionName)
    X <- array(0, c(length(file_list), ncol(df)-1, nrow(df)))
  } else {
    df <- left_join(region_name, df, by = "RegionName")
  }
  
  # Subject id
  subject_id[i] <- strsplit(file_list[i], ".csv") %>% 
    unlist() %>% 
    as.integer()
  
  # Brain signals  
  X[i, , ] <- t( as.matrix(df[, -1]) )
}
dim(X)

# Remove fail from quality control
pheno_sub <- tibble(`ScanDir ID` = subject_id) %>% 
  left_join(data_pheno, by = "ScanDir ID")
idx_qc <- which(pheno_sub$QC_Athena == 1)
pheno_sub <- pheno_sub[idx_qc, ]
X <- X[idx_qc, , ]
dim(X)

# ADHD index (0 is control; 1 is ADHD)
y <- ifelse(pheno_sub$DX == 0, 0, 1)

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables


save(X, y, file = "RData/ADHD-200_Schaefer17_400regions.RData")

# matplot(t(X[1:30, , 10]), type = "l", col = ifelse(y == 1, 2, "gray"), lty = 1)

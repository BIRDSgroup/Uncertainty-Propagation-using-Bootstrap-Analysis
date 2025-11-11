library(data.table)
library(foreach)
library(doParallel)

# ==================================================
# INPUT DATA
cores <- 5

in_file <- 'Muscle_Skeletal_recount normalized.csv'
cov_file <- 'covariates.csv'
gene_file <- 'genes.csv'

out_file <- 'Muscle_Skeletal_recount.csv'

# # ==================================================

df <- read.csv(in_file)
cov_data <- read.csv(cov_file,header = TRUE)

# Gene-wise Covariate analysis
cluster <- makeCluster(cores)
registerDoParallel(cluster)
df_adjusted <- foreach(i = 1:nrow(df), .combine = 'rbind') %dopar% {
  lm_model <- lm( as.numeric(df[i, -1]) ~ cov_data[,2] + cov_data[,3])
  residuals(lm_model)
}
stopCluster(cluster)
df_adjusted <- as.data.frame(df_adjusted)
df_adjusted <- cbind(df[, 1], df_adjusted)
colnames(df_adjusted) <- colnames(df)
write.csv(df_adjusted, out_file, row.names = FALSE, sep = ",")

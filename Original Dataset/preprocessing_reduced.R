# ==================================================
# INPUT DATA
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]

dir_path <- paste('Preprocessed Files/', tissue, '/', sep = "")
in_file <- paste(dir_path, tissue, '.csv', sep = '')
# # ==================================================

seed <- as.integer(Sys.time())
write(seed, paste(dir_path, '/seed.txt', sep = ''))
set.seed(seed)

df_adjusted <- read.csv(in_file, header = TRUE, sep = ",", row.names = 1)
# Sampling
for (s in c(73, 237)) {
  cat('For', s, 'samples\n')
  
  sampled_cols <- sample(colnames(df_adjusted), s, replace = FALSE)
  df_sampled <- df_adjusted[,sampled_cols]
  out_file <- paste(dir_path, tissue, '_', s, '.csv', sep = '')
  write.csv(df_sampled, out_file, row.names = TRUE, sep = ",")
}

library(data.table)
library(foreach)
library(doParallel)

# ==================================================
# INPUT DATA
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
cores <- as.numeric(args[2])

in_file <- paste('Gene Expression Matrices/', tissue, '.v8.normalized_expression.bed.gz', sep = '')
cov_file <- paste('Covariates/', tissue, '.v8.covariates.txt', sep = '')
gene_file <- paste('Genes/', tissue, '.v8.egenes.txt.gz', sep = '')

dir_path <- paste('Preprocessed Files/', tissue, '/', sep = "")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
out_file <- paste(dir_path, tissue, '.csv', sep = '')
outg_file <- paste(dir_path, 'genes.csv', sep = '')
outr_file <- paste(dir_path, tissue, '_orig_filtered.csv', sep = '')

# # ==================================================

# Extract Protein-Coding Genes
print('Extracting Protein-Coding Genes')
gencode <- rtracklayer::import("gencode.v26.GRCh38.genes.gtf")
gencode <- as.data.frame(gencode)
gencode <- gencode[gencode$type == "gene" & gencode$gene_type == "protein_coding",]

df <- fread(in_file, sep = '\t', header = TRUE)  
genes <- fread(gene_file, sep = '\t', header = TRUE)
genes_filtered <- genes[which(paste(genes$gene_name, genes$gene_id) %in% paste(gencode$gene_name, gencode$gene_id)), ]
df_filtered <- df[df$gene_id %in% genes_filtered$gene_id,]

write.csv(df_filtered, outr_file, row.names = FALSE, sep = ",")
write.csv(genes_filtered, outg_file, row.names = FALSE, sep = ",")

# Covariates
print('Covariate Analysis')
cov_data <- read.table(cov_file,header = TRUE)
cov_data <- cov_data[,2:ncol(cov_data)]
df_filtered <- df_filtered[, 5:ncol(df_filtered)]

if (ncol(df_filtered) != ncol(cov_data)){
  stop("Error: covariate samples mismatch")
}

# Gene-wise Covariate analysis
cluster <- makeCluster(cores)
registerDoParallel(cluster)
df_adjusted <- foreach(i = 1:nrow(df_filtered), .combine = 'rbind') %dopar% {
  lm_model <- lm( as.numeric(df_filtered[i, ]) ~ t(cov_data))
  residuals(lm_model)
}
stopCluster(cluster)
df_adjusted <- as.data.frame(df_adjusted)
row.names(df_adjusted) <- genes_filtered$gene_id

write.csv(df_adjusted, out_file, row.names = TRUE, sep = ",")

cat(tissue, '\t', nrow(df_adjusted), '\t', ncol(df_adjusted), '\n')

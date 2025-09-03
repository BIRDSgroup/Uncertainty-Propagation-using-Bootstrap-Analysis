# ==================================================
# INPUT
args = commandArgs(trailingOnly = TRUE)
tissue = args[1]
B = as.numeric(args[2])
# ==================================================

# Seed Extraction
seed = as.integer(Sys.time())
write(seed, paste(tissue, '/seed.txt', sep = ''))
set.seed(seed)

# Seed Generation for all samples
sample_seeds = sample(x=.Machine$integer.max, size=B, replace=T)
for (i in 1:B) {
  write(sample_seeds[i], paste( 'seeds/seed_', i-1, sep = '') )
}
# Extract Protein-Coding Genes
print('Extracting Protein-Coding Genes')
gencode <- rtracklayer::import("../Original Dataset/gencode.v26.GRCh38.genes.gtf")
gencode <- as.data.frame(gencode)
gencode <- gencode[gencode$type == "gene" & gencode$gene_type == "protein_coding",]
write.csv(gencode, 'gencode.csv', row.names = FALSE, sep = ',')

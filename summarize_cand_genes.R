# This script aggregates candidate gene analysis data from all QTL 
# into one table for publication 
suppressMessages(library(readr))
suppressMessages(library(stringi))

summaries <- list.files(path = "results/cand-gene-analysis/summary")
results <- list.files(path = "results/cand-gene-analysis/vardata")

numQTL <- 8
# dataframe: num genes, num SNPs, num genes with reg variants, 
# num genes with protein-coding variants, num genes with both, 
# num genes related to allergy (GSEA)? 
summary <- data.frame('QTL' = vector(length = numQTL),
                      'cand_genes' = vector(length = numQTL),
                      'num_variants' = vector(length = numQTL),
                      'pc_count' = vector(length = numQTL),
                      'reg_count' = vector(length = numQTL))

for (i in 1:numQTL){
  sum <- readLines(paste("results/cand-gene-analysis/summary/", summaries[i], sep=""))
  num_genes <- stri_split_charclass(sum[length(sum)-1], '\\p{WHITE_SPACE}')[[1]][5]
  num_vars <- stri_split_charclass(sum[length(sum)], '\\p{WHITE_SPACE}')[[1]][5]
  
  vardata <- readr::read_csv(paste("results/cand-gene-analysis/vardata/", results[i], sep=""))
  pc_count <- sum(vardata$pc_var_count>0) # genes with only pc variants or pc + reg variants 
  reg_count <- sum((vardata$reg_var_count>0)&(vardata$pc_var_count==0)) # genes with only reg variants 
  
  summary[i,] <- c(paste("Qpa", i, sep=""), num_genes, num_vars, pc_count, reg_count)
}

# include total gene counts per QTL 
gene_counts <- readr::read_csv('results/qtl-gene-counts.csv')
summary$total_genes <- gene_counts$gene_count
# re-order columns 
summary <- summary[,c('QTL', 'total_genes', 'cand_genes', 'num_variants', 'pc_count', 'reg_count')]

write_csv(summary, file = "results/cand-genes-summary.csv")


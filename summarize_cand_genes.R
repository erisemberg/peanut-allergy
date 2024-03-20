# This script aggregates candidate gene analysis data from all QTL 
# into one table for publication 
library(readr)
library(readxl)
library(stringi)
library(dplyr)

summaries <- list.files(path = "results/summary")
results <- list.files(path = "results/vardata")

# load mouse2human gene mape
m2hgenemap <- read_csv('derived_data/mouse_human_gene_map.csv')
m2hgenemap <- m2hgenemap[,2:3]
names(m2hgenemap)[1] <- 'gene_name' # for merging with vardata 

# load list of peanut genes 
pnutgenes <- read_xlsx('source_data/peanut-genes.xlsx')
high_ev_pnutgenes <- pnutgenes$Gene[which(pnutgenes$`Level of evidence` %in% c(1,2,3))]
names(pnutgenes)[1] <- 'human'

numQTL <- 8
# dataframe: num genes, num SNPs, num genes with reg variants, 
# num genes with protein-coding variants, num genes with both, 
# num genes related to allergy (GSEA)? 
summary <- data.frame('QTL' = vector(length = numQTL),
                      'cand_genes' = vector(length = numQTL),
                      'num_variants' = vector(length = numQTL),
                      'pc_count' = vector(length = numQTL),
                      'reg_count' = vector(length = numQTL),
                      'mapped_human_genes' = vector(length = numQTL))

all_pnut_genes <- data.frame(QTLid = character(), pnutgene = character())

for (i in 1:numQTL){
  sum <- readLines(paste("results/summary/", summaries[i], sep=""))
  num_genes <- stri_split_charclass(sum[length(sum)-1], '\\p{WHITE_SPACE}')[[1]][5]
  num_vars <- stri_split_charclass(sum[length(sum)], '\\p{WHITE_SPACE}')[[1]][5]
  
  varfname <- paste("results/vardata/", results[i], sep="")
  vardata <- readr::read_csv(varfname)
  pc_count <- sum(vardata$pc_var_count>0) # genes with only pc variants or pc + reg variants 
  reg_count <- sum((vardata$reg_var_count>0)&(vardata$pc_var_count==0)) # genes with only reg variants 
  
  vardata <- left_join(vardata, m2hgenemap, by = 'gene_name') # map to human gene names 
  write_csv(vardata, file = varfname) 
  
  # are there mapped human genes that are also high evidence pnut genes?
  mapd <- unique(vardata$human[which(vardata$human %in% high_ev_pnutgenes)])
  
  # are there mapped human genes that are pnut genes of any evidence? 
  mapd2 <- vardata[which(vardata$human %in% pnutgenes$human),c('gene_name','human')]
  mapd2 <- left_join(mapd2, pnutgenes, by = 'human')
  mapd2 <- mapd2[,c(1:3)]
  qtl_name <- paste("Qpa", i, sep = "")
  mapd2$QTL = rep(qtl_name, nrow(mapd2))
  mapd2 <- mapd2[,c('QTL', 'gene_name', 'human', 'Source')]
  mapd2 <- mapd2 %>% distinct()
  
  # save high evidence genes to summary table 
  summary[i,] <- c(paste("Qpa", i, sep=""), num_genes, num_vars, pc_count, reg_count, paste(mapd, collapse = ", "))
  
  # save all mapped peanut genes to supp table
  all_pnut_genes <- rbind(all_pnut_genes, mapd2)
}

# include total gene counts per QTL 
gene_counts <- readr::read_csv('results/qtl-gene-counts.csv')
summary$total_genes <- gene_counts$gene_count
# re-order columns 
summary <- summary[,c('QTL', 'total_genes', 'cand_genes', 'num_variants', 'pc_count', 'reg_count', 'mapped_human_genes')]

write_csv(summary, file = "results/cand-genes-summary.csv")
write_csv(all_pnut_genes, file = 'results/QTL-pnut-genes.csv')

library(rtracklayer)
library(readr)
#library(RCurl)

# Download from FTP 
# ftp_url <- "https://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/Mus_musculus.GRCm38.97.gtf.gz"
  
gtf <- rtracklayer::import('Mus_musculus.GRCm38.79.gtf')
genemap <- mcols(gtf)[,c('gene_id', 'gene_name')]
genemap <- as.data.frame(genemap) 

genemap <- genemap[!duplicated(genemap[1:2]),] # 1869634 -> 55487 rows 

write_csv(genemap, "GRCm38-gtf-genemap.csv")

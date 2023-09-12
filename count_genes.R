### This script counts genes within the QTL (without candidate gene analysis) 
suppressMessages(library(rtracklayer))
suppressMessages(library(readr))
suppressMessages(library(GenomicFeatures))

#' BED to GRanges
#'
#' This function loads a BED-like file and stores it as a GRanges object.
#' The tab-delimited file must be ordered as 'chr', 'start', 'end', 'id', 'score', 'strand'.
#' The minimal BED file must have the 'chr', 'start', 'end' columns.
#' Any columns after the strand column are ignored.
#' 
#' @param file Location of your file
#' @keywords BED GRanges
#' @export
#' @examples
#' bed_to_granges('my_bed_file.bed')

bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  library("GenomicRanges")
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}




#------------------------------------------------------------------------------#
# Download from FTP 
# ftp_url <- "https://ftp.ensembl.org/pub/release-79/gtf/mus_musculus/Mus_musculus.GRCm38.79.gtf.gz"

gtf <- rtracklayer::import('source_data/Mus_musculus.GRCm38.79.gtf')
print('Imported Mus_musculus.GRCm38.79.gtf')
gtf_genes <- gtf[gtf$type == 'gene'] # filter down to only genes  

qtls <- bed_to_granges('source_data/qtl-regions.bed')
print('Imported qtl-regions.bed')

hits <- findOverlaps(query = qtls, subject = gtf_genes)

qtlids <- paste0('Qpa', seq(1,8))
gene_counts <- data.frame(QTL = qtlids, gene_count = as.table(hits))
print('Genes counted')
write_csv(gene_counts, file = 'results/qtl-gene-counts.csv')
print('Results saved to qtl-gene-counts.csv')






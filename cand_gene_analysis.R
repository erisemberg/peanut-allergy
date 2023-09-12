# This script performs candidate gene filtering 
suppressMessages(library(tidyverse))
library(readr)
library(stringi)
source('code-dependencies/qtl_functions.R')
source('code-dependencies/cmdline.R')

#-----------------------------------Setup--------------------------------------#
# input
vcffiles <- cmdline.option("vcf")
# split up input by newline character if more than 1 passed in 
vcffiles <- stri_split_fixed(vcffiles, ',', omit_empty = TRUE)[[1]]

# output
prefix <- cmdline.option("prefix")

ensure_directory('results/cand-gene-analysis')
ensure_directory('results/cand-gene-analysis/vardata')
outfile <- paste('results/cand-gene-analysis/vardata/', prefix, '-vardata.csv', sep="")

ensure_directory('results/cand-gene-analysis/summary')
logfile <- paste('results/cand-gene-analysis/summary/', prefix, '.txt', sep="")
log <- make_logger(filename = logfile)

log(paste("VCF files: ", vcffiles))

pc_var_types <- c('missense_variant', 'splice_region_variant', 'stop_gained',
                  'stop_lost', 'splice_donor_variant', 'splice_acceptor_variant',
                  'inframe_deletion', 'initiator_codon_variant',  'coding_sequence_variant', 
                  'incomplete_terminal_codon_variant', 'inframe_insertion', 'frameshift_variant',
                  'splice_')
reg_var_types <- c('downstream_gene_variant', 'intron_variant', 'upstream_gene_variant',
                   'non_coding_transcript_exon_variant', 'non_coding_transcript_variant',
                   '3_prime_UTR_variant', '5_prime_UTR_variant', 'synonymous_variant',
                   'NMD_transcript_variant', 'mature_miRNA_variant', 'stop_retained_variant')


#----------------------------------Analysis------------------------------------#
vardata_by_gene <- data.frame('gene_id' = vector(), 'var_count' = vector(), 
                              'pc_var_count' = vector(), 'reg_var_count' = vector())
vars_counted <- vector()

for (vcffile in vcffiles){ # for each region in QTL 
  vcf <- read.csv(file = vcffile, header = TRUE)
  
  for (i in 1:nrow(vcf)){ # for each variant 
    #variant_id <- vcf$ID[i]
    chrpos <- paste(vcf$CHROM[i], vcf$POS[i], sep="-")
    # if variant already analyzed (likely in prior region of same QTL), don't need 
    # to repeat analysis 
    if (chrpos %in% vars_counted){ next } 
    
    annotation <- vcf$ANN[i] 
    
    # if annotation column doesn't start with an ensembl ID, then it's probably an
    # intergenic variant, skip it 
    if (!(startsWith(annotation, 'ENSMUSG'))) { next }
    
    ann_entries <- stri_split_regex(annotation, pattern = '(?=ENSMUSG)', omit_empty = TRUE)[[1]]
    newdata <- as_tibble(stri_split_regex(ann_entries, pattern = '\\|', simplify = TRUE),
                         .name_repair = 'minimal')[,c(1,4)]
    colnames(newdata) <- c('gene_id', 'variant_type')
    
    newdata <- newdata %>% 
      separate(col = 'variant_type', into = c('type1', 'type2', 'type3'), sep = '&', extra = 'warn', 
               fill = 'right') %>%
      pivot_longer(cols = c(2:4), names_to = NULL, values_to = 'variant_type')
    
    # remove rows with variant_type == NA (extra rows produced by separate)
    newdata <- newdata[!is.na(newdata$variant_type),]
    # remove intergenic variants (we don't care about those)
    newdata <- newdata[newdata$variant_type != 'intergenic_variant',]
    # remove duplicate gene/variant_type combinations (we just want a list of variant types by gene)
    newdata <- newdata[!duplicated(newdata[c(1:2)]),] 
    
    # add variant to variant counts by gene. Each row in VCF file is one variant; 
    # if that variant affects multiple genes (unc), should count for each gene 
    unique_genes <- unique(newdata$gene_id)
    for (j in 1:length(unique_genes)){ # for each (unique) gene in newdata (usually 1)
      cur_gene <- unique_genes[j]
      if (!(cur_gene %in% vardata_by_gene$gene_id)){ # first time gene has appeared 
        # add new row with gene ID, variant count of 1 and pc/reg variant counts of 0 
        vardata_by_gene[nrow(vardata_by_gene)+1,] <- c(cur_gene, 1, 0, 0) 
      } else if (cur_gene %in% vardata_by_gene$gene_id){ # gene already in dataframe; increment var count  
        vardata_by_gene[vardata_by_gene$gene_id == cur_gene, 'var_count'] <- 
          as.integer(vardata_by_gene[vardata_by_gene$gene_id == cur_gene, 'var_count']) + 1
      }
      
      # for each variant type listed for each gene, increment either pc_var_count or 
      # reg_var_count depending on group variant type belongs to. Some genes have 
      # more than 1 variant type associated with a particular SNP - count will be 
      # incremented once for each. For example, if there are 3 regulatory variants 
      # and 1 protein-affecting variant, pc_var_count and reg_var_count will each be 
      # incremented by 1. 
      vartmp <- newdata[newdata$gene_id == cur_gene,]
      pc_var_count_inc <- FALSE # indicator for whether pc_var_count has been incremented for this variant 
      reg_var_count_inc <- FALSE # indicator for whether reg_var_count has been incremented for this variant
      
      for (k in 1:nrow(vartmp)){
        cur_var_type <- vartmp$variant_type[k]
        
        if (!(cur_var_type %in% pc_var_types) & !(cur_var_type %in% reg_var_types)){
          log(paste("variant type not categorized:", cur_var_type))
          next
        }
        
        if ((cur_var_type %in% pc_var_types) & !(pc_var_count_inc)){ # if protein-affecting variant type, increment pc var count 
          vardata_by_gene[vardata_by_gene$gene_id == cur_gene, 'pc_var_count'] <- 
            as.integer(vardata_by_gene[vardata_by_gene$gene_id == cur_gene, 'pc_var_count']) + 1
          pc_var_count_inc <- TRUE
        } else if ((cur_var_type %in% reg_var_types) & !(reg_var_count_inc)){ # if regulatory variant type, increment reg var count 
          vardata_by_gene[vardata_by_gene$gene_id == cur_gene, 'reg_var_count'] <- 
            as.integer(vardata_by_gene[vardata_by_gene$gene_id == cur_gene, 'reg_var_count']) + 1
          reg_var_count_inc <- TRUE 
        }
      }
    }
    
    # list variant as counted so we don't recount it in subsequent regions in same QTL 
    vars_counted <- append(vars_counted, chrpos) 
  }
}

# map ensembl IDs to gene names 
genemap <- readr::read_csv('derived_data/GRCm38-gtf-genemap.csv', show_col_types = FALSE)

vardata_by_gene <- merge(vardata_by_gene, genemap, all.x = TRUE, by = 'gene_id')
vardata_by_gene <- vardata_by_gene[,c(5,1:4)]


#----------------------------------Output--------------------------------------#
write_csv(vardata_by_gene, file = outfile)

log(paste("Number of genes: ", length(unique(vardata_by_gene$gene_id)))) # total genes in QTL
log(paste("Number of variants: ", length(unique(vars_counted)))) # total variants in all genes in QTL 


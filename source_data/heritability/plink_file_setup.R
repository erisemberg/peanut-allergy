# This script processes genotype data into PLINK format for estimation of 
# genomic relationship matrix 

library(tidyverse)
library(readxl)
library(stringi)
source("code-dependencies/qtl_functions.R")

#---------------------------------Parameters-----------------------------------#
geno_file <- "source_data/Cr_WB02_miniMUGA-06242021.csv"
pheno_file <- "source_data/CC27xC3H_backcross_pheno_clean.xlsx"
ensure_directory("derived_data")
out_file <- "derived_data/genos_for_plink.csv"
#Y_MT_out_file <- "derived_data/Geno_CC27xC3H_Y_MT.csv"

num_F2s = 365 # number of genotyped F2 (or BC) mice
crosstype = "bc" # cross type (f2 or bc)
A.ref = "CC027.GeniUnc" # parent mouse representing A genotype 
B.ref = "C3H.HeJ" # parent mouse representing B genotype 

#------------------------------------------------------------------------------#
#-------------------------------Genotype data----------------------------------#
#------------------------------------------------------------------------------#
ped <- read.csv(geno_file)

F2_start = ncol(ped)-num_F2s+1 # column index of first F2 mouse column
F2_end = ncol(ped) # column index of last F2 mouse column 

ped <- ped[ped$Chromosome != 0, ] # Remove chr0 rows
ped$Position_b38 <- ped$Position_b38/1e6 # convert bp position to Mb (2 Mb = 1 cM)

#----------------------------Metric calculation--------------------------------#
# Get ref, alt, het and N calls in the F2 mice at each marker
ped$num_ref <- rowSums(ped[,F2_start:F2_end] == ped$reference)
ped$num_alt <- rowSums(ped[,F2_start:F2_end] == ped$alternate)
ped$num_N <- rowSums(ped[,F2_start:F2_end] == "N")
ped$pct_N <- ped$num_N/num_F2s
ped$num_H <- rowSums(ped[,F2_start:F2_end] == "H")

ped$ref_alt <- ped$num_ref/(ped$num_ref + ped$num_alt) # Calculate ref as proportion of ref+alt
ped$het_all <- ped$num_H/(ped$num_ref + ped$num_alt + ped$num_H) # Calculate het relative to good (non-N) calls

#--------------------------------Filtering-------------------------------------#
### Separate out MT and Y markers // what to do with PAR? 
ped.Y.MT <- ped[which((ped$Chromosome == "Y") | ped$Chromosome == "MT"), ]
### Filter Y chromosome markers - all females should have N calls for each Y marker
#ped.Y.MT <- ped[which((ped$Chromosome == "Y") & (ped$pct_N == 0))]
ped <- ped[which((ped$Chromosome != "Y") & (ped$Chromosome != "MT")),]

# Filter out markers with failure rate > 5% N  
ped <- ped[which(ped$pct_N <= 0.05),]

# Filter out bad and uninformative markers
# Anything with almost all ref or alt is uninformative - remove markers where ref/alt is > 95% 
ped <- ped[which((ped$num_alt/num_F2s) <= 0.95),] 
ped <- ped[which((ped$num_ref/num_F2s) <= 0.95),] 
# Remove markers with almost all het calls (or all H/N calls)
ped <- ped[which(ped$het_all != 1),] 

if (crosstype == "f2"){
  # For all chromosomes (autosomes and X): looking for ref:alt of ~1:1
  ped <- ped[which((ped$ref_alt >= 0.4) & (ped$ref_alt <= 0.6)),] # Looking for ref:alt of ~1:1
  # For autosomes, looking for het_all to be ~0.5. For X-chr, het_all should be ~0.25
  ped <- ped[which(((ped$Chromosome != "X") & (ped$het_all >= 0.4) & (ped$het_all <= 0.6)) | 
                       ((ped$Chromosome == "X") & (ped$het_all >= 0.15) & (ped$het_all <= 0.35))),]
} else if (crosstype == "bc") {
  ped <- ped[which((ped$het_all >= 0.4) & (ped$het_all <= 0.6)),] # Looking for het:hom ratio of ~1:1
  ped <- ped[which((ped$ref_alt == 0) | (ped$ref_alt == 1)),] # Need hard cutoff to satisfy R/qtl 
}

# remove marker stats 
ped <- ped[,-c(376:382)]

#--------------------------------Formatting------------------------------------#
map <- ped # save for map file

rownames(ped) <- ped$Marker
ped <- ped[,-c(1:3)]
ped <- ped[,-c(3:7)] # change to 3:7
ped <- t(ped)

ped <- do.call(cbind, apply(ped, 2, function(x) data.frame(x,x)))

for (i in 3:nrow(ped)){ # for each row starting at row 3
  for (j in 1:ncol(ped)){ # for each column 
    if (ped[i,j] == 'H'){
      if (j%%2 == 1){ # odd row
        ped[i,j] <- ped[1,j] # assign reference allele 
      } else if (j%%2 == 0){ # even row 
        ped[i,j] <- ped[2,j] # assign alternate allele 
      }
    }
  }
}

ped <- ped[-c(1:2),] # remove ref/alt rows 
rownames(ped) <- toupper(rownames(ped)) # capitalize mouse IDs

#----------------------------------Add sex-------------------------------------#
ped <- ped %>% add_column('sex' = rep(NA, nrow(ped)), .before = 1)
 
### All females 
ped$sex <- rep(2, nrow(ped)) 
ped <- ped %>% add_column('mouseID' = rownames(ped), .before=1) # bc rownames aren't saved to file

ensure_directory("derived_data/grm-files")
write_delim(ped, col_names = FALSE, delim = " ",
            file = "derived_data/grm-files/pnut.ped")

#------------------------------Create map file---------------------------------#
map$Position_b38 <- map$Position_b38*1e6 # convert back to bp

map <- map[,c(2,1,3)] # only need chr, snp identifier, and bp position (in that order)

write_delim(map, col_names = FALSE, delim = " ", 
            file = "derived_data/grm-files/pnut.map")





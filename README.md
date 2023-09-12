# A mutation in *Themis* contributes to peanut-induced oral anaphylaxis in CC027 mice 

This document describes how to reproduce heritability and QTL analyses of various peanut allergy-related phenotypes in a backcross between CC0027 (reacts orally to peanut) and C3H (does not react orally to peanut). 

Environment prep
----------------

This project uses a Docker container to produce an environment similar to that used in the original analysis (e.g. R v4.2.1 and R package versions available on August 1, 2023). In order to run this container you will need [Docker](https://docs.docker.com/get-docker/) installed. 

Build the docker container:

```
docker build . -t pnut 
```

Run the docker container, opening a terminal session within the container:

```
docker run -e PASSWORD=pw123 --rm -v $(pwd):/home/rstudio/work -p 8787:8787 -it pnut /bin/bash
```

Navigate to the working directory: 

```
cd home/rstudio/work 
```

Prepare R/qtl file 
-----------------------

Run the following code to produce a file in the format required by R/qtl.

```
Rscript rqtl-file-gen-pnutbc.R
```
This script will generate the following .csv file, which can be imported into R and analyzed using the `Rqtl` package: `derived_data/Rqtl_CC27xC3H_BCv2.csv` 

Inbred analysis
-----------------------

Analyze inbred parent data:

```
Rscript parent-analysis.R
```

This script will generate:
* `figs/parent-temps.png`: plot of temperature trajectories in inbred CC027 and C3H mice (Fig. 1A in manuscript) 
* `results/inbred-parent-data-summary.csv`: means and ranges for each phenotype by strain (Supp. Table 1 in manuscript)
* `results/h2-from-inbred-parents.csv`: heritability estimates from inbred parent data (Supp. Table 1 in manuscript)

Backcross analysis 
-----------------------

**Note**: Add steps to produce genomic relationship matrix for backcross mice (will be used in heritability calculations). For now, plink files are included in `source_data` directory (`.map` and `.bed` produced by `plink_file_setup.R`, everything else produced by commands in `plink.sh`).

Generate a report of analysis performed on backcross mice, including heritability estimation and QTL mapping. 

```
Rscript -e 'library(rmarkdown); rmarkdown::render("qtl-analysis.Rmd", "html_document")'
```

This script will generate a report, `qtl-analysis.html`, as well as:
* several figures in the `figs` directory (including figures that compose Figures 2-4 and Supplemental Figures 1-5 in manuscript)
* `results/h2-from-bc-mice.csv`: heritability estimates from backcross mice (Supp. Table 1 in manuscript)
* `results/QTLsummary.csv`: summary of QTL mapping results (Table 1 in manuscript)


Candidate gene analysis
-----------------------

First, unzip the compressed directory of VCF files for each QTL region:

```
unzip source_data/VCFs.zip -d source_data/
```

Then unzip the GTF file representing all genes in the mouse genome:

```
gunzip source_data/Mus_musculus.GRCm38.97.gtf.gz --keep
```

Create a map between Ensembl IDs and gene names to be used by later scripts:

```
Rscript make_gene_map.R
```

Count total genes in each QTL region:

```
Rscript count_genes.R
```

Count *candidate* genes in each QTL region:  

```
bash cand_gene.sh
```

This analysis can take several hours. To run the above command on a SLURM-based high-performance computing cluster, run:

```
bash cand_gene.sh -m slurm 
```

This will generate two directories:
* `results/cand-gene-analysis/summary`: summary of candidate gene analysis for the region 
* `results/cand-gene-analysis/vardata`: data on candidate genes in the region, including the number of total variants, regulatory variants and protein-modifying variants 

Summarize candidate gene data in table:  

```
Rscript summarize_cand_genes.R
```

This script will produce a summary of the candidate gene analysis results, in `results/cand-gene-summary.csv` (Table 2 in manuscript). 


Themis analysis
-----------------------

This script will generate a report of analysis (linear regression, Tukey's honest significant difference, and PCA) on flow phenotypes from *Themis*-variant CC strains. 

```
Rscript -e 'library(rmarkdown); rmarkdown::render("themis-analysis.Rmd", "html_document")'
```

This script will generate a report, `themis-analysis.html`, along with several figures in `figs/fig5` (Fig. 5 in manuscript) and `figs/supplemental/tukey/` (Supp. Fig. 6 in manuscript).
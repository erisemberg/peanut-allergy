module load plink

plink --file pnut --make-bed --no-fid --no-parents --no-pheno --out pnut

module load gcta 

gcta64 --bfile pnut --make-grm --autosome-num 19 --autosome --out pnut

gcta64 --bfile pnut --make-grm-xchr --out pnut_xchr

gcta64 --grm pnut --pheno pheno_aac.phen --reml --out pnut_aac --thread-num 10 
#!/bin/bash
# Loop over QTL and perform candidate gene analysis  

while getopts "m:" OPTION
do 
	case "$OPTION" in 
		m) mode=$OPTARG
	esac
done

if [ -z "$mode" ]; then mode="bash"; fi # default mode is bash 

if [ "$mode" != "bash" ] && [ "$mode" != "slurm" ]; then
	echo "ERROR: -m flag must be set to \"bash\" or \"slurm\"" 
	exit 
fi 

mkdir -p results/cand-gene-analysis
mkdir -p results/cand-gene-analysis/summary
mkdir -p results/cand-gene-analysis/vardata
mkdir -p logs/cand-gene-analysis 

qtls=("Qpa1" "Qpa2" "Qpa3" "Qpa4" "Qpa5" "Qpa6" "Qpa7" "Qpa8")

for qtl in ${qtls[@]}; do
	vcf=$(find source_data/VCFs -name "$qtl*" | tr "\n" ",")

	if [ $mode == "slurm" ]; then 
		logfile="logs/${qtl}.out"
		sbatch --mem=10G -t 4:00:00 --output=${logfile} --wrap="module add r; Rscript cand_gene_analysis.R --args --vcf=${vcf} --prefix=${qtl}"
	fi
	if [ $mode == "bash" ]; then 
		logfile="logs/${qtl}.log"
		nohup Rscript cand_gene_analysis.R --args --vcf=${vcf} --prefix=${qtl} > ${logfile} 2>&1 < /dev/null &
	fi
done





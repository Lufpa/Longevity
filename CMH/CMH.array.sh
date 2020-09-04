#!/bin/bash
#SBATCH --mem=10000
#SBATCH --time=24:00:00 --qos=1day
#SBATCH --job-name=CMH_all
#SBATCH --cpus-per-task=1 
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --array=1-10 #adjust according to the nrow of the "listtables" file

## Requires: listtables tab separated file with 5 columns: altfile reffile treatment chr replicates(2reps,3reps)

## Requires: alt and ref tables with SNPs that overlap all cages.
## The tables HAVE to be organized in the following order:
## a1 a2 at01 at02 dt0 d1

## Requires: a tab delimited file with the column # of the first and last fly for each cage in the alt/ref tables
## The rows must correspond to the order of the actual alt and ref tables. e.g:
## a1	1	500
## a2	501	1000
## t01	1001	1500
## t02	1501	200
## filename should be "indpercage_arg[3]_arg[4].txt - arg3-4 are col3-4 from the listtables table

## This code downsamples read counts to 1 read count per individual per SNP to control for differential coverage per individual sample. 
## The downsampling is done 1000 times to get an good average estimate of the real ref/alt counts per SNP
## A file with the average allele count per SNP is generated and formated to run CMH test in PoPoolation2 
 
date
inDIR=/Genomics/ayroleslab2/shared/longevity/CMH/Permutation_indir
outDIR=/Genomics/ayroleslab2/shared/longevity/CMH/PoPoo_outdir
poppath=/Genomics/ayroleslab2/shared/longevity/CMH/popoolation2
perm=1000
list=listtables_realdata

source ~/scripts/Longevity/CMH.pipeline.sh  ${inDIR} ${outDIR} ${poppath} ${perm} ${list}

echo "Done!"
date

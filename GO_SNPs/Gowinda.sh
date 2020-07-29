#!/bin/bash
#SBATCH --mem=40000
#SBATCH --time=1:00:00 --qos=1hr
#SBATCH --job-name=GOl_gene
#SBATCH --cpus-per-task=8
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --array=1-4  #if running this for the 4 groups (shared SNPs, ctrl only, hs only, all) set 1-4

#Requires:
	# --snp-file: file with all SNPs used in the analysis, chromosome should not contain "chr"
	# --candidate-snp-file: file with the SNPs to be tested for enrichment, chrm should not contain "chr". Tto use with array, just keep track of what 1-4 means. Here, 1 = shared SNPs, 2 = ctrl only, 3 = hs only, 4 = all SNPs.
	# --gene-set-file: file with the GO terms downloaded from FuncAssociate2
	# --annotation-file: we're using the genome gtf version 6.23



date

##### Calculate enrichment assuming all SNPs are independent + extends the gene 1kb up and downstream
##### This fixes number of SNPs to resample, and as a result the total number of genes per permutation varies

#       java -Xmx4g -jar ~/bin/Gowinda-1.12.jar --snp-file TotalSNPs.txt --candidate-snp-file CandidateSNPs.${SLURM_ARRAY_TASK_ID} --gene-set-file funcassociate_go_associations.txt --annotation-file dmel-all-r6.23.gtf --simulations 1000000 --min-significance 1 --gene-definition updownstream1000 --threads 8 --output-file results_${SLURM_ARRAY_TASK_ID}_snp_updown1k.txt --mode SNP --min-genes 5


#### Calculate enrichment assuming all SNPs in a gene are haplogroup + extends the gene 1kb up and downstream 
#### This option fixes the number of genes to resample, and therefore the number of SNPs per permutation might vary

        java -Xmx4g -jar ~/bin/Gowinda-1.12.jar --snp-file TotalSNPs.txt --candidate-snp-file CandidateSNPs.${SLURM_ARRAY_TASK_ID} --gene-set-file funcassociate_go_associations.txt --annotation-file dmel-all-r6.23.gtf --simulations 100000 --min-significance 1 --gene-definition updownstream1000 --threads 8 --output-file results_${SLURM_ARRAY_TASK_ID}_gene_updown1k.txt --mode gene --min-genes 3


echo 'done!'

date


#/bin/bash

inDIR=$1
outDIR=$2
poppath=$3
perm=$4
list=$5

tableREF=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $1 }' $inDIR/$list`
tableALT=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $2 }' $inDIR/$list`
echo "tables to downsample  " $tableALT $tableREF >&2
treatment=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $3 }' $inDIR/$list`
chr=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $4 }' $inDIR/$list`
reps=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $5 }' $inDIR/$list`
part=`awk -v file=$SLURM_ARRAY_TASK_ID '{if (NR==file) print $6 }' $inDIR/$list`
echo "treatment and chr " $treatment $chr >&2

Rscript ~/scripts/Longevity/CMH.downsampling.prob.R $tableALT $tableREF $perm $treatment $chr $reps $inDIR $part

inputfile=syncfile.${treatment}.${chr}.${reps}.${part}.sync
echo $inputfile 
outputfile=${outDIR}/${inputfile%.sync}.cmh
echo $outputfile

if [ $reps == '2reps' ]
        then
                perl $poppath/cmh-test.pl --input $inputfile --output $outputfile --min-count 12 --min-coverage 50 --max-coverage 1000 --population 1-3,2-4
        else
perl $poppath/cmh-test.pl --input $inputfile --output $outputfile --min-count 12 --min-coverage 50 --max-coverage 1000 --population 1-3,2-4,5-6
fi



date

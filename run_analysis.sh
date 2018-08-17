#!/usr/bin/env bash

#################################################
# Check if script is run from the exomeseq folder
# If yes, move to the parent (main project) folder
##################################################
folder=${PWD##*/}
if [ ${folder} = 'exomeseq' ]
then
    cd ..
fi

#################################################
# ANALYSIS:
# Run the analysis for each sample:
# alignment, variant calling, annotation,
# output-parsing (after all samples are analysed)
#################################################
for SAMPLE in `cat samples.paired.list samples.single.list`;
do
    ${SAMPLE}/Scripts/${SAMPLE}_GATK.sh
done
exomeseq/parsevcf.py -o full_report.txt

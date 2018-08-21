#!/bin/bash

set -e

#################################################
# SETUP Part0:
# Check if script is run from the exomeseq folder
# If yes, move to the parent (main project) folder
##################################################
folder=${PWD##*/}
if [ ${folder} = 'exomeseq' ]
then
    cd ..
fi
######################################################
# SETUP Part1:
# Create folders and configuration file
# Sort and place fastq files, collect sample names
######################################################
exomeseq/setup.py

#############################################################
# SETUP Part2:
# Download reference and database files
# Populate the resources folder with the required files
#############################################################

# NB! This link is temporary and is not guaranteed to work outside of the time for the workshop, Sept 10-11.
wget https://si-medbif-exomeseq.sgp1.digitaloceanspaces.com/exome_pipeline_resources.tar.gz
tar xzvf exome_pipeline_resources.tar.gz

#cd resources/hg38bundle
#wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.*
#docker run --rm -v $PWD:/data biocontainers/bwa bwa index /data/Homo_sapiens_assembly38.fasta.gz
#gunzip -c Homo_sapiens_assembly38.fasta.gz > Homo_sapiens_assembly38.fasta
#wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz*
#wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz*
#cd ../snpeff_db
## dbNSFP
#wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.5a.zip
## Preparing the dbNSFP database: (from https://github.com/bcbio/bcbio-nextgen/issues/1571)
#unzip dbNSFPv3.5a.zip
#head -n1 dbNSFP3.5a_variant.chr1 > dbNSFPv3.5a.txt
#cat dbNSFP3.5a_variant.chr* | grep -v ^# >> dbNSFPv3.5a.txt
#bgzip dbNSFPv3.5a.txt
#tabix -s 1 -b 2 -e 2 dbNSFPv3.5a.txt.gz
## PhastCons
#wget --recursive --no-parent -l2 http://hgdownload-test.cse.ucsc.edu/goldenPath/hg38/phastCons100way/
#mv hgdownload-test.cse.ucsc.edu/goldenPath/hg38/phastCons100way/ phastCons100way
#rm -rf hgdownload-test.cse.ucsc.edu/
## ClinVar
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20180729.vcf.gz*
## GWASCat
#wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations.tsv
#cd ../..
##################################
# SETUP Part3:
# Download required docker images
##################################
docker pull biocontainers/fastqc
docker pull biocontainers/bwa
docker pull broadinstitute/picard
docker pull broadinstitute/gatk3:3.8-1
docker pull fjukstad/trimmomatic
docker build --rm -t "snpeff38:v1" exomeseq/snpeff4.3/.
#####################################
# SETUP Part4:
# Create runscripts for each sample
#####################################
for SAMPLE in `cat samples.paired.list`;
do
    exomeseq/make_scripts.py ${SAMPLE}
    sh ${SAMPLE}/Scripts/0_${SAMPLE}_fastqc.sh
done
#################################################
# ANALYSIS:
# Run the analysis for each sample:
# QC,
# alignment, variant calling, annotation,
# output-parsing (after all samples are analysed)
#################################################
for SAMPLE in `cat samples.paired.list samples.single.list`;
do
    ${SAMPLE}/Scripts/${SAMPLE}_GATK.sh
done
exomeseq/parsevcf.py -o full_report.txt


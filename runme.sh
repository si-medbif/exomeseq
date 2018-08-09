#!/bin/bash
# Create folders and configuration file
# Sort and place fastq files, collect sample names
# Download reference and database files
exomeseq/setup.py
# Populate the resources folder with the required references
cd resources/hg38bundle
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.*
docker run --rm -v $PWD:/data biocontainers/bwa bwa index /data/Homo_sapiens_assembly38.fasta.gz
unzip -c Homo_sapiens_assembly38.fasta.gz > Homo_sapiens_assembly38.fasta
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz*
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz*
cd ../snpeff_db
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.5c.zip
# Add more databases here according to need.
# wget phastcons
# wget clinvar
# wget gwascatalog
cd ../..
# Create runscripts for each sample
for SAMPLE in `cat samples.paired.list`;
do
    exomeseq/make_scripts.py ${SAMPLE}
done
# Run the analysis for each sample: QC, alignment, variant calling, annotation, output-parsing
for SAMPLE in `cat samples.paired.list`;
do
    exomeseq/${SAMPLE}_GATK.sh
    exomeseq/${SAMPLE}_annotate.sh
    exomeseq/parsevcf.py
done

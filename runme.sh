#!/bin/bash
# Create folders and configuration file
# Sort and place fastq files, collect sample names
# Download reference and database files
exomeseq/setup.py
###########
# Populate the resources folder with the required references
###########
cd resources/hg38bundle
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.*
docker run --rm -v $PWD:/data biocontainers/bwa bwa index /data/Homo_sapiens_assembly38.fasta.gz
gunzip -c Homo_sapiens_assembly38.fasta.gz > Homo_sapiens_assembly38.fasta
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz*
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz*
cd ../snpeff_db
# Download databases for annotation
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.5a.zip
# Preparing the dbNSFP database: (from https://github.com/bcbio/bcbio-nextgen/issues/1571)
unzip dbNSFPv3.5a.zip
head -n1 dbNSFP3.5a_variant.chr1 > dbNSFPv3.5a.txt
cat dbNSFP3.5a_variant.chr* | grep -v ^# >> dbNSFPv3.5a.txt
bgzip dbNSFPv3.5a.txt > tabix -s 1 -b 2 -e 2 dbNSFPv3.5a.txt.gz
# More databases
wget --recursive --no-parent -l2 http://hgdownload-test.cse.ucsc.edu/goldenPath/hg38/phastCons100way/
#*Should not be necessary if the above command works properly.
mv hgdownload-test.cse.ucsc.edu/goldenPath/hg38/phastCons100way/ phastCons100way
rm -rf hgdownload-test.cse.ucsc.edu/
#****
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20180729.vcf.gz*
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations.tsv
cd ../..
#########
# Get docker images
#########
#docker pull biocontainers/fastqc
#docker pull biocontainers/bwa
#docker pull broadinstitute/picard
#docker pull broadinstitute/gatk3:3.8-1
#docker pull biocontainers/snpeff
#########
# Create runscripts for each sample
#########
for SAMPLE in `cat samples.paired.list`;
do
    exomeseq/make_scripts.py ${SAMPLE}
done
##########
# Run the analysis for each sample: QC, alignment, variant calling, annotation, output-parsing
##########
for SAMPLE in `cat samples.paired.list`;
do
    ./${SAMPLE}_GATK.sh
    ./${SAMPLE}_annotate.sh
    exomeseq/parsevcf.py
done

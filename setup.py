#!/usr/bin/env python3
"""
Module Docstring
"""

__author__ = "Harald Grove"
__version__ = "1.0.0"
__license__ = "MIT"

import argparse
import os
import sys
import shutil
import subprocess
import glob


def do_setup(args):
    """
    This file should be located in the exomeseq sub-folder to the project folder
    :param args:
    :return:
    """
    try:
        os.mkdir("{}/fastq_paired".format(args.name))
        os.chmod("{}/fastq_paired".format(args.name), 0o777)
        os.mkdir("{}/fastq_single".format(args.name))
        os.chmod("{}/fastq_single".format(args.name), 0o777)
        os.mkdir("{}/resources".format(args.name))
        os.chmod("{}/resources".format(args.name), 0o777)
        os.mkdir("{}/resources/hg38bundle".format(args.name))
        os.chmod("{}/resources/hg38bundle".format(args.name), 0o777)
        os.mkdir("{}/resources/snpeff_db".format(args.name))
        os.chmod("{}/resources/snpeff_db".format(args.name), 0o777)
    except FileExistsError:
        sys.stderr.write("WARNING! Project folder already exists...\n")
        # While testing, just keep this as a warning.
        #sys.exit(1)


def make_cfg(args):
    """ Create the configuration file with location of required files and folders """
    cfg_file = "{}/exome.cfg".format(args.name)
    with open(cfg_file, "w") as fout:
        fout.write("##################\n")
        fout.write("# Folders containing input and output files\n")
        fout.write("##################\n")
        fout.write("out_dir={}\n".format(args.name[1:]))
        fout.write("fastq_paired_dir=fastq_paired\n")
        fout.write("fastq_single_dir=fastq_single\n")
        fout.write("##################\n")
        fout.write("# Project specific reference files\n")
        fout.write("##################\n")
        fout.write("exon_bed={}\n".format(args.regions38))
        # These are files with SNP allele frequencies, will not be universally applicable
        if args.afdb:
            fout.write("freq_main=resources/allelefreqs/1KG/variants_hg38.frq\n")
            fout.write("exac_freqfile=resources/allelefreqs/ExAC/ExAC_exome_SNV_freq.txt\n")
            fout.write("hgvd_freqfile=resources/allelefreqs/HGVD/HGVD_freq.txt\n")
            fout.write("esp6500_freqfile=resources/allelefreqs/ESP6500/ESP6500_freq.txt\n")
            fout.write("gonl_freqfile=resources/allelefreqs/GoNL/gonl_freq.txt\n")
            fout.write("clinvar_freq=resources/allelefreqs/clinvar/clinvar.vcf\n")
            fout.write("mutationtaster=resources/allelefreqs/mutationtaster/mutationtaster.list\n")
        fout.write("##################\n")
        fout.write("# Generic reference files\n")
        fout.write("##################\n")
        fout.write("ref_dir=resources/hg38bundle\n")
        fout.write("ref_genome=resources/hg38bundle/Homo_sapiens_assembly38.fasta\n")
        fout.write("indel_1=resources/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz\n")
        fout.write("DBSNP=resources/hg38bundle/dbsnp_146.hg38.vcf.gz\n")
        fout.write("##################\n")
        fout.write("# Parameter settings\n")
        fout.write("##################\n")
        fout.write("cores=4\n")
        fout.write("gatk_num_threads=1\n")
        fout.write("gatk_num_cpu_threads=4\n")
        fout.write("java_mem=7G\n")
        fout.write("##################\n")
        fout.write("# SNPEff databases\n")
        fout.write("##################\n")
        fout.write("snpeff_dbver=GRCh38.86\n")
        fout.write("DBNSFP=resources/snpeff_db/dbNSFPv3.5a.txt.gz\n")
        fout.write("GWASCATALOG=resources/snpeff_db/gwas-catalog-associations.tsv\n")
        fout.write("PHASTCONS=resources/snpeff_db/phastCons100way/hg38.100way.phastCons\n")
        fout.write("CLINVAR=resources/snpeff_db/clinvar_20180729.vcf.gz\n")
        fout.write("##################\n")
        fout.write("# Docker images\n")
        fout.write("##################\n")
        fout.write("FASTQC=biocontainers/fastqc\n")
        fout.write("BWA=biocontainers/bwa\n")
        fout.write("PICARD=broadinstitute/picard\n")
        fout.write("GATK=broadinstitute/gatk3:3.8-1\n")
        fout.write("TRIMMOMATIC=fjukstad/trimmomatic\n")
        fout.write("SNPEFF=snpeff38:v1\n")
        fout.write("##################\n")
        fout.write("# Trimmomatic settings\n")
        fout.write("##################\n")
        fout.write("TRIM_window=4\n")
        fout.write("TRIM_score=15\n")
        fout.write("TRIM_minlen=36\n")

def process_fastq(args):
    """ Collects fastq files and moves them to either pared or single folder
        Also looks for a single bed file specifying the target regions
    """
    file_list = glob.glob("*.*")
    samples = {}
    count_bed = 0
    args.regions38 = 'NA'
    for filename in file_list:
        # names of fastq files are expected to match "*.[12].fastq.gz"
        lf = filename.rsplit('.',3)
        if lf[-1] == 'bed':
            count_bed += 1
            continue
        elif len(lf) < 4 or lf[1] not in ['0','1','2'] or lf[2] not in ['fastq'] or lf[3] not in ['gz']:
            sys.stdout.write('WARNING: Unexpected file {}.\n'.format(filename))
            continue
        samples[lf[0]] = samples.get(lf[0], 0) + 1
    with open('samples.paired.list', 'w') as fout2, open('samples.single.list', 'w') as fout1:
        for filename in file_list:
            lf = filename.rsplit('.', 3)
            if lf[-1] == 'bed' and count_bed == 1:
                sys.stdout.write('Found target regions bed file: {}.\n'.format(filename))
                args.regions38 = filename
                continue
            elif len(lf) < 4 or lf[1] not in ['0', '1', '2'] or lf[2] not in ['fastq'] or lf[3] not in ['gz']:
                continue
            if lf[1] in ['1','2']:
                if samples[lf[0]] == 2:
                    shutil.move(filename,"fastq_paired/")
                    if lf[1] == '1':
                        fout2.write('{}\n'.format(lf[0]))
                    continue
                else:
                    sys.stderr.write(
                        'ERROR: File {} is indicated as one of a pair but the matching file is missing.\n'.format(filename))
                    continue
            if lf[1] in ['0']:
                if samples[lf[0]] == 1:
                    shutil.move(filename,"fastq_single/")
                    fout1.write('{}\n'.format(lf[0]))
                else:
                    sys.stderr.write(
                        'ERROR: File {} is indicated as single, but there is another file with same sample name.\n'.format(filename))
                    continue

def main(args):
    # Determine current folder, this should be the top project folder
    # Everything else should go below this
    # Adjust location if this script is run from the 'exomeseq' folder
    args.name = os.getcwd()
    ln = args.name.split('/')
    if ln[-1] == 'exomeseq':
        args.name = '/'.join(ln[0:-1])
    # Create folder structure, download references, make configuration
    # file and populate it with the available references.
    do_setup(args)
    process_fastq(args)
    make_cfg(args)


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--afdb", action="store_true", default=False, help="Include allele frequency databases.")
    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v", "--verbose", action="count", default=0, help="Verbosity (-v, -vv, etc)"
    )

    # Specify output of '--version'
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    args = parser.parse_args()
    main(args)

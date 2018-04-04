#!/usr/bin/env python3
"""
Module Docstring
"""

__author__ = "Harald Grove"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import os
import sys
import shutil

def do_setup(args):
    try:
        os.mkdir('{}'.format(args.name))
        os.chmod('{}'.format(args.name),0o777)
        os.mkdir('{}/fastq'.format(args.name))
        os.chmod('{}/fastq'.format(args.name),0o777)
    except FileExistsError:
        sys.stderr.write('ERROR! Project folder already exists...\n')
        # While testing, just keep this as a warning.
        sys.exit(1)
    shutil.copy2('make_scripts.py', args.name)
    shutil.copy2('parsevcf.py', args.name)
    shutil.copy2('vcf_header.txt',args.name)
    gene_file = '{}/genes.list'.format(args.name)
    with open(gene_file, 'w') as fout:
        fout.write('# List all genes that should be included in the report (if different from all).\n')
    sample_file = '{}/samples.list'.format(args.name)
    with open(sample_file, 'w') as fout:
        fout.write('# List of samples that will be included in the generated report\n')

def check_args(args):
    if not args.name.startswith('/'):
        sys.stderr.write('ERROR! Invalid path\n')
        sys.stderr.write('Please provide full path, starting with root folder.\n')
        sys.exit(1)
    if args.fastq is not None and not args.fastq.startswith('/'):
        sys.stderr.write('Fastq folder must be located somewhere under root ("/"), and must contain the full path.\n')
        sys.exit(1)

def make_cfg(args):
    """ Create the configuration file with location of required files and folders """
    cfg_file = '{}/exome.cfg'.format(args.name)
    with open(cfg_file, 'w') as fout:
        fout.write('# Folders containing input and output files\n')
        if args.fastq is not None:
            fout.write('fastq_dir={}/fastq\n'.format(args.fastq[1:]))
        else:
            fout.write('fastq_dir={}/fastq\n'.format(args.name[1:]))
        fout.write('out_dir={}\n'.format(args.name[1:]))
        fout.write('# Project specific reference files\n')
        if args.regions38 is not None:
            fout.write('exon_bed={}\n'.format(args.regions38.strip('/')))
        else:
            fout.write('exon_bed=NA\n')
        fout.write('freq_main=tiger/resources/allelefreqs/1KG/variants_hg38.frq\n')
        fout.write('exac_freqfile=tiger/resources/allelefreqs/ExAC/ExAC_exome_SNV_freq.txt\n')
        fout.write('hgvd_freqfile=tiger/resources/allelefreqs/HGVD/HGVD_freq.txt\n')
        fout.write('esp6500_freqfile=tiger/resources/allelefreqs/ESP6500/ESP6500_freq.txt\n')
        fout.write('gonl_freqfile=tiger/resources/allelefreqs/GoNL/gonl_freq.txt\n')
        fout.write('clinvar_freq1=tiger/resources/allelefreqs/clinvar/clinvar.vcf\n')
        fout.write('mutationtaster=tiger/resources/allelefreqs/mutationtaster/mutationtaster.list\n')
        fout.write('# Generic reference files\n')
        fout.write('ref_dir=tiger/resources/hg38bundle\n')
        fout.write('ref_genome=tiger/resources/hg38bundle/Homo_sapiens_assembly38.fasta\n')
        fout.write('indel_1=tiger/resources/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz\n')
        fout.write('indel_2=tiger/resources/hg38bundle/Homo_sapiens_assembly38.known_indels.vcf.gz\n')
        fout.write('DBSNP=tiger/resources/hg38bundle/dbsnp_144.hg38.vcf.gz\n')
        fout.write('# Parameter settings\n')
        fout.write('cores=4\n')
        fout.write('gatk_num_threads=1\n')
        fout.write('gatk_num_cpu_threads=4\n')
        fout.write('java_mem=8G\n')
        fout.write('# DO NOT CHANGE! The following directories refer to locations within the docker image for snpeff\n')
        fout.write('snpeff_dir=/home/snpeff/snpEff\n')
        fout.write('DBNSFP=/home/snpeff/snpEff/data/dbNSFP3.4a.txt.gz\n')
        fout.write('GWASCATALOG=/home/snpeff/snpEff/data/gwas_catalog_v1.0.1-associations_e84_r2016-07-10.tsv\n')
        fout.write('PHASTCONS=/home/snpeff/snpEff/data/phastCons100way\n')
        fout.write('CLINVAR=/home/snpeff/snpEff/data/clinvar.vcf\n')

def main(args):
    check_args(args)
    do_setup(args)
    make_cfg(args)

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("name", action="store", help="Main folder for project (full path)")

    # Optional argument flag which defaults to False

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-f", "--fastq", help="Folder containing fastq files, if files are to be kept separate.")
    parser.add_argument("-r", "--regions38", help="Bed file with hg38 region information (full path).")
    #parser.add_argument("-s", "--regions19", help="Bed file with hg19 region information.")

    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help="Verbosity (-v, -vv, etc)")

    # Specify output of '--version'
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s (version {version})'.format(version=__version__))

    args = parser.parse_args()
    main(args)

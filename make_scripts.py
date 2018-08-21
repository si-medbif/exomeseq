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


def read_config(args):
    db = {}
    with open(args.config, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            key, value = line.strip().split("=")
            db[key] = value
        if "exon_bed" in db and not db["exon_bed"].endswith('NA'):
            db["bed_argument"] = "-L /data/{}".format(db["exon_bed"])
        else:
            db["bed_argument"] = ""
    return db


def do_setup(args, db):
    try:
        os.mkdir("{}".format(args.name))
        os.chmod("{}".format(args.name), 0o777)
        subdirs = [
            "Scripts",
            "LOG",
            "TEMP",
            "SAM",
            "BAM",
            "BQSR",
            "GVCF",
            "VCF",
            "FastQC_pre",
            "FastQC_post",
            "QC",
            "QC/FILTERED",
            "Report",
        ]
        for sd in subdirs:
            os.mkdir("{}/{}".format(args.name, sd))
            os.chmod("{}/{}".format(args.name, sd), 0o777)
    except FileExistsError:
        sys.stderr.write("WARNING! Output folders already exists for this sample...\n")
        # sys.exit(1)


def make_preQC(args, db):
    """ Creates the script for running FastQC after setup but before any other processing"""
    script_file = "/{}/{}/Scripts/0_{}_fastqc.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step0: FastQC on raw data\n")
        fout.write("##-------------\n")
        fout.write("docker run --rm -v {}:/data {} fastqc ".format(db["out_dir"], db["FASTQC"]))
        fout.write("-o /data/{}/FastQC_pre ".format(args.name, args.name))
        fout.write("-t 2 ")
        if not args.single:
            fout.write("/data/{}/{}.1.fastq.gz ".format(db["fastq_paired_dir"], args.name))
            fout.write("/data/{}/{}.2.fastq.gz\n".format(db["fastq_paired_dir"], args.name))
            fout.write(
                "unzip -j -d /{}/{}/FastQC_pre/{}.1_fastqc /{}/{}/FastQC_pre/{}.1_fastqc.zip {}.1_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name, args.name
                )
            )
            fout.write(
                "unzip -j -d /{}/{}/FastQC_pre/{}.2_fastqc /{}/{}/FastQC_pre/{}.2_fastqc.zip {}.2_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name, args.name
                )
            )
            fout.write("# Check FastQC output files\n")
            fout.write("python /{}/exomeseq/collect_fastqc_data.py -o /{}/{}/Report/fastqc_pre_report.txt /{}/{}/FastQC_pre/{}.1_fastqc/fastqc_data.txt /{}/{}/FastQC_pre/{}.2_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], db["out_dir"], args.name,
                    db["out_dir"], args.name, args.name,
                    db["out_dir"], args.name, args.name
                )
            )
        else:
            fout.write("/data/{}/{}.0.fastq.gz\n".format(db["fastq_single_dir"], args.name))
            fout.write(
                "unzip -j -d /{}/{}/FastQC_pre/{}.0_fastqc /{}/{}/FastQC_pre/{}.0_fastqc.zip {}.0_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name, args.name
                )
            )
            fout.write("# Check FastQC output files\n")
            fout.write(
                "python /{}/exomeseq/collect_fastqc_data.py -o /{}/{}/Report/fastqc_pre_report.txt /{}/{}/FastQC_pre/{}.0_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], db["out_dir"], args.name,
                    db["out_dir"], args.name, args.name,
                    db["out_dir"], args.name, args.name
                )
            )


def make_postQC(args, db):
    """ Creates the script for running FastQC after optional Trimmomatic"""
    script_file = "/{}/{}/Scripts/0b_{}_fastqc.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step0b: FastQC on cleaned data\n")
        fout.write("##-------------\n")
        fout.write("docker run --rm -v {}:/data {} fastqc ".format(db["out_dir"], db["FASTQC"]))
        fout.write("-o /data/{}/FastQC_post ".format(args.name, args.name))
        fout.write("-t 2 ")
        if not args.single:
            fout.write("/data/{}/{}.1.fastq.gz ".format(db["fastq_paired_dir"], args.name))
            fout.write("/data/{}/{}.2.fastq.gz\n".format(db["fastq_paired_dir"], args.name))
            fout.write(
                "unzip -j -d /{}/{}/FastQC_post/{}.1_fastqc /{}/{}/FastQC_post/{}.1_fastqc.zip {}.1_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name, args.name
                )
            )
            fout.write(
                "unzip -j -d /{}/{}/FastQC_post/{}.2_fastqc /{}/{}/FastQC_post/{}.2_fastqc.zip {}.2_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name, args.name
                )
            )
            fout.write("# Check FastQC output files\n")
            fout.write("python /{}/exomeseq/collect_fastqc_data.py -o /{}/{}/Report/fastqc_post_report.txt /{}/{}/FastQC_post/{}.1_fastqc/fastqc_data.txt /{}/{}/FastQC_post/{}.2_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], db["out_dir"], args.name,
                    db["out_dir"], args.name, args.name,
                    db["out_dir"], args.name, args.name
                )
            )
        else:
            fout.write("/data/{}/{}.0.fastq.gz\n".format(db["fastq_single_dir"], args.name))
            fout.write(
                "unzip -j -d /{}/{}/FastQC_post/{}.0_fastqc /{}/{}/FastQC_post/{}.0_fastqc.zip {}.0_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name, args.name
                )
            )
            fout.write("# Check FastQC output files\n")
            fout.write(
                "python /{}/exomeseq/collect_fastqc_data.py -o /{}/{}/Report/fastqc_post_report.txt /{}/{}/FastQC_post/{}.0_fastqc/fastqc_data.txt\n".format(
                    db["out_dir"], db["out_dir"], args.name,
                    db["out_dir"], args.name, args.name,
                    db["out_dir"], args.name, args.name
                )
            )


def make_trimfastq(args, db):
    """ Creates the script for running Trimmomatic """
    script_file = "/{}/{}/Scripts/0a_{}_trimfastq.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step0a: Trimmomatic\n")
        fout.write("##-------------\n")
        fout.write("docker run --rm -v {}:/data {} ".format(db["out_dir"], db["TRIMMOMATIC"]))
        if not args.single:
            fout.write("PE ")
        else:
            fout.write("SE ")
        fout.write("-phred33 ")
        fout.write("-threads {} ".format(db["cores"]))
        fout.write("-trimlog /data/{}/LOG/0a_{}_trimfastq.log ".format(args.name,args.name))
        if not args.single:
            fout.write("/data/{}/{}.1.fastq.gz ".format(db["fastq_paired_dir"], args.name))
            fout.write("/data/{}/{}.2.fastq.gz ".format(db["fastq_paired_dir"], args.name))
            fout.write("/data/{}/{}_paired.1.fastq.gz ".format(db["fastq_paired_dir"], args.name))
            fout.write("/data/{}/{}_unpaired.1.fastq.gz ".format(db["fastq_paired_dir"], args.name))
            fout.write("/data/{}/{}_paired.2.fastq.gz ".format(db["fastq_paired_dir"], args.name))
            fout.write("/data/{}/{}_unpaired.2.fastq.gz ".format(db["fastq_paired_dir"], args.name))
        else:
            fout.write("/data/{}/{}.0.fastq.gz ".format(db["fastq_single_dir"], args.name))
            fout.write("/data/{}/{}_single.0.fastq.gz ".format(db["fastq_single_dir"], args.name))
        fout.write("LEADING:3 ")
        fout.write("TRAILING:3 ")
        fout.write("SLIDINGWINDOW:{}:{} ".format(db["TRIM_window"], db["TRIM_score"]))
        fout.write("MINLEN:{}\n".format(db["TRIM_minlen"]))
        if not args.single:
            fout.write("mv /{}/{}.1.fastq.gz /{}/{}.RAW.1.fastq.gz\n".format(db["fastq_paired_dir"], args.name, db["fastq_paired_dir"], args.name))
            fout.write("mv /{}/{}.2.fastq.gz /{}/{}.RAW.2.fastq.gz\n".format(db["fastq_paired_dir"], args.name, db["fastq_paired_dir"], args.name))
            fout.write("mv /{}/{}_paired.1.fastq.gz /{}/{}.1.fastq.gz\n".format(db["fastq_paired_dir"], args.name, db["fastq_paired_dir"], args.name))
            fout.write("mv /{}/{}_paired.2.fastq.gz /{}/{}.2.fastq.gz\n".format(db["fastq_paired_dir"], args.name,
                                                                                            db["fastq_paired_dir"],
                                                                                            args.name)
                       )
        else:
            fout.write("mv /{}/{}.0.fastq.gz /{}/{}.RAW.0.fastq.gz\n".format(db["fastq_single_dir"], args.name, db["fastq_single_dir"], args.name))
            fout.write("mv /{}/{}_single.0.fastq.gz /{}/{}.0.fastq.gz\n".format(db["fastq_single_dir"], args.name, db["fastq_single_dir"], args.name))

def make_align(args, db):
    """ Creates the script for aligning reads """
    script_file = "/{}/{}/Scripts/1_{}_align.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step1: Align\n")
        fout.write("##-------------\n")
        fout.write("docker run --rm -v {}:/data {} bwa mem ".format(db["out_dir"], db["BWA"]))
        fout.write("-t {} ".format(db["cores"]))
        fout.write(
            '-R "@RG\\tID:DM_{}\\tSM:{}\\tPL:Illumina\\tLB:WES\\tPU:unit1" '.format(
                args.name, args.name
            )
        )
        fout.write("/data/{}.gz ".format(db["ref_genome"]))
        if not args.single:
            fout.write("/data/{}/{}.1.fastq.gz ".format(db["fastq_paired_dir"], args.name))
            fout.write("/data/{}/{}.2.fastq.gz ".format(db["fastq_paired_dir"], args.name))
        else:
            fout.write("/data/{}/{}.0.fastq.gz ".format(db["fastq_single_dir"], args.name))
        fout.write(
            "> /{}/{}/SAM/{}_aligned.sam\n".format(db["out_dir"], args.name, args.name)
        )


def make_sort(args, db):
    """ Creates the script for sorting the aligned SAM file """
    script_file = "/{}/{}/Scripts/2_{}_sort.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step2: Sort\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} ".format(db["out_dir"], db["PICARD"])
        )
        fout.write("SortSam ")
        fout.write("INPUT=/data/{}/SAM/{}_aligned.sam ".format(args.name, args.name))
        fout.write("OUTPUT=/data/{}/BAM/{}_sorted.bam ".format(args.name, args.name))
        fout.write("SORT_ORDER=coordinate ")
        fout.write("TMP_DIR=/tmp\n")


def make_deduplicate(args, db):
    """ Creates the script for removing duplicate reads """
    script_file = "/{}/{}/Scripts/3_{}_mark_duplicates.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step3: Mark duplicates\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} ".format(db["out_dir"], db["PICARD"])
        )
        fout.write("MarkDuplicates ")
        fout.write("INPUT=/data/{}/BAM/{}_sorted.bam ".format(args.name, args.name))
        fout.write("OUTPUT=/data/{}/BAM/{}_deduplicated.bam ".format(args.name, args.name))
        fout.write("METRICS_FILE=/data/{}/BAM/{}_deduplication_metrics.txt ".format(args.name, args.name))
        fout.write("CREATE_INDEX=TRUE\n")
        fout.write("TMP_DIR=/tmp\n")


def make_index(args, db):
    """
    Create the script for building the index
    Unnecessary, picard can generate index as part of analysis.
    """
    script_file = "/{}/{}/Scripts/4_{}_build_index.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step4: Build Index\n")
        fout.write(
            "docker run --rm -v {}:/data {} ".format(db["out_dir"], db["PICARD"])
        )
        fout.write("BuildBamIndex ")
        fout.write("INPUT=/data/{}/BAM/{}_sorted.bam ".format(args.name, args.name))
        fout.write("TMP_DIR=/tmp\n")

def make_realign(args, db):
    """ Create the script for performing Indel Realingment """
    script_file = "/{}/{}/Scripts/5_{}_realign_indels.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step5-1: Create aligner target\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T RealignerTargetCreator ")
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        try:
            fout.write("-known /data/{} ".format(db["indel_1"]))
        except KeyError:
            pass
        try:
            fout.write("-known /data/{} ".format(db["indel_2"]))
        except KeyError:
            pass
        fout.write("{} ".format(db["bed_argument"]))
        fout.write(
            "-I /data/{}/BAM/{}_deduplicated.bam ".format(
                args.name, args.name
            )
        )
        fout.write("-nt {} ".format(db["gatk_num_threads"]))
        fout.write("-dt NONE ")
        fout.write(
            "-o /data/{}/BAM/{}_indel_target_intervals.list ".format(
                args.name, args.name
            )
        )
        fout.write(
            "-log /data/{}/LOG/5-1_{}_indel_target_intervals.log ".format(
                args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write("##-------------\n")
        fout.write("##Step5-2: Realign indels\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T IndelRealigner ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write(
            "-I /data/{}/BAM/{}_deduplicated.bam ".format(
                args.name, args.name
            )
        )
        fout.write(
            "-targetIntervals /data/{}/BAM/{}_indel_target_intervals.list ".format(
                args.name, args.name
            )
        )
        try:
            fout.write("-known /data/{} ".format(db["indel_1"]))
        except KeyError:
            pass
        try:
            fout.write("-known /data/{} ".format(db["indel_2"]))
        except KeyError:
            pass
        fout.write("-dt NONE ")
        fout.write("-o /data/{}/BAM/{}_realigned.bam ".format(args.name, args.name))
        fout.write("-log /data/{}/LOG/5-2_{}_realign_indels.log\n".format(args.name, args.name))

def make_BQSR(args, db):
    """ Create the script for performing Base Quality Score Recalibration """
    script_file = "/{}/{}/Scripts/6_{}_recalibrate_base.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step6-1: Perform Base Recalibration\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T BaseRecalibrator ")
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        try:
            fout.write("-knownSites /data/{} ".format(db["indel_1"]))
        except KeyError:
            pass
        try:
            fout.write("-knownSites /data/{} ".format(db["indel_2"]))
        except KeyError:
            pass
        try:
            fout.write("-knownSites /data/{} ".format(db["DBSNP"]))
        except KeyError:
            pass
        fout.write("{} ".format(db["bed_argument"]))
        fout.write("--interval_padding 100 ")
        fout.write("-I /data/{}/BAM/{}_realigned.bam ".format(args.name, args.name))
        fout.write("-nct {} ".format(db["gatk_num_cpu_threads"]))
        fout.write("-o /data/{}/BQSR/{}_perform_bqsr.table ".format(args.name, args.name))
        fout.write("-log /data/{}/LOG/6-1_{}_perform_bqsr.log ".format(args.name, args.name))
        # This line is needed when the FASTQ file contains the 'other' quality score format: Q64/phred64
        fout.write("#--fix_misencoded_quality_scores")
        fout.write("\n\n")
        fout.write("##-------------\n")
        fout.write("##Step6-4: Print Reads\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T PrintReads ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write("-I /data/{}/BAM/{}_realigned.bam ".format(args.name, args.name))
        fout.write("-BQSR /data/{}/BQSR/{}_perform_bqsr.table ".format(args.name, args.name))
        fout.write("-dt NONE ")
        fout.write("-EOQ ")
        fout.write("-nct {} ".format(db["gatk_num_cpu_threads"]))
        fout.write("-o /data/{}/BAM/{}_GATK.bam ".format(args.name, args.name))
        fout.write("-log /data/{}/LOG/6-4_{}_perform_bqsr.log\n".format(args.name, args.name))


def make_call_haplotype(args, db):
    """ Create script for calling haplotypes """
    script_file = "/{}/{}/Scripts/7_{}_call_haplotype.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step7: Call Haplotype\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T HaplotypeCaller ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        fout.write(
            "--input_file /data/{}/BAM/{}_GATK.bam ".format(args.name, args.name))
        fout.write("--emitRefConfidence GVCF ")
        fout.write("--genotyping_mode DISCOVERY ")
        fout.write("{} ".format(db["bed_argument"]))
        fout.write("--interval_padding 100 ")
        fout.write("-o /data/{}/GVCF/{}_GATK.g.vcf ".format(args.name, args.name))
        fout.write("-log /data/{}/LOG/7_{}_haplotype_caller.log ".format(args.name, args.name))
        fout.write("-A DepthPerSampleHC ")
        fout.write("-pairHMM VECTOR_LOGLESS_CACHING ")
        fout.write("-nct {}\n".format(db["gatk_num_cpu_threads"]))


def make_genotype(args, db):
    """ Create script for calling genotypes """
    script_file = "/{}/{}/Scripts/8_{}_genotype_gvcf.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step8: Genotype\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T GenotypeGVCFs ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        fout.write("--variant /data/{}/GVCF/{}_GATK.g.vcf ".format(db["out_dir"], args.name, args.name))
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write("-nt 1 ")
        fout.write("-o /data/{}/VCF/{}_RAW.vcf ".format(args.name, args.name)
        )
        fout.write(
            "-log /data/{}/LOG/8_{}_genotype_gvcf.log".format(args.name, args.name))
        fout.write("\n\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T VariantAnnotator ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        fout.write("--variant /data/{}/VCF/{}_RAW.vcf ".format(args.name, args.name))
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write("-I /data/{}/BAM/{}_GATK.bam ".format(args.name, args.name))
        fout.write("--dbsnp /data/{} ".format(db["DBSNP"]))
        fout.write("-L /data/{}/VCF/{}_RAW.vcf ".format(args.name, args.name))
        fout.write("-dt NONE ")
        fout.write("-nt 1 ")
        fout.write("-o /data/{}/VCF/{}_RAW_ANNOTATED.vcf ".format(args.name, args.name))
        fout.write("-log /data/{}/LOG/8-1_{}_QC_snv_annotation.log ".format(args.name, args.name))
        fout.write("-A ClippingRankSumTest ")
        fout.write("-A ReadPosRankSumTest ")
        fout.write("-A MappingQualityRankSumTest ")
        fout.write("-A GCContent ")
        fout.write("-A AlleleBalanceBySample ")
        fout.write("-A AlleleBalance ")
        fout.write("-A VariantType\n")


def make_SNV_QC(args, db):
    """ Create the script for filtering and QC of SNVs """
    script_file = "/{}/{}/Scripts/9_{}_SNV_quality_control.sh".format(
        db["out_dir"], args.name, args.name
    )
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##Step9-1-1: Extract SNPs\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T SelectVariants ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        fout.write("--variant /data/{}/VCF/{}_RAW_ANNOTATED.vcf ".format(args.name, args.name))
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write("-selectType SNP ")
        fout.write("--excludeFiltered ")
        fout.write("-nt 1 ")
        fout.write("-o /data/{}/VCF/{}_RAW_SNV.vcf ".format(args.name, args.name))
        fout.write("-log /data/{}/LOG/9-1-1_{}_QC_select_snv.log".format(args.name, args.name))
        fout.write("\n\n")
        fout.write("##-------------\n")
        fout.write("##Step9-1-2: Filter SNPs\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T VariantFiltration ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        fout.write("--variant /data/{}/VCF/{}_RAW_SNV.vcf ".format(args.name, args.name))
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write('--filterExpression "QD < 2.0" ')
        fout.write('--filterName "QD" ')
        fout.write('--filterExpression "DP < 8.0" ')
        fout.write('--filterName "DP" ')
        fout.write("--logging_level ERROR ")
        fout.write("-o /data/{}/VCF/{}_FILTERED_SNV.vcf ".format(args.name, args.name))
        fout.write("-log /data/{}/LOG/9-1-2_{}_QC_filter_snv.log\n".format(args.name, args.name))
        fout.write("\n")
        fout.write("##-------------\n")
        fout.write("##Step9-2-1: Extract INDELs\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T SelectVariants ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        fout.write("--variant /data/{}/VCF/{}_RAW_ANNOTATED.vcf ".format(args.name, args.name))
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write("-selectType INDEL ")
        fout.write("-selectType MNP ")
        fout.write("-selectType MIXED ")
        fout.write("-selectType SYMBOLIC ")
        fout.write("--excludeFiltered ")
        fout.write("-nt 1 ")
        fout.write("-o /data/{}/VCF/{}_RAW_INDEL.vcf ".format(args.name, args.name))
        fout.write("-log /data/{}/LOG/9-2-1_{}_QC_select_INDEL.log\n".format(args.name, args.name))
        fout.write("\n")
        fout.write("##-------------\n")
        fout.write("##Step9-2-2: Filter INDELs\n")
        fout.write("##-------------\n")
        fout.write(
            "docker run --rm -v {}:/data {} java -jar GenomeAnalysisTK.jar ".format(db["out_dir"], db["GATK"])
        )
        fout.write("-T VariantFiltration ")
        fout.write("-R /data/{} ".format(db["ref_genome"]))
        fout.write("--variant /data/{}/VCF/{}_RAW_INDEL.vcf ".format(args.name, args.name))
        fout.write("--disable_auto_index_creation_and_locking_when_reading_rods ")
        fout.write('--filterExpression "QD < 2.0" ')
        fout.write('--filterName "QD" ')
        fout.write('--filterExpression "DP < 8.0" ')
        fout.write('--filterName "DP" ')
        fout.write("--logging_level ERROR ")
        fout.write("-o /data/{}/VCF/{}_FILTERED_INDEL.vcf ".format(args.name, args.name))
        fout.write("-log /data/{}/LOG/9-2-2_{}_QC_filter_indel.log\n".format(args.name, args.name))

def make_QC_script(args, db):
    """ Create the script to run QC (FastQC and Trimmomatic)"""
    script_file = "/{}/{}/Scripts/{}_doQC.sh".format(db["out_dir"], args.name, args.name)
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("## Quality control of FASTQ files")
        fout.write("##-------------\n")
        fout.write(
            "(bash /{}/{}/Scripts/0a_{}_trimfastq.sh) 2>&1 | tee /{}/{}/LOG/0a_{}_trimfastqc.log\n".format(
                db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "(bash /{}/{}/Scripts/0b_{}_fastqc.sh) 2>&1 | tee /{}/{}/LOG/0b_{}_fastqc.log\n".format(
                db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name
            )
        )
    os.chmod(script_file, 0o777)


def make_masterscript(args, db):
    """ Create the master script used for running all the scripts """
    script_file = "/{}/{}/Scripts/{}_GATK.sh".format(db["out_dir"], args.name, args.name)
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("set -e\n")
        fout.write("##-------------\n")
        fout.write("##{}'s Variant Calling\n".format(args.name))
        fout.write("## Script version: {}\n".format(__version__))
        fout.write("##-------------\n")
        fout.write(
            "(bash /{}/{}/Scripts/1_{}_align.sh) 2>&1 | tee /{}/{}/LOG/1_{}_alignment.log\n".format(
                db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "(bash /{}/{}/Scripts/2_{}_sort.sh) 2>&1 | tee /{}/{}/LOG/2_{}_sort.log\n".format(
                db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "if [[ -e /{}/{}/BAM/{}_sorted.bam ]] ; then\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "    rm /{}/{}/SAM/{}_aligned.sam\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("fi\n")
        fout.write(
            "(bash /{}/{}/Scripts/3_{}_mark_duplicates.sh) 2>&1 | tee /{}/{}/LOG/3_{}_sort.log\n".format(
                db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "#(bash /{}/{}/Scripts/4_{}_build_index.sh) 2>&1 | tee /{}/{}/LOG/4_{}_building_index.log\n".format(
                db["out_dir"], args.name, args.name, db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "bash /{}/{}/Scripts/5_{}_realign_indels.sh\n".format(
                db["out_dir"], args.name, args.name
            )
        )

        fout.write(
            "bash /{}/{}/Scripts/6_{}_recalibrate_base.sh\n".format(
                db["out_dir"], args.name, args.name
            )
        )

        fout.write(
            "if [[ -e /{}/{}/BAM/{}_GATK.bam ]]; then\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "    rm -f /{}/{}/BAM/{}_{{sorted,realigned,deduplicated}}.{{bam,bai}}\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("fi\n")

        fout.write(
            "bash /{}/{}/Scripts/7_{}_call_haplotype.sh\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "bash /{}/{}/Scripts/8_{}_genotype_gvcf.sh\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "bash /{}/{}/Scripts/9_{}_SNV_quality_control.sh\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "bash /{}/{}/Scripts/10_{}_annotate.sh\n".format(
                db["out_dir"], args.name, args.name
            )
        )
    os.chmod(script_file, 0o777)


def make_annotate(args, db):
    """ Create the script for performing annotation on the VCF file """
    script_file = "/{}/{}/Scripts/10_{}_annotate.sh".format(
        db["out_dir"], args.name, args.name
    )
    #script_file = "{}_annotate.sh".format(args.name)
    with open(script_file, "w") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("# set -e\n")
        fout.write("\n")
        fout.write("##-------------\n")
        fout.write("##Step1: dbSNP\n")
        fout.write("##-------------\n")
        fout.write('echo "1/6 dbSNP Annotation Started"\n')
        fout.write("START_TIME=$SECONDS\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift annotate ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("/data/{} ".format(db["DBSNP"]))
        fout.write("/data/{}/VCF/{}_FILTERED_SNV.vcf ".format(args.name, args.name))
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_SNV.temp1.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift annotate ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("/data/{} ".format(db["DBSNP"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_INDEL.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_INDEL.temp1.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write("ELAPSED_TIME=$(($SECONDS - $START_TIME))\n")
        fout.write('echo "1/6 dbSNP Annotation Completed, Time: $ELAPSED_TIME sec"\n')
        fout.write("##-------------\n")
        fout.write("##Step2: dbNSFP\n")
        fout.write("##-------------\n")
        fout.write('echo "2/6 dbNSFP Annotation Started"\n')
        fout.write("START_TIME=$SECONDS\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift dbnsfp ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("-db /data/{} ".format(db["DBNSFP"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_SNV.temp1.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_SNV.temp2.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift dbnsfp ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("-db /data/{} ".format(db["DBNSFP"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_INDEL.temp1.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_INDEL.temp2.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write("ELAPSED_TIME=$(($SECONDS - $START_TIME))\n")
        fout.write('echo "2/6 dbNSFP Annotation Completed, Time: $ELAPSED_TIME sec"\n')
        fout.write(" ##-------------\n")
        fout.write("##Step3: gwasCat\n")
        fout.write("##-------------\n")
        fout.write('echo "3/6 gwasCat Annotation Started"\n')
        fout.write("START_TIME=$SECONDS\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift gwasCat ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("-db /data/{} ".format(db["GWASCATALOG"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_SNV.temp2.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_SNV.temp3.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift gwasCat ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("-db /data/{} ".format(db["GWASCATALOG"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_INDEL.temp2.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_INDEL.temp3.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write("ELAPSED_TIME=$(($SECONDS - $START_TIME))\n")
        fout.write('echo "3/6 gwasCat Annotation Completed, Time: $ELAPSED_TIME sec"\n')
        fout.write("##-------------\n")
        fout.write("##Step4: PhastCons\n")
        fout.write("##-------------\n")
        fout.write('echo "4/6 PhastCons Annotation Started"\n')
        fout.write("START_TIME=$SECONDS\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift phastCons ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("/data/{} ".format(db["PHASTCONS"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_SNV.temp3.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_SNV.temp4.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift phastCons ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("/data/{} ".format(db["PHASTCONS"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_INDEL.temp3.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_INDEL.temp4.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write("ELAPSED_TIME=$(($SECONDS - $START_TIME))\n")
        fout.write(
            'echo "4/6 PhastCons Annotation Completed, Time: $ELAPSED_TIME sec"\n'
        )
        fout.write("##-------------\n")
        fout.write("##Step5: ClinVar\n")
        fout.write("##-------------\n")
        fout.write('echo "5/6 ClinVar Annotation Started"\n')
        fout.write("START_TIME=$SECONDS\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift annotate ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("/data/{} ".format(db["CLINVAR"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_SNV.temp4.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_SNV.temp5.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpsift annotate ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("/data/{} ".format(db["CLINVAR"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_INDEL.temp4.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_INDEL.temp5.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write("ELAPSED_TIME=$(($SECONDS - $START_TIME))\n")
        fout.write('echo "5/6 ClinVar Annotation Completed, Time: $ELAPSED_TIME sec"\n')
        fout.write("##-------------\n")
        fout.write("##Step6: SnpEff\n")
        fout.write("##-------------\n")
        fout.write('echo "6/6 SnpEff Annotation Started"\n')
        fout.write("START_TIME=$SECONDS\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpeff ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("{} -t ".format(db["snpeff_dbver"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_SNV.temp5.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_SNV.temp6.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write(
            "docker run --rm -v {}:/data {} snpeff ".format(db["out_dir"], db["SNPEFF"])
        )
        fout.write("{} -t ".format(db["snpeff_dbver"]))
        fout.write(
            "/data/{}/VCF/{}_FILTERED_INDEL.temp5.vcf ".format(
                args.name, args.name
            )
        )
        fout.write(
            "> /{}/{}/VCF/{}_FILTERED_INDEL.temp6.vcf".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("\n\n")
        fout.write("ELAPSED_TIME=$(($SECONDS - $START_TIME))\n")
        fout.write('echo "6/6 SnpEff Annotation Completed, Time: $ELAPSED_TIME sec"\n')
        fout.write("##-------------\n")
        fout.write("##Step8: Clean-up\n")
        fout.write("##-------------\n")
        fout.write(
            "mv /{}/{}/VCF/{}_FILTERED_SNV.temp6.vcf ".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "/{}/{}/VCF/{}_SNV.vcf\n".format(db["out_dir"], args.name, args.name)
        )
        fout.write(
            "mv /{}/{}/VCF/{}_FILTERED_INDEL.temp6.vcf ".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "/{}/{}/VCF/{}_INDEL.vcf\n".format(db["out_dir"], args.name, args.name)
        )
        fout.write("\n\n")
        fout.write(
            "if [[ -e /{}/{}/VCF/{}_INDEL.vcf ]] ; then\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "    rm /{}/{}/VCF/{}_FILTERED_INDEL.temp?.vcf\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("fi\n")
        fout.write("\n\n")
        fout.write(
            "if [[ -e /{}/{}/VCF/{}_SNV.vcf ]] ; then\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write(
            "    rm /{}/{}/VCF/{}_FILTERED_SNV.temp?.vcf\n".format(
                db["out_dir"], args.name, args.name
            )
        )
        fout.write("fi\n")
        fout.write("\n\n")
        if args.replacement:
            fout.write(
                "rm /{}/{}/VCF/{}_FILTERED_INDEL.vcf\n".format(
                    db["out_dir"], args.name, args.name
                )
            )
            fout.write(
                "rm /{}/{}/VCF/{}_FILTERED_SNV.vcf\n".format(
                    db["out_dir"], args.name, args.name
                )
            )
    os.chmod(script_file, 0o777)


def variant_call_pipeline(args):
    print("Reading config file")
    cfg_db = read_config(args)
    print("Creating folder structure")
    do_setup(args, cfg_db)
    print("Generating scripts")
    make_preQC(args, cfg_db)
    make_postQC(args, cfg_db)
    make_trimfastq(args, cfg_db)
    make_align(args, cfg_db)
    make_sort(args, cfg_db)
    make_deduplicate(args, cfg_db)
    make_index(args, cfg_db)
    make_realign(args, cfg_db)
    make_BQSR(args, cfg_db)
    make_call_haplotype(args, cfg_db)
    make_genotype(args, cfg_db)
    make_SNV_QC(args, cfg_db)
    make_QC_script(args, cfg_db)
    make_masterscript(args, cfg_db)


def annotate_pipeline(args):
    ant_db = read_config(args)
    make_annotate(args, ant_db)


def main(args):
    variant_call_pipeline(args)
    annotate_pipeline(args)


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument(
        "name", action="store", help="Sample name (everything before '.1.fastq.gz')"
    )

    # Dependant on the user providing a region bed-file. Optional argument flag which defaults to True
    #parser.add_argument(
    #    "-e", "--exome", action="store_false", default=True, help="If not using Exome data"
    #)
    parser.add_argument(
        "-r",
        "--replacement",
        action="store_true",
        default=False,
        help="Keep intermediary annotation files",
    )

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-c", "--config", action="store", default="exome.cfg")
    parser.add_argument("-s", "--single", action="store_true", default=False, help="If single end reads")

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

#!/usr/bin/env python3
"""
Parses annotated VCF files from WGS/WES pipeline
"""

__author__ = "Harald Grove"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import sys


def read_genes(args, db):
    """ Read a list of genes to include in report, gene list can also contain
        information about transcript_ID, chrom, start, end, offset for conversion to hg19 """
    db["genes"] = {}
    with open(args.genes, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            line_l = line.strip().split()
            id_, chrom, start, end, conv = "NA", "NA", "NA", "NA", "NA"
            if len(line_l) < 1:
                continue
            name = line_l[0]
            if len(line_l) > 1:
                id_ = line_l[1]
            if len(line_l) > 2:
                chrom = line_l[2].strip("chr")
            if len(line_l) > 4:
                start, end = int(line_l[3]), int(line_l[4])
            if len(line_l) > 5:
                conv = int(line_l[5])
            db["genes"][name] = [id_, chrom, start, end, conv]


def read_samples(args, db):
    """ Read a list of samples to include in report """
    db["samples"] = []
    for sample_file in args.samples:
        with open(sample_file, "r") as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                newsample = line.strip()
                if len(newsample) == 0:
                    continue
                db["samples"].append(newsample)


def read_variants(args, db):
    """ Quick scan of which variants are present """
    db["scan"] = {}
    for sid in db["samples"]:
        for mode in ["SNV", "INDEL"]:
            vcf_file = "{}/VCF/{}_{}.vcf".format(sid, sid, mode)
            try:
                open(vcf_file, "r")
            except:
                continue
            with open(vcf_file, "r") as fin:
                for line in fin:
                    if line.startswith("#"):
                        continue
                    chrom, pos, id_, ref, alt, qual, filter_, info, format_, *samples = line.strip().split(
                        "\t"
                    )
                    chrom = chrom.strip("chr")
                    if filter_ != "PASS":
                        continue
                    db["scan"][chrom, pos] = [id_, ref, alt]


def read_vcfheader(args, db):
    """ Read the information fields that are present in the vcf file and should be included in the report """
    db["header_d"] = {}
    db["header_l"] = []
    db["ANN_header_l"] = []
    vcf_header_file = "exomeseq/vcf_header.txt"
    with open(vcf_header_file, "r") as fin:
        for line in fin:
            try:
                head, temp = line.split("=<")
            except:
                continue
            if head == "##INFO":
                try:
                    ID, Number, Type, Description = temp.strip(">").split(",", 3)
                except ValueError:
                    print(temp)
                    sys.exit()
                ID1, ID2 = ID.split("=")
                Number1, Number2 = Number.split("=")
                Type1, Type2 = Type.split("=")
                try:
                    Description1, Description2 = Description.split("=", 1)
                except ValueError:
                    print(Description)
                    sys.exit()
                if ID2 != "ANN":
                    db["header_l"].append(ID2)
                    db["header_d"][ID2] = {
                        "Number": Number2,
                        "Type": Type2,
                        "Description": Description2,
                    }
                else:
                    ann_header = Description2.strip('"').split("'")[1]
                    ann_header_l = ann_header.split("|")
                    for ahl in ann_header_l:
                        newkey = "ANN_{}".format(ahl.strip())
                        # header_l.append(newkey)
                        # header_d[newkey] = {'Number':'.','Type':'.','Description':'.'}
                        db["ANN_header_l"].append(newkey)


def read_dbSNP(args, db):
    """ Read information from dbSNP and 1KG extract """
    db["dbsnp"] = {}
    dbsnpfiles = ["/" + db["freq_main"]]
    for dbsnpfile in dbsnpfiles:
        with open(dbsnpfile, "r") as fin:
            for line in fin:
                allele = {}
                line_l = line.strip().split()
                chrom, pos, rs, chrom19, pos19, allelelist = line_l
                chrom = chrom.strip("chr")
                if (chrom, pos) not in db["scan"]:
                    continue
                if allelelist != "NA":
                    for al in allelelist.split(","):
                        # al = population:allele:frequency
                        p, a, f = al.split(":")
                        if a not in allele:
                            allele[a] = {}
                        allele[a][p] = float(f)
                db["dbsnp"][chrom, pos] = [rs, allele, chrom19, pos19]


def read_ExAC(args, db):
    """ Read information from ExAC """
    db["exac"] = {}
    dbsnpfiles = ["/" + db["exac_freqfile"]]
    for dbsnpfile in dbsnpfiles:
        with open(dbsnpfile, "r") as fin:
            for line in fin:
                allele = {}
                line_l = line.strip().split()
                chrom, pos, rs, chrom19, pos19, allelelist = line_l
                chrom = chrom.strip("chr")
                if (chrom, pos) not in db["scan"]:
                    continue
                if allelelist != "NA":
                    for al in allelelist.split(","):
                        # al = population:allele:frequency
                        p, a, f = al.split(":")
                        if a not in allele:
                            allele[a] = {}
                        allele[a][p] = float(f)
                db["exac"][chrom, pos] = [rs, allele, chrom19, pos19]


def read_ESP6500(args, db):
    """ Read information from ESP6500 extract """
    db["esp6500"] = {}
    dbsnpfiles = ["/" + db["esp6500_freqfile"]]
    for dbsnpfile in dbsnpfiles:
        with open(dbsnpfile, "r") as fin:
            for line in fin:
                allele = {}
                line_l = line.strip().split()
                chrom, pos, rs, chrom19, pos19, allelelist = line_l
                chrom = chrom.strip("chr")
                if (chrom, pos) not in db["scan"]:
                    continue
                if allelelist != "NA":
                    for al in allelelist.split(","):
                        # al = population:allele:frequency
                        try:
                            p, a, f = al.split(":")
                        except ValueError:
                            print("ERROR:", al)
                            continue
                        if a not in allele:
                            allele[a] = {}
                        allele[a][p] = float(f)
                db["esp6500"][chrom, pos] = [rs, allele, chrom19, pos19]


def read_HGVD(args, db):
    """ Read information from HGVD extract """
    db["hgvd"] = {}
    dbsnpfiles = ["/" + db["hgvd_freqfile"]]
    for dbsnpfile in dbsnpfiles:
        with open(dbsnpfile, "r") as fin:
            for line in fin:
                allele = {}
                line_l = line.strip().split()
                chrom, pos, rs, chrom19, pos19, allelelist = line_l
                if (chrom, pos) not in db["scan"]:
                    continue
                if allelelist != "NA":
                    for al in allelelist.split(","):
                        # al = population:allele:frequency
                        p, a, f = al.split(":")
                        if a not in allele:
                            allele[a] = {}
                        allele[a][p] = float(f)
                db["hgvd"][chrom19, pos19] = [rs, allele, chrom19, pos19]


def read_GONL(args, db):
    """ Read information from GoNL extract """
    db["gonl"] = {}
    dbsnpfiles = ["/" + db["gonl_freqfile"]]
    for dbsnpfile in dbsnpfiles:
        with open(dbsnpfile, "r") as fin:
            for line in fin:
                allele = {}
                line_l = line.strip().split()
                chrom, pos, rs, chrom19, pos19, allelelist = line_l
                if (chrom, pos) not in db["scan"]:
                    continue
                if allelelist != "NA":
                    for al in allelelist.split(","):
                        # al = population:allele:frequency
                        p, a, f = al.split(":")
                        if a not in allele:
                            allele[a] = {}
                        allele[a][p] = float(f)
                db["gonl"][chrom19, pos19] = [rs, allele, chrom19, pos19]


def read_clinvar(args, db):
    """ Read information from ClinVar extract """
    db["dbclinvar"] = {}
    clinvarfiles = ["/" + db["clinvar_freq"]]
    for name in clinvarfiles:
        with open(name, "r") as fin, open(args.logfile, "a") as flog:
            for line in fin:
                if line.startswith("#"):
                    continue
                chrom, pos, id_, ref, alt, qual, filter_, info = line.strip().split(
                    "\t"
                )
                if (chrom, pos) not in db["scan"]:
                    continue
                info_l = info.split(";")
                if (chrom, pos) in db["dbclinvar"]:
                    flog.write(
                        "WARNING: Multiple ClinVar records for {}:{}\n".format(
                            chrom, pos
                        )
                    )
                    continue
                db["dbclinvar"][chrom, pos] = {}
                for cv in info_l:
                    try:
                        key, value = cv.split("=")
                    except ValueError:
                        key, value = cv, True
                    db["dbclinvar"][chrom, pos][key] = value


def read_mutationtaster(args, db):
    """ Read information from MutationTaster database """
    db["dbmutationtaster"] = {}
    mutfiles = ["/" + db["mutationtaster"]]
    for name in mutfiles:
        with open(name, "r") as fin, open(args.logfile, "a") as flog:
            for line in fin:
                if line.startswith("#"):
                    continue
                chrom, pos, rs, chrom19, pos19, gene, pheno, allele = line.strip().split(
                    "\t"
                )
                if (chrom, pos) not in db["scan"]:
                    continue
                db["dbmutationtaster"][chrom19, pos19, gene] = [pheno, allele]


def read_config(args, db):
    """ Locate file containing regions of interest """
    with open(args.config, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            try:
                key, value = line.strip().split("=")
            except ValueError:
                continue
            db[key] = value


def read_regions(ars, db):
    region_file = "/{}/{}".format(db["out_dir"], db["exon_bed"])
    db["regions"] = {}
    with open(region_file, "r") as fin:
        for reg_line in fin:
            l = reg_line.strip().split()
            if len(l) == 4:
                chrom, start, end, name = l
            elif len(l) == 3:
                chrom, start, end = l
                name = "NA"
            else:
                continue
            chrom = chrom.strip("chr")
            if chrom not in db["regions"]:
                db["regions"][chrom] = []
            db["regions"][chrom].append([start, end, name])


def do_setup(args, db):
    """ Create files and check that all is in order """
    with open(args.outfile, "w") as fout:
        fout.write(
            "SAMPLE\tCHR\tPOS\tID\tREF\tALT\t{}\t{}\n".format(
                "\t".join(db["header_l"]), "\t".join(db["ANN_header_l"])
            )
        )
    with open(args.logfile, "w") as fout:
        fout.write("# Logfile\n")
    if args.reportfile is None:
        return
    with open(args.reportfile, "w") as fout:
        fout.write(
            (
                "SAMPLE\tCHR\tPOS\tID\tCHR19\tPOS19\tREF\tALT\tREF_AD\tALT_AD\tGENE\t"
                "KHV_ALT-AF\tCDX_ALT-AF\tAFR_ALT-AF\tAMR_ALT-AF\tEAS_ALT-AF\tEUR_ALT-AF\tSAS_ALT-AF\t"
                "ExAC_raw_AF\tExAC_Female_AF\tExAC_Male_AF\tExAC_AFR_AF\tExAC_AMR_AF\tExAC_ASJ_AF\t"
                "ExAC_EAS_AF\tExAC_FIN_AF\tExAC_NFE_AF\tExAC_OTH_AF\tExAC_SAS_AF\tExAC_POPMAX_AF\t"
                "ExAC_AFR_Male_AF\tExAC_AMR_Male_AF\tExAC_ASJ_Male_AF\tExAC_EAS_Male_AF\t"
                "ExAC_FIN_Male_AF\tExAC_NFE_Male_AF\tExAC_OTH_Male_AF\tExAC_SAS_Male_AF\tExAC_AFR_Female_AF\t"
                "ExAC_AMR_Female_AF\tExAC_ASJ_Femlae_AF\tExAC_EAS_Female_AF\t"
                "ExAC_FIN_Female_AF\tExAC_NFE_Female_AF\tExAC_OTH_Female_AF\tExAC_SAS_Female_AF\t"
                "EA-AF\tAA-AF\tHGVD\tGONL-AF\t"
                "CLNSIG\tdbNSFP_MutationTaster_pred\tdbNSFP_SIFT_pred\t"
                "dbNSFP_phastCons100way_vertebrate\tdbNSFP_FATHMM_pred\tdbNSFP_MetaSVM_pred\n"
            )
        )
    db["poplist"] = ["KHV", "CDX", "AFR", "AMR", "EAS", "EUR", "SAS"]
    db["exacpoplist"] = [
        "raw",
        "Female",
        "Male",
        "AFR",
        "AMR",
        "ASJ",
        "EAS",
        "FIN",
        "NFE",
        "OTH",
        "SAS",
        "POPMAX",
        "AFR_Male",
        "AMR_Male",
        "ASJ_Male",
        "EAS_Male",
        "FIN_Male",
        "NFE_Male",
        "OTH_Male",
        "SAS_Male",
        "AFR_Female",
        "AMR_Female",
        "ASJ_Female",
        "EAS_Female",
        "FIN_Female",
        "NFE_Female",
        "OTH_Female",
        "SAS_Female",
    ]


def parse_vcf(args, db, sample, mode):
    """ Read the vcf files for one sample """
    vcf_file = "{}/VCF/{}_{}.vcf".format(sample, sample, mode)
    try:
        open(vcf_file, "r")
    except:
        sys.stderr.write("Skipping {} {}, missing files\n".format(sample, mode))
        return
    if args.reportfile is not None:
        ffreq = open(args.reportfile, "a")
    with open(vcf_file, "r") as fin, open(args.outfile, "a") as fout,\
            open(args.logfile, "a") as flog:
        for line in fin:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                title = line.strip("#").strip().split("\t")
                continue
            chrom, pos, id_, ref, alt, qual, filter_, info, format_, *samples = line.strip().split(
                "\t"
            )
            if filter_ != "PASS":
                continue
            chrom = chrom.strip("chr")
            # Update rs-info and collect extra allelefrequency information
            if (chrom, pos) in db["dbsnp"]:
                rs, allelefreq, chrom19, pos19 = db["dbsnp"][(chrom, pos)]
                if id_ == ".":
                    id_ = rs
                if rs not in id_.split(";"):
                    flog.write(
                        "Non-matching rs number: {}\t{}\t{}\t{}\n".format(
                            chrom, pos, id_, rs
                        )
                    )
            else:
                rs, allelefreq, chrom19, pos19 = ".", {}, ".", "."
            if (chrom, pos) in db["esp6500"]:
                esp6500_frq = db["esp6500"][(chrom, pos)][1]
            else:
                esp6500_frq = {}
            if (chrom, pos) in db["exac"]:
                exac_frq = db["exac"][(chrom, pos)][1]
            else:
                exac_frq = {}
            if (chrom19, pos19) in db["hgvd"]:
                hgvd_frq = db["hgvd"][("" + chrom19, pos19)][1]
            else:
                hgvd_frq = {}
            if (chrom19, pos19) in db["gonl"]:
                gonl_frq = db["gonl"][(chrom19, pos19)][1]
            else:
                gonl_frq = {}
            # Gather information about target region
            target_region = "NA"
            if chrom in db["regions"]:
                for reg_start, reg_end, reg_name in db["regions"][chrom]:
                    if int(pos) >= int(reg_start) and int(pos) <= int(reg_end):
                        target_region = reg_name
                        break
            # Collect information about the genotype
            refAD, altAD = "NA", "NA"
            if len(samples) == 1:
                # Locate the allele depth (AD in the format-field)
                AD_ind = format_.split(":").index("AD")
                try:
                    refAD, altAD = samples[0].split(":")[AD_ind].split(",")
                except ValueError:
                    flog.write(
                        "WARNING:\tUnexpected AD field: [{},{}] [{}]\n".format(
                            ref, alt, samples[0]
                        )
                    )
            # Split up info field and store all the elements
            info_d = {}
            info_l = info.split(";")
            for lel in info_l:
                if "=" in lel:
                    key, value = lel.split("=")
                    info_d[key] = value
                else:
                    info_d[lel] = True
            # Select the right ClinVar record (Clinvar records are separated by ',')
            if "CLNHGVS" in info_d:
                cvlist = info_d["CLNHGVS"].split(",")
                CLNEntry = None
                for i, e in enumerate(cvlist):
                    if ">" in e:
                        try:
                            e1, e2 = e.split(">")
                        except ValueError:
                            flog.write(
                                'INFO:\tToo many elements in CLNHGVS:\t"{}"\t{}\t{}\t{}\n'.format(
                                    e, sample, chrom, pos
                                )
                            )
                            continue
                        if len(e2) != 1:
                            flog.write(
                                'INFO:\tUnexpected format of CLNHGVS:\t"{}"\t{}\t{}\t{}\n'.format(
                                    e, sample, chrom, pos
                                )
                            )
                            continue
                        if e2 != alt:
                            continue
                        CLNEntry = i
                        break
                    elif "del" in e:
                        try:
                            e1, e2 = e.split("del")
                        except ValueError:
                            flog.write(
                                'INFO:\tToo many elements in CLNHGVS:\t"{}"\t{}\t{}\t{}\n'.format(
                                    e, sample, chrom, pos
                                )
                            )
                            continue
                        if ref != alt + e2:
                            continue
                        CLNEntry = i
                        break
                    elif "dup" in e:
                        try:
                            e1, e2 = e.split("dup")
                        except ValueError:
                            flog.write(
                                'INFO:\tToo many elements in CLNHGVS:\t"{}"\t{}\t{}\t{}\n'.format(
                                    e, sample, chrom, pos
                                )
                            )
                            continue
                        if ref + e2 != alt:
                            continue
                        CLNEntry = i
                        break
                else:
                    flog.write(
                        'INFO:\tNo match to ALT in CLNHGVS:\t{}\t{}\t"{}"\t{}\t{}\t{}\n'.format(
                            ref, alt, cvlist, sample, chrom, pos
                        )
                    )
                    # sys.exit(1)
                # Remove CLINVAR entries that doesn't match the actual mutation
                if CLNEntry is not None:
                    for key in info_d:
                        if not key.startswith("CLN"):
                            continue
                        value = info_d[key]
                        info_d[key] = value.split(",")[CLNEntry]
            else:
                # flog.write('INFO:\tNo entry named CLNHGVS:\t{}\t{}\t{}\n'.format(sample,chrom,pos))
                pass
            # Split ANN field, this part is the slowest section.
            ann_l = info_d["ANN"].split(",")  # Each transcript
            for transcript in ann_l:
                ANN_info_d = {}
                t_el = transcript.split("|")
                for key, value in zip(db["ANN_header_l"], t_el):
                    ANN_info_d[key] = value
                # Decide if variant should be kept or skipped
                if len(db["genes"]) > 0:
                    # 1. Variant is not in a relevant gene
                    if ANN_info_d["ANN_Gene_Name"] not in db["genes"]:
                        continue
                    # 2. Variant is not in the right transcript
                    if (
                        db["genes"][ANN_info_d["ANN_Gene_Name"]][0] != "NA"
                        and db["genes"][ANN_info_d["ANN_Gene_Name"]][0]
                        not in ANN_info_d["ANN_Feature_ID"]
                    ):
                        continue
                # A colon in transcript name is ?
                if ":" in ANN_info_d["ANN_Feature_ID"]:
                    continue
                # If information is not known about hg19, estimate it
                if pos19 == "." and chrom19 == "." and len(db["genes"]) > 0:
                    pos19 = int(pos) + db["genes"][ANN_info_d["ANN_Gene_Name"]][4]
                    chrom19 = db["genes"][ANN_info_d["ANN_Gene_Name"]][1]
                # Write full report file
                fout.write(
                    "{}\t{}\t{}\t{}\t{}\t{}".format(sample, chrom, pos, id_, ref, alt)
                )
                for col in db["header_l"]:
                    val = info_d.get(col, "NA")
                    fout.write("\t{}".format(val))
                for col in db["ANN_header_l"]:
                    val = ANN_info_d.get(col, "NA")
                    fout.write("\t{}".format(val))
                fout.write("\n")
            if args.reportfile is None:
                continue
            # Write variant file for preliminary filtering
            # Fields: (Read-depth, frequency, clinvar-sig,
            ffreq.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    sample,
                    chrom,
                    pos,
                    id_,
                    chrom19.strip("chr"),
                    pos19,
                    ref,
                    alt,
                    refAD,
                    altAD,
                    target_region,
                )
            )
            for pop in db["poplist"]:
                try:
                    ffreq.write("\t{}".format(allelefreq[alt][pop]))
                except KeyError:
                    ffreq.write("\tNA")
            for pop in ["EA", "AA"]:
                try:
                    ffreq.write("\t{}".format(esp6500_frq[alt][pop]))
                except KeyError:
                    ffreq.write("\tNA")
            for pop in db["exacpoplist"]:
                try:
                    ffreq.write("\t{}".format(exac_frq[alt][pop]))
                except KeyError:
                    ffreq.write("\tNA")
            try:
                ffreq.write("\t{}".format(hgvd_frq[alt]["HGVD"]))
            except KeyError:
                ffreq.write("\tNA")
            try:
                ffreq.write("\t{}".format(gonl_frq[alt]["GONL"]))
            except KeyError:
                ffreq.write("\tNA")
            try:
                ffreq.write("\t{}".format(db["dbclinvar"][chrom, pos]["CLNSIG"]))
            except KeyError:
                ffreq.write("\tNA")
            try:
                pheno, allele = db["dbmutationtaster"][chrom, pos19, target_region]
                if ref == allele[0] and alt == allele[2]:
                    ffreq.write("\t{}".format(pheno))
                else:
                    ffreq.write("\tNA")
            except KeyError:
                ffreq.write("\tNA")
            try:
                ffreq.write("\t{}".format(info_d["dbNSFP_SIFT_pred"]))
            except KeyError:
                ffreq.write("\tNA")
            try:
                ffreq.write("\t{}".format(info_d["dbNSFP_phastCons100way_vertebrate"]))
            except KeyError:
                ffreq.write("\tNA")
            try:
                ffreq.write("\t{}".format(info_d["dbNSFP_FATHMM_pred"]))
            except KeyError:
                ffreq.write("\tNA")
            try:
                ffreq.write("\t{}".format(info_d["dbNSFP_MetaSVM_pred"]))
            except KeyError:
                ffreq.write("\tNA")

            ffreq.write("\n")


def parse_vcfs(args, db):
    """ Read all vcf files for all the samples """
    for sid in db["samples"]:
        for mode in ["SNV", "INDEL"]:
            parse_vcf(args, db, sid, mode)


def main(args):
    """ """
    import time

    db = {}
    try:
        t = time.time()
        read_genes(args, db)
        print("{}: {:.3f}s".format(args.genes, time.time() - t))
    except FileNotFoundError:
        db["genes"] = {}
    try:
        t = time.time()
        read_samples(args, db)
        print("Read {} samples: {:.3f}s".format(len(db["samples"]), time.time() - t))
    except FileNotFoundError:
        db["samples"] = []
    t = time.time()
    read_variants(args, db)
    print("Scanned VCF file: {:.3f}s".format(time.time() - t))
    try:
        t = time.time()
        read_config(args, db)
        print("{}: {:.3f}s".format(args.config, time.time() - t))
    except FileNotFoundError:
        sys.stderr.write("Missing configuration file: {}\n".format(args.config))
        return
    try:
        t = time.time()
        read_regions(args, db)
        print("Read regions file: {:.3f}s".format(time.time() - t))
    except FileNotFoundError:
        db["regions"] = {}
    try:
        t = time.time()
        read_vcfheader(args, db)
        print("VCFheaders: {:.3f}s".format(time.time() - t))
    except FileNotFoundError:
        db["header_d"] = {}
        db["header_l"] = []
        db["ANN_header_l"] = []
    try:
        t = time.time()
        do_setup(args, db)
        print("Setup: {:.3f}s".format(time.time() - t))
    except:
        sys.stderr.write("Error performing setup\n")
        return
    db["dbsnp"] = {}
    db['dbmutationtaster'] = {}
    db['dbclinvar'] = {}
    db["esp6500"] = {}
    db["exac"] = {}
    db["hgvd"] = {}
    db["gonl"] = {}
    if args.reportfile is not None:
        try:
            t = time.time()
            print("dbSNP: ", end="")
            read_dbSNP(args, db)
            print("{:.3f}s".format(time.time() - t))
        except FileNotFoundError:
            print("dbSNP database not found.")

        try:
            t = time.time()
            read_mutationtaster(args, db)
            print("mutationtaster: {:.3f}s".format(time.time() - t))
        except FileNotFoundError:
            print("Mutationtaster database not found.")
        try:
            t = time.time()
            read_clinvar(args, db)
            print("ClinVar: {:.3f}s".format(time.time() - t))
        except FileNotFoundError:
            print("ClinVar database not found.")
        try:
            t = time.time()
            read_ESP6500(args, db)
            print("ESP6500: {:.3f}s".format(time.time() - t))
        except FileNotFoundError:
            print("ESP6500 database not found.")
        try:
            t = time.time()
            read_ExAC(args, db)
            print("ExAc: {:.3f}s".format(time.time() - t))
        except FileNotFoundError:
            print("ExAc database not found.")
        try:
            t = time.time()
            read_HGVD(args, db)
            print("HGVD: {:.3f}s".format(time.time() - t))
        except FileNotFoundError:
            print("HGVD database not found.")
        try:
            t = time.time()
            print("GONL: ", end="")
            read_GONL(args, db)
            print("{:.3f}s".format(time.time() - t))
        except FileNotFoundError:
            print("GoNL database not found.")
    # read_params(args)
    t = time.time()
    parse_vcfs(args, db)
    print("Applying updates: {:.3f}s".format(time.time() - t))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o",
        "--outfile",
        action="store",
        default="full_report.txt",
        help="Name of full output",
    )
    parser.add_argument(
        "-r",
        "--reportfile",
        action="store",
        help="Name of allelel frequency report",
    )
    parser.add_argument(
        "-t",
        "--transcripts",
        action="store",
        help="File with list of the transcripts to include in report",
    )
    parser.add_argument(
        "-c", "--config", action="store", help="Configuration file", default="exome.cfg"
    )
    parser.add_argument(
        "-g", "--genes", action="store", help="Genes of interest", default="genes.list"
    )
    parser.add_argument(
        "-s",
        "--samples",
        action="store",
        help="Samples to include in report",
        nargs='*',
        default=["samples.paired.list", "samples.single.list"],
    )
    # parser.add_argument("-i", "--info", action="store", help="List of INFO fields to keep", default="info_fields.txt")
    parser.add_argument(
        "-l", "--logfile", action="store", help="Logfile", default="logfile.log"
    )

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

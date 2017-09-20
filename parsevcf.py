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
    db['genes'] = {}
    with open(args.genes, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            line_l = line.strip().split()
            id_, chrom, start, end, conv = 'NA', 'NA', 'NA', 'NA', 'NA'
            if len(line_l) < 1:
                continue
            name = line_l[0]
            if len(line_l) > 1:
                id_ = line_l[1]
            if len(line_l) > 2:
                chrom = line_l[2]
            if len(line_l) > 4:
                start, end = int(line_l[3]), int(line_l[4])
            if len(line_l) > 5:
                conv = int(line_l[5])
            db['genes'][name] = [id_, chrom, start, end, conv]


def read_samples(args, db):
    """ Read a list of samples to include in report """
    db['samples'] = []
    with open(args.samples, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            db['samples'].append(line.strip())

def read_vcfheader(args, db):
    """ Read the information fields that are present in the vcf file and should be included in the report """
    db['header_d'] = {}
    db['header_l'] = []
    db['ANN_header_l'] = []
    vcf_header_file = 'vcf_header.txt' #TODO: think about how to handle this file
    with open(vcf_header_file,'r') as fin:
        for line in fin:
            try:
                head, temp = line.split('=<')
            except:
                continue
            if head == '##INFO':
                try:
                    ID,Number,Type,Description = temp.strip('>').split(',',3)
                except ValueError:
                    print(temp)
                    sys.exit()
                ID1,ID2 = ID.split('=')
                Number1,Number2 = Number.split('=')
                Type1,Type2 = Type.split('=')
                try:
                    Description1,Description2 = Description.split('=',1)
                except ValueError:
                    print(Description)
                    sys.exit()
                if ID2 != 'ANN':
                    db['header_l'].append(ID2)
                    db['header_d'][ID2] = {'Number':Number2, 'Type':Type2, 'Description':Description2}
                else:
                    ann_header = Description2.strip('"').split("'")[1]
                    ann_header_l = ann_header.split('|')
                    for ahl in ann_header_l:
                        newkey = 'ANN_{}'.format(ahl.strip())
                        #header_l.append(newkey)
                        #header_d[newkey] = {'Number':'.','Type':'.','Description':'.'}
                        db['ANN_header_l'].append(newkey)

def read_dbSNP(args, db):
    """ Read information from dbSNP and 1KG extract """
    db['dbsnp'] = {}
    dbsnpfiles = ['/'+db['freq_main'], '/'+db['freq_new']]
    for dbsnpfile in dbsnpfiles:
        with open(dbsnpfile, 'r') as fin:
            for line in fin:
                allele = {}
                line_l = line.strip().split()
                chrom, pos, rs, chrom19, pos19, allelelist = line_l
                if allelelist != 'NA':
                    for al in allelelist.split(','):
                        # al = population:allele:frequency
                        p, a, f = al.split(':')
                        if a not in allele:
                            allele[a] = {}
                        allele[a][p] = float(f)
                db['dbsnp'][chrom,pos] = [rs, allele, chrom19, pos19]

def read_ESP6500(args, db):
    """ Read information from ESP6500 extract """
    db['esp6500'] = {}
    dbsnpfiles = ['/'+db['esp6500_freqfile']]
    for dbsnpfile in dbsnpfiles:
        with open(dbsnpfile, 'r') as fin:
            for line in fin:
                allele = {}
                line_l = line.strip().split()
                chrom, pos, rs, chrom19, pos19, allelelist = line_l
                if allelelist != 'NA':
                    for al in allelelist.split(','):
                        # al = population:allele:frequency
                        try:
                            p, a, f = al.split(':')
                        except ValueError:
                            print('ERROR:',al)
                            continue
                        if a not in allele:
                            allele[a] = {}
                        allele[a][p] = float(f)
                db['esp6500'][chrom,pos] = [rs, allele, chrom19, pos19]

def read_HGVD(args, db):
    """ Read information from HGVD extract """
    db['hgvd'] = {}
    dbsnpfiles = ['/'+db['hgvd_freqfile']]
    for dbsnpfile in dbsnpfiles:
        with open(dbsnpfile, 'r') as fin:
            for line in fin:
                allele = {}
                line_l = line.strip().split()
                chrom, pos, rs, chrom19, pos19, allelelist = line_l
                if allelelist != 'NA':
                    for al in allelelist.split(','):
                        # al = population:allele:frequency
                        p, a, f = al.split(':')
                        if a not in allele:
                            allele[a] = {}
                        allele[a][p] = float(f)
                db['hgvd'][chrom19,pos19] = [rs, allele, chrom19, pos19]

def read_GONL(args, db):
    """ Read information from GoNL extract """
    db['gonl'] = {}
    dbsnpfiles = ['/'+db['gonl_freqfile']]
    for dbsnpfile in dbsnpfiles:
        with open(dbsnpfile, 'r') as fin:
            for line in fin:
                allele = {}
                line_l = line.strip().split()
                chrom, pos, rs, chrom19, pos19, allelelist = line_l
                if allelelist != 'NA':
                    for al in allelelist.split(','):
                        # al = population:allele:frequency
                        p, a, f = al.split(':')
                        if a not in allele:
                            allele[a] = {}
                        allele[a][p] = float(f)
                db['gonl'][chrom19,pos19] = [rs, allele, chrom19, pos19]

def read_clinvar(args, db):
    """ Read information from ClinVar extract """
    db['dbclinvar'] = {}
    clinvarfiles = ['/'+db['clinvar_freq1'], '/'+db['clinvar_freq2'], '/'+db['clinvar_freq3']]
    for name in clinvarfiles:
        with open(name, 'r') as fin, open(args.logfile, 'a') as flog:
            for line in fin:
                if line.startswith('#'):
                    continue
                chrom,pos,id_,ref,alt,qual,filter_,info = line.strip().split('\t')
                info_l = info.split(';')
                if (chrom, pos) in db['dbclinvar']:
                    flog.write('WARNING: Multiple ClinVar records for {}:{}\n'.format(chrom,pos))
                    continue
                db['dbclinvar'][chrom,pos] = {}
                for cv in info_l:
                    try:
                        key, value = cv.split('=')
                    except ValueError:
                        key, value = cv, True
                    db['dbclinvar'][chrom,pos][key] = value

def read_mutationtaster(args, db):
    """ Read information from MutationTaster database """
    db['dbmutationtaster'] = {}
    mutfiles = ['/'+db['mutationtaster']]
    for name in mutfiles:
        with open(name, 'r') as fin, open(args.logfile, 'a') as flog:
            for line in fin:
                if line.startswith('#'):
                    continue
                chrom,pos,gene,pheno,allele = line.strip().split('\t')
                db['dbmutationtaster'][chrom,pos,gene] = [pheno, allele]

def read_regions(args, db):
    """ Locate file containing regions of interest """
    with open(args.config, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            try:
                key, value = line.strip().split('=')
            except ValueError:
                continue
            db[key] = value

def do_setup(args, db):
    """ Create files and check that all is in order """
    with open(args.outfile, 'w') as fout:
        fout.write('SAMPLE\tCHR\tPOS\tID\tREF\tALT\t{}\t{}\n'.format('\t'.join(db['header_l']),'\t'.join(db['ANN_header_l'])))
    with open(args.logfile, 'w') as fout:
        fout.write('# Logfile\n')
    with open(args.reportfile, 'w') as fout:
        fout.write(('SAMPLE\tCHR\tPOS\tID\tCHR19\tPOS19\tREF\tALT\tREF_AD\tALT_AD\tGENE\t'
                   'KHV_ALT-AF\tCDX_ALT-AF\tAFR_ALT-AF\tAMR_ALT-AF\tEAS_ALT-AF\tEUR_ALT-AF\tSAS_ALT-AF\t'
                   'EA-AF\tAA-AF\tHGVD\tGONL-AF\t'
                   'CLNSIG\tdbNSFP_MutationTaster_pred\tdbNSFP_SIFT_pred\t'
                   'dbNSFP_phastCons100way_vertebrate\tdbNSFP_FATHMM_pred\tdbNSFP_MetaSVM_pred\n'))
    db['poplist'] = ['KHV', 'CDX', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']

def parse_vcf(args, db, sample, mode):
    """ Read the vcf files for one sample """
    vcf_file = '{}/VCF/{}_{}.vcf'.format(sample, sample, mode)
    try:
        open(vcf_file,'r')
    except:
        sys.stderr.write('Skipping {} {}, missing files\n'.format(sample,mode))
        return
    with open(vcf_file, 'r') as fin, open(args.outfile, 'a') as fout, open(args.logfile, 'a') as flog, open(args.reportfile, 'a') as ffreq:
        for line in fin:
            if line.startswith('##'): continue
            if line.startswith('#'):
                title = line.strip('#').strip().split('\t')
                continue
            chrom,pos,id_,ref,alt,qual,filter_,info,format_,*samples = line.strip().split('\t')
            if filter_ != 'PASS':
                continue
            # Update rs-info and collect extra allelefrequency information
            if (chrom[3:], pos) in db['dbsnp']:
                rs, allelefreq, chrom19, pos19 = db['dbsnp'][(chrom[3:],pos)]
                if id_ == '.':
                    id_ = rs
                if rs not in id_.split(';'):
                    flog.write('Non-matching rs number: {}\t{}\t{}\t{}\n'.format(chrom,pos,id_,rs))
            else:
                rs, allelefreq, chrom19, pos19 = '.', {}, '.', '.'
            if (chrom[3:], pos) in db['esp6500']:
                esp6500_frq = db['esp6500'][(chrom[3:], pos)][1]
            else:
                esp6500_frq = {}
            if (chrom19, pos19) in db['hgvd']:
                hgvd_frq = db['hgvd'][(''+chrom19, pos19)][1]
            else:
                hgvd_frq = {}
            if (chrom19, pos19) in db['gonl']:
                gonl_frq = db['gonl'][(chrom19, pos19)][1]
            else:
                gonl_frq = {}
            # Gather information about target region
            target_region = 'NA'
            region_file = '/'+db['exon_bed']
            with open(region_file, 'r') as reg_fin:
                for reg_line in reg_fin:
                    try:
                        reg_chrom, reg_start, reg_end, reg_name = reg_line.strip().split()
                    except ValueError:
                        continue
                    if chrom == reg_chrom and int(pos) > int(reg_start) and int(pos) < int(reg_end):
                        target_region = reg_name
                        break
            # Collect information about the genotype
            refAD, altAD = 'NA','NA'
            if len(samples) == 1:
                # Locate the allele depth (AD in the format-field)
                AD_ind = format_.split(':').index('AD')
                try:
                    refAD, altAD = samples[0].split(':')[AD_ind].split(',')
                except ValueError:
                    flog.write('WARNING:\tUnexpected AD field: [{},{}] [{}]\n'.format(ref,alt,samples[0]))
            # Split up info field and store all the elements
            info_d = {}
            info_l = info.split(';')
            for lel in info_l:
                if '=' in lel:
                    key,value = lel.split('=')
                    info_d[key] = value
                else:
                    info_d[lel] = True
            # Select the right ClinVar record (Clinvar records are separated by ',')
            if 'CLNHGVS' in info_d:
                cvlist = info_d['CLNHGVS'].split(',')
                CLNEntry = None
                for i,e in enumerate(cvlist):
                    if '>' in e:
                        try:
                            e1,e2 = e.split('>')
                        except ValueError:
                            flog.write('INFO:\tToo many elements in CLNHGVS:\t"{}"\t{}\t{}\t{}\n'.format(e,sample,chrom,pos))
                            continue
                        if len(e2) != 1:
                            flog.write('INFO:\tUnexpected format of CLNHGVS:\t"{}"\t{}\t{}\t{}\n'.format(e,sample,chrom,pos))
                            continue
                        if e2 != alt:
                            continue
                        CLNEntry = i
                        break
                    elif 'del' in e:
                        try:
                            e1,e2 = e.split('del')
                        except ValueError:
                            flog.write('INFO:\tToo many elements in CLNHGVS:\t"{}"\t{}\t{}\t{}\n'.format(e,sample,chrom,pos))
                            continue
                        if ref != alt + e2:
                            continue
                        CLNEntry = i
                        break
                    elif 'dup' in e:
                        try:
                            e1,e2 = e.split('dup')
                        except ValueError:
                            flog.write('INFO:\tToo many elements in CLNHGVS:\t"{}"\t{}\t{}\t{}\n'.format(e,sample,chrom,pos))
                            continue
                        if ref + e2 != alt:
                            continue
                        CLNEntry = i
                        break
                else:
                    flog.write('INFO:\tNo match to ALT in CLNHGVS:\t{}\t{}\t"{}"\t{}\t{}\t{}\n'.format(ref,alt,cvlist,sample,chrom,pos))
                    #sys.exit(1)
                # Remove CLINVAR entries that doesn't match the actual mutation
                if CLNEntry is not None:
                    for key in info_d:
                        if not key.startswith('CLN'):
                            continue
                        value = info_d[key]
                        info_d[key] = value.split(',')[CLNEntry]
            else:
                #flog.write('INFO:\tNo entry named CLNHGVS:\t{}\t{}\t{}\n'.format(sample,chrom,pos))
                pass
            # Split ANN field
            ann_l = info_d['ANN'].split(',') # Each transcript
            for transcript in ann_l:
                ANN_info_d = {}
                t_el = transcript.split('|')
                for key,value in zip(db['ANN_header_l'],t_el):
                    ANN_info_d[key] = value
                # Decide if variant should be kept or skipped
                # 1. Variant is not in a relevant gene
                if len(db['genes']) > 0 and ANN_info_d['ANN_Gene_Name'] not in db['genes']:
                    continue
                # 2. Variant is not in the right transcript
                if db['genes'][ANN_info_d['ANN_Gene_Name']][0] != 'NA' and \
                   db['genes'][ANN_info_d['ANN_Gene_Name']][0] not in ANN_info_d['ANN_Feature_ID']:
                    continue
                if ':' in ANN_info_d['ANN_Feature_ID']:
                    continue
                # If information is not known about hg19, estimate it
                if pos19 == '.' and chrom19 == '.':
                    pos19 = int(pos) + db['genes'][ANN_info_d['ANN_Gene_Name']][4]
                    chrom19 = db['genes'][ANN_info_d['ANN_Gene_Name']][1]
                # Write full report file
                fout.write('{}\t{}\t{}\t{}\t{}\t{}'.format(sample,chrom,pos,id_,ref,alt))
                for col in db['header_l']:
                    val = info_d.get(col,'NA')
                    fout.write('\t{}'.format(val))
                for col in db['ANN_header_l']:
                    val = ANN_info_d.get(col,'NA')
                    fout.write('\t{}'.format(val))
                fout.write('\n')
            # Write variant file for preliminary filtering
            # Fields: (Read-depth, frequency, clinvar-sig,
            ffreq.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.\
                        format(sample,chrom,pos,id_,chrom19.strip('chr'),pos19,ref,alt,refAD,altAD, target_region))
            for pop in db['poplist']:
                try:
                    ffreq.write('\t{}'.format(allelefreq[alt][pop]))
                except KeyError:
                    ffreq.write('\tNA')
            for pop in ['EA', 'AA']:
                try:
                    ffreq.write('\t{}'.format(esp6500_frq[alt][pop]))
                except KeyError:
                    ffreq.write('\tNA')
            try:
                ffreq.write('\t{}'.format(hgvd_frq[alt]['HGVD']))
            except KeyError:
                ffreq.write('\tNA')
            try:
                ffreq.write('\t{}'.format(gonl_frq[alt]['GONL']))
            except KeyError:
                ffreq.write('\tNA')
            try:
                ffreq.write('\t{}'.format(db['dbclinvar'][chrom[3:],pos]['CLNSIG']))
            except KeyError:
                ffreq.write('\tNA')
            try:
                pheno, allele = db['dbmutationtaster'][chrom[3:],pos19,target_region]
                if ref == allele[0] and alt == allele[2]:
                    ffreq.write('\t{}'.format(pheno))
                else:
                    ffreq.write('\tNA')
            except KeyError:
                ffreq.write('\tNA')
            try:
                ffreq.write('\t{}'.format(info_d['dbNSFP_SIFT_pred']))
            except KeyError:
                ffreq.write('\tNA')
            try:
                ffreq.write('\t{}'.format(info_d['dbNSFP_phastCons100way_vertebrate']))
            except KeyError:
                ffreq.write('\tNA')
            try:
                ffreq.write('\t{}'.format(info_d['dbNSFP_FATHMM_pred']))
            except KeyError:
                ffreq.write('\tNA')
            try:
                ffreq.write('\t{}'.format(info_d['dbNSFP_MetaSVM_pred']))
            except KeyError:
                ffreq.write('\tNA')

            ffreq.write('\n')

def parse_vcfs(args, db):
    """ Read all vcf files for all the samples """
    for sid in db['samples']:
        for mode in ['SNV', 'INDEL']:
            parse_vcf(args, db, sid, mode)

def main(args):
    """ """
    import time
    db = {}
    t = time.time()
    read_genes(args, db)
    print('Genes: {:.3f}s'.format(time.time()-t))
    t = time.time()
    read_samples(args, db)
    print('Samples: {:.3f}s'.format(time.time()-t))
    t = time.time()
    read_regions(args, db)
    print('Regions: {:.3f}s'.format(time.time()-t))
    t = time.time()
    read_vcfheader(args, db)
    print('VCFheaders: {:.3f}s'.format(time.time()-t))
    t = time.time()
    do_setup(args, db)
    print('Setup: {:.3f}s'.format(time.time()-t))
    t = time.time()
    read_dbSNP(args, db)
    print('dbSNP: {:.3f}s'.format(time.time()-t))
    t = time.time()
    read_mutationtaster(args, db)
    print('mutationtaster: {:.3f}s'.format(time.time()-t))
    t = time.time()
    read_clinvar(args, db)
    print('ClinVar: {:.3f}s'.format(time.time()-t))
    # The following three calls are slow due to loading the whole genome
    t = time.time()
    read_ESP6500(args, db)
    print('ESP6500: {:.3f}s'.format(time.time()-t))
    t = time.time()
    read_HGVD(args, db)
    print('HGVD: {:.3f}s'.format(time.time()-t))
    t = time.time()
    read_GONL(args, db)
    print('GONL: {:.3f}s'.format(time.time()-t))
    #read_params(args)
    parse_vcfs(args, db)

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    #parser.add_argument("arg", help="Required positional argument")

    # Optional argument flag which defaults to False
    parser.add_argument('-f', '--flag', action="store_true", default=False)

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-o", "--outfile", action="store", help="Full output", default="annotationreport.txt")
    parser.add_argument("-r", "--reportfile", action="store", help="Smaller report file", default="annotationreport_small.txt")
    parser.add_argument("-c", "--config", action="store", help="Config file", default="exome.cfg")
    parser.add_argument("-g", "--genes", action="store", help="List of genes", default="genes.list")
    parser.add_argument("-s", "--samples", action="store", help="List of samples", default="samples.list")
    #parser.add_argument("-i", "--info", action="store", help="List of INFO fields to keep", default="info_fields.txt")
    parser.add_argument("-l", "--logfile", action="store", help="Logfile", default="logfile.log")

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

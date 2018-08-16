#!/usr/bin/env python3
"""
Module Docstring
"""

__author__ = "Harald Grove"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import time
import sys


def check_file(outputfile, qc_file):
    try:
        open(qc_file, "r")
    except:
        sys.stderr.write("ERROR: File not found [{}]\n".format(qc_file))
        sys.exit(1)
    with open(qc_file, "r") as fin, open(outputfile, "a") as fout:
        fout.write(
            "Sample\t \
                    Total Sequences\t \
                    Per base sequence quality\t \
                    Per tile sequence quality\t \
                    Per sequence quality scores\t \
                    Per base sequence content\t \
                    Per sequence GC content\t \
                    Per base N content\t \
                    Sequence Length Distribution\t \
                    Sequence Duplication Levels\t \
                    Overrepresented sequences\t \
                    Adapter Content\t \
                    Kmer Content\n"
        )
        overrep = False
        issues = [0, 0]  # all, serious
        for line in fin:
            l = line.strip().split("\t")
            if l[0] == "Filename":
                print("Filename: {}".format(qc_file.split("/")[-2]))
                fout.write("{}".format(qc_file))
            elif l[0] == "Total Sequences":
                print("Total sequences:\t{}".format(l[1]))
                fout.write("\t{}".format(l[1]))
            elif l[0] == ">>Per base sequence quality":
                if l[1] == "fail":
                    print("Per base sequence quality:\tFAIL")
                    issues[0] += 1
                    issues[1] += 1
                fout.write("\t{}".format(l[1]))
            elif l[0] == ">>Per tile sequence quality":
                if l[1] == "fail":
                    print("Per tile sequence quality:\tFAIL")
                    issues[0] += 1
                    issues[1] += 1
                fout.write("\t{}".format(l[1]))
            elif l[0] == ">>Per sequence quality scores":
                if l[1] == "fail":
                    print("Per sequence quality scores:\tFAIL")
                    issues[0] += 1
                    issues[1] += 1
                fout.write("\t{}".format(l[1]))
            elif l[0] == ">>Per base sequence content":
                if l[1] == "fail":
                    print("Per base sequence content:\tFAIL")
                    issues[0] += 1
                    issues[1] += 1
                fout.write("\t{}".format(l[1]))
            elif l[0] == ">>Per sequence GC content":
                if l[1] == "fail":
                    print("Per base sequence GC content:\tFAIL")
                    issues[0] += 1
                    issues[1] += 0
                fout.write("\t{}".format(l[1]))
            elif l[0] == ">>Per base N content":
                if l[1] == "fail":
                    print("Per base N content:\tFAIL")
                    issues[0] += 1
                    issues[1] += 1
                fout.write("\t{}".format(l[1]))
            elif l[0] == ">>Sequence Length Distribution":
                if l[1] == "fail":
                    print("Sequence length distribution:\tFAIL")
                    issues[0] += 1
                    issues[1] += 1
                fout.write("\t{}".format(l[1]))
            elif l[0] == ">>Sequence Duplication Levels":
                if l[1] == "fail":
                    print("Sequence duplication levels:\tFAIL")
                    issues[0] += 1
                    issues[1] += 1
                fout.write("\t{}".format(l[1]))
            elif l[0] == "#Total Deduplicated Percentage":
                fout.write("({})".format(l[1][0:5]))
            elif l[0] == ">>Overrepresented sequences":
                if l[1] == "fail":
                    print("Overrepresented sequences:\tFAIL")
                    issues[0] += 1
                    issues[1] += 1
                fout.write("\t{}".format(l[1]))
                overrep = True
            elif l[0][0:4] == "NNNN" and overrep:
                fout.write("({})".format(l[2][0:5]))
                overrep = False
            elif l[0] == ">>Adapter Content":
                if l[1] == "fail":
                    print("Adapter content:\tFAIL")
                    issues[0] += 1
                    issues[1] += 1
                fout.write("\t{}".format(l[1]))
            elif l[0] == ">>Kmer Content":
                if l[1] == "fail":
                    print("Kmer content:\tFAIL")
                    issues[0] += 1
                    issues[1] += 0
                fout.write("\t{}".format(l[1]))
                fout.write("\n")
    print(
        "QC finished. There were {} issues, {} of which were serious.".format(
            issues[0], issues[1]
        )
    )


def main(args):
    """ Main entry point of the app """
    for qc_file in args.infile:
        check_file(args.outfile, qc_file)
    if args.log:
        with open("README.txt", "a") as fout:
            fout.write("[{}]\t[{}]\n".format(time.asctime(), " ".join(sys.argv)))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("infile", nargs="*", help="Input file(s)")

    # Optional argument flag which defaults to False
    parser.add_argument("-o", "--outfile", help="Output report file")
    parser.add_argument(
        "-l",
        "--log",
        action="store_true",
        default=False,
        help="Save command to 'README.txt'",
    )

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-n", "--name", action="store", dest="name")

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

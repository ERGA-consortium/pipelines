#!/usr/bin/env python3

import argparse
import sys
import re

parser = argparse.ArgumentParser(
    description='Analyse gff output file for ANNOTATO with respect of the number of exons per mRNA')
parser.add_argument('-f', '--gff_file', type=str, required=True,
                    help="Input gff file")
args = parser.parse_args()

# read gff file into a dictionary that has transcript name as key and number of exons as value
mRNA_dict = {}
tRNA_dict = {}
try:
    with open(args.gff_file, "r") as gff_file:
        for line in gff_file:
            if re.search(r'\tmRNA\t', line):
                if re.search(r'ID=([^;]+)', line).group(1) not in mRNA_dict:
                    mRNA_dict[re.search(r'ID=([^;]+)', line).group(1)] = 0
            if re.search(r'\ttRNA\t', line):
                if re.search(r'ID=([^;]+)', line).group(1) not in tRNA_dict:
                    tRNA_dict[re.search(r'ID=([^;]+)', line).group(1)] = 0
            if re.search(r'exon', line):
                key = re.search(r'ID=([^;]+)', line).group(1).split(".exon")[0]
                if key in mRNA_dict:
                    mRNA_dict[key] += 1
                else:
                    tRNA_dict[key] += 1
except IOError:
    print("Could not read file:", args.gff_file)
    sys.exit(1)

# output some information about the dictionary
reports = [mRNA_dict]
if len(tRNA_dict) > 0:
    reports.append(tRNA_dict)

for info_dict in reports:
    if info_dict == mRNA_dict:
        print("INFORMATION REGARDING mRNA")
    else:
        print("INFORMATION REGARDING tRNA")
    print("Number of transcripts:", len(info_dict))
    print("Largest number of exons in all transcripts:", max(info_dict.values()))
    print("Monoexonic transcripts:", list(info_dict.values()).count(1))
    print("Multiexonic transcripts:", len(info_dict) - list(info_dict.values()).count(1))
    if (len(info_dict) - list(info_dict.values()).count(1)) == 0:
        print("No multiexonic transcripts, unable to calculate Mono:Mult Ratio")
    else:
        print("Mono:Mult Ratio:", round(list(info_dict.values()).count(1) / (len(info_dict) - list(info_dict.values()).count(1)), 2))

    # print an ASCII boxplot of data in tx_dict
    # get the largest number of exons in all transcripts
    max_exons = max(info_dict.values())
    # get the smallest number of exons in all transcripts
    min_exons = min(info_dict.values())
    # get the 25% quartile of values in tx_dict
    q25 = sorted(info_dict.values())[int(len(info_dict) * 0.25)]
    # get the 50% quartile of values in tx_dict
    q50 = sorted(info_dict.values())[int(len(info_dict) * 0.5)]
    # get the 75% quartile of values in tx_dict
    q75 = sorted(info_dict.values())[int(len(info_dict) * 0.75)]


    # print the boxplot
    print("Boxplot of number of exons per transcript:")
    print("Min:", min_exons)
    print("25%:", q25)
    print("50%:", q50)
    print("75%:", q75)
    print("Max:", max_exons)
    print("Mean:", sum(info_dict.values())/len(info_dict))
    print("="*50)

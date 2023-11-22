import os
import argparse
import pandas as pd

def main(lookupTab, inputGff, outputGff):
    df_lookup = pd.read_csv(lookupTab, sep = "\t", names=["original", "scaffold"])
    df_inGff = pd.read_csv(inputGff, sep = "\t", names=["scaffold", "tool", "type", "start", "end", "score", "strand", "phase", "attributes"])
    df_merged = df_inGff.merge(df_lookup, on="scaffold")[["original", "tool", "type", "start", "end", "score", "strand", "phase", "attributes"]]
    df_merged.to_csv(outputGff, sep="\t", header=False, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename gff file contig back to the original one")
    parser.add_argument("-l", "--lookup", help="Lookup table for converting contig name")
    parser.add_argument("-i", "--inGff", help="Gff file to be converted")
    parser.add_argument("-o", "--outGff", help="Name of the output Gff file")
    args = parser.parse_args()
    main(os.path.abspath(args.lookup), os.path.abspath(args.inGff), os.path.abspath(args.outGff))
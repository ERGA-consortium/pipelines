#!/usr/bin/env python
import argparse
import sys
import os
import re
import textwrap

# external dependencies
# import pysam
from Bio import SeqIO
from Bio.Seq import Seq


# -----------------------------------------
class GFF_line:
    def __init__(self, sequence_name, source, feature, start, end, score, strand, frame, attribute):
        self.sequence_name = sequence_name
        self.source = source
        self.feature = feature
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame
        attlist = re.split("[=;]", attribute)
        self.attribute = {x: y for x,y in zip(attlist[0::2], attlist[1::2])}

    def print_line(self):
        last_column = [key + '=' + value for key, value in self.attribute.items()]
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.sequence_name, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, ';'.join(last_column)))
    
    def ofprint_line(self, fileName):
        last_column = [key + '=' + value for key, value in self.attribute.items()]
        fileName.write(self.sequence_name+"\t"+self.source+"\t"+self.feature+"\t"+str(self.start)+"\t" +
                       str(self.end)+"\t"+str(self.score)+"\t"+self.strand+"\t"+self.frame+"\t" +
                       ';'.join(last_column)+"\n")

    def get_seqname(self):
        return str(self.sequence_name)

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end
    
    def get_strand(self):
        return str(self.strand)

    def get_feature(self):
        return str(self.feature)

    def has_id(self): 
        if "ID" in self.attribute:
            return True
        else:
            return False

    def has_parent(self): 
        if "Parent" in self.attribute:
            return True
        else:
            return False

    def get_id(self): 
        if self.has_id():
            return str(self.attribute["ID"])
        else:
            return None

    def get_parent(self): 
        if self.has_parent():
            return str(self.attribute["Parent"])
        else:
            return None # il faudrait lever une exception et quitter prog

#-----------------------------------------

def printseq(fh, nom, struct, record_dict, to_translate, gencode):
    if len(struct) > 0:
        frag = ''
        for segment in struct:
            subseq = record_dict[segment.get_seqname()].seq[segment.get_start()-1:segment.get_end()]
            frag += str(subseq)
            #frag = frag + seqfile.fetch(segment.get_seqname(), segment.get_start()-1, segment.get_end()) # 0-based for start
            #frag = frag + str(fetch_sequence(fasta, segment.get_seqname(), segment.get_start(), segment.get_end()))

        seq = Seq(frag)
        if struct[0].get_strand() == "-":
            seq = seq.reverse_complement() # try

        if to_translate: # si pas multiple de 3, exception --> try
            seq = seq.translate(table=gencode) # mettre un warning si des stop en phase ?

        wrapper = textwrap.TextWrapper(width=60)
        sublines = wrapper.fill(text=str(seq))
        fh.write(">"+nom+"\n")
        fh.write(sublines+"\n")

#-----------------------------------------
def __main__() :
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", action="store", dest="gff_file", help="GFF3 file containing locations to extract. It have to be sorted by identifier (ID tag's value in the last field). Default = STDIN", default=sys.stdin, required=False)
    parser.add_argument("--seq", action="store", dest="seq_file", help="Fasta file that contain sequences describe in the GFF3 file", required=True)
    parser.add_argument("--out", action="store", dest="out_file", help="Extracted sequence from GFF3 file (Fasta). Default = STDOUT", default=sys.stdout, required=False)
    parser.add_argument("--concat", action="store_true", dest="concat_line", help="Concatenate lines by ID and tagged as type defined by --type option", required=False)
    parser.add_argument("--type", action="store", dest="type_to_print", help="Type of line to print (fied 2 of GFF). Default = CDS", default="CDS", required=False)
    parser.add_argument("--tr", action="store_true", dest="to_translate", help="Translate extracted sequences using the genetic code given by --gc option.", required=False)
    parser.add_argument("--gc", action="store", dest="genetic_code", help="Genetic code used to translate sequences (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). Default = 1 (standard code)", default="1", required=False)
    args = parser.parse_args()

    # open gff input file
    if args.gff_file != sys.stdin : # file
        fileIn = open(args.gff_file, 'r')
    else : # STDIN
        fileIn = args.gff_file

    # open faidx file
    #try:
        #seqbk = pysam.FastaFile(args.seq_file) # si erreur: retourne ValueError ou IOError
        #seqbk = SeqIO.parse(args.seq_file, "fasta") 
    #except Exception as E:
        #print(E)
        #raise E
        #exit(1) 

    # indexing fasta file with Biopython
    try:
        record_dict = SeqIO.index(args.seq_file, "fasta")
    except Exception as E:
        print(E)
        raise E
        exit(1)

    # output fasta file
    if args.out_file != sys.stdout : # file
        fileOut = open(args.out_file, 'w') # mettre dans un try
    else : # STDOUT
        fileOut = args.out_file

    # read gff input file
    current_id = ''
    struct_to_extract = list()
    gff_line = fileIn.readline().strip()
    while len(gff_line) > 1:
        gff_line = gff_line.split('\t')
        while "#" in gff_line[0]: # delete comment lines
            gff_line = fileIn.readline().strip().split("\t")
 
        record = GFF_line(*gff_line)
        
        if record.get_feature() != args.type_to_print: # next line if the type is not the choosen one (--type)
            gff_line = fileIn.readline().strip() # read new gff line
            continue

        # select id
        name = ''
        if record.has_parent():
            name = record.get_parent()
        else:
            name = record.get_id()

        if args.concat_line:
            # si nouveau id, alors print seq, flush de la stuct et mettre a jour id
            if current_id != name:
                printseq(fileOut, current_id, struct_to_extract, record_dict, args.to_translate, args.genetic_code) # dans la fonction, tester si la structure est vide ou non
                struct_to_extract = list()
                current_id = name
            struct_to_extract.append(record)
 
        else:
            struct_to_extract.append(record)
            printseq(fileOut, name, struct_to_extract, record_dict, args.to_translate, args.genetic_code) # dans la fonction, tester si la structure est vide ou non
            struct_to_extract = list()

        gff_line = fileIn.readline().strip() # read gff new line
    
    if args.concat_line:
        printseq(fileOut, current_id, struct_to_extract, record_dict, args.to_translate, args.genetic_code)

    fileOut.close()
    fileIn.close()
    #seqbk.close()
    record_dict.close()

__main__()

#!/bin/bash

fasta1=$1
fasta2=$2
type=$3

if [[ $type = "uniprot" ]]
then
zcat $fasta1 | awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
       (FNR==1){next}
       { name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
       !(seen[seq]++){ print ">" $0 }' - $fasta2 | gzip > combined.fasta.gz
else
awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
       (FNR==1){next}
       { name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
       !(seen[seq]++){ print ">" $0 }' $fasta1 $fasta2 > combined.fasta
fi

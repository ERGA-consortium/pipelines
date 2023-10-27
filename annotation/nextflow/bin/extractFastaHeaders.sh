#!/bin/bash

fasta1=$1
fasta2=$2
out=$3

paste -d '\t' <(grep '^>' $fasta1 | sed 's/>//') <(grep '^>' $fasta2 | sed 's/>//') > $out
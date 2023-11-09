#!/bin/bash

infasta=$1
outfasta=$2

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $infasta | \
	awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' | \
	sort -k1,1nr | cut -f 2- | tr "\t" "\n" > $outfasta
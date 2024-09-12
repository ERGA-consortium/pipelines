#!/bin/bash

# Check if two arguments are provided
if [ $# -ne 2 ]; then
  echo "Usage: $0 <file1.tsv> <file2.tsv>"
  exit 1
fi

# Define files
file1="$1"
file2="$2"

# Use awk for processing
awk -v file1="$file1" '
FNR == NR {
  key1[$1] = $2
  key2[$2] = $1
}

FNR != NR {
  if ($2 == key2[$1] && $1 == key1[$2]) {
    print $0
  }
}
' "$file1" "$file2" >> matching_orthologs.txt
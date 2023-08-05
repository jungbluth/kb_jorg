#!/usr/bin/env bash

binID=${1}
outfile=${2}

awk '{print $1}' list.txt | tail -n+2 > temp.list && mv temp.list list.txt
seqtk subseq ${binID}_MIRA.fasta list.txt > "$outfile"

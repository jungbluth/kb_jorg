#!/usr/bin/env bash

binID=${1}

seqtk comp ${binID}.out.fasta | sort -k 2,2nr | awk '{printf "%-40s %8d %7.2f %8d\n", $1, $2, 100.0 * ($4 + $5)/($3 + $4 + $5 + $6), sum += $2}' >> iterations.txt


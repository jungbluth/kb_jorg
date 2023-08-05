#!/usr/bin/env bash

min_coverage=${1}
binID=${2}

awk -v coverage_min="$min_coverage" '{if($2 >= 1000 && $6 >= coverage_min || NR == 1) print}' "$binID"_assembly/"$binID"_d_info/"$binID"_info_contigstats.txt > list.txt

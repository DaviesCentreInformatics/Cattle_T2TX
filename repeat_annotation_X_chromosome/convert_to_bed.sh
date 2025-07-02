#!/bin/bash

trf_output="trf_chrX_2000.txt"
output_bed="trf_chrX_2000.bed"

awk 'BEGIN {OFS="\t"} {
    if ($1 ~ /^[0-9]+$/ && NF >= 14) {
        chromosome="X:38000000-50000000"
        start=$1 + 1
        end=$2 + 1
        name=$14
        print chromosome, start, end, name
    }
}' $trf_output > $output_bed

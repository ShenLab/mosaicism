#!/bin/bash
## USAGE: sh runfilters_sam.sh <FILTERED DENOVOS CLEANED>
## Purpose: Preprocessing filter to handle the following:
##			mappability (300bp), segmental duplication, LCR, strand bias, etc
## Dependencies: LCR, mapapbility 300bp, and segmental duplication tracks from UCSC genome browser

file=$1
base=${file%.*}
echo $base

## Generate bed file from candidates file
echo "Generating candidates bed file"
awk -F '\t' '{if($1!="id") print "chr"$3"\t"$4"\t"$4}' "$file" | sort -k1,1 -k2,2n > $base".bed"
echo "done - see "$base".bed"
echo "====================================="
echo " "
## run repeat filters
echo "Running intersectBed"
intersectBed -wa -wb -a "$base".bed -b ~/reference_files/LCR.addchr.bed ~/reference_files/mappability_lt1_300mer.bed  ~/reference_files/hg19_segdup.bed -filenames > "rfilt."$file
echo "done - see rfile."$file
echo "====================================="
echo " "
## run pileup
echo "Running pileup.py"
cp "$file" "$base".dp4.txt
## run strand_bias.R
echo "Running strand_bias.R"
Rscript ~/mosaic/strand_bias.R "$base".dp4.txt "$base"
echo "done - see "$base".sb.txt"
## run parsefilters.py
echo "Running parsefilters.py"
python ~/mosaic/parsefilters.py "$file" rfilt."$file" "$file" "$base".sb.txt> "$base".filtered.txt
echo "done - see "$base".filtered.txt"
wc -l "$base".filtered.txt
echo "====================================="
echo " "
## done message
echo "DONE ... see $base".filtered.txt""
python ~/mosaic/getfiltct.py "$base".filtered.txt

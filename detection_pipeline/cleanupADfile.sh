#!/bin/bash
## USAGE: sh cleanupAFfile.sh <ADfile> <complex variant distance cutoff>
## Purpose: remove complex variants (as defined by distance cutoff) = variants found within $distance bp of each other

file=$1
base=${file%.*}
distance=$2
wc -l "$file"
python ~/mosaic/remove_complex.py "$file" "$distance" > "$base".clean.txt
wc -l "$base".clean.txt


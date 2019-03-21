#!/bin/bash

VAR=$1
OUTANNO=$2
OUTPREFIX=$3

if [[ ! "$VAR" ]] || [[ ! "$OUTANNO" ]] || [[ ! "$OUTPREFIX" ]]; then
 echo "Missing/Incorrect required arguments"
 echo "Required: <variants file> <annotation output filename> <output file prefix>"
 exit
fi

echo "# Input File: $VAR"
echo "# Annotation Output File: $OUTANNO"
echo "# Output File Prefix: $OUTPREFIX"

echo "Annotating Variants ..."
python ~/mosaic/production_code_06052018/pcgc_annotate_variants.py $VAR > $OUTANNO
echo "done"

echo "Running Enrichment Analysis ... "
Rscript ~/mosaic/production_code_06052018/pcgc_enrichment_analysis.CHD_risk.R $OUTANNO $OUTPREFIX
echo "done"

echo "Running Case-Only Analysis ... "
Rscript ~/mosaic/production_code_06052018/pcgc_case_only_analysis.CHD_risk.R $OUTANNO $OUTPREFIX
echo "done"

echo "Running VAF vs. vclass Analysis ... "
Rscript ~/mosaic/production_code_06052018/VAF_vs_vclass.PCGC.CHD_risk.R $OUTANNO $OUTPREFIX
#Rscript ~/mosaic/production_code_06052018/pcgc_vaf_analysis.CHD_risk.R $OUTANNO $OUTPREFIX
echo "done"

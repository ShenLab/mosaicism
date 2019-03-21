#!/bin/bash

# This script takes a parsed denovo variants file (samtools or GATK) as input
# and filters the variants before running mosaic scoring and detection

#############################################################################

# usage message
usage="
	## Usage:
	-v (required) - .txt file containing unfiltered denovos
	-c (required) - for now, either 'samtools' or 'gatk'
	-H (flag) - echo this message and exit
"

# parse arguments
while getopts v:c:H opt; do 
	case "$opt" in 
		v) NOFILT="$OPTARG";;
		c) CALLER="$OPTARG";;
		H) echo "$usage"; exit;;
		\?) echo "$usage"; exit;;
		:) echo "$uage"; exit;;
	esac
done 

# check for arguments
if [[ -z $NOFILT ]] | [[ -z $CALLER ]]; then
	echo "$usage"; exit
fi

## get filename bases
BASE0=`basename $NOFILT | sed 's/.no_filt.txt//g'` # strip .no_filt.txt
BASE1=`basename $NOFILT | sed 's/.txt//g'` # strip .txt

## set directory containing all necessary scripts
SDIR=/home/local/ARCS/alh2194/mosaic/code_update_11Jan2019

#===============================================================================
##
## STEP (1) Run ANNOVAR on unfiltered denovos
##
echo "
##### (STEP 1) Running ANNOVAR ... "

#sh $SDIR/run_annovar.sh $NOFILT

## check output
ANNOUT=`basename $NOFILT | sed 's/.txt/.avinput.hg19_multianno.txt/g'`

if [[ !(-s $ANNOUT) ]]; then # check that $ANNOUT exists and has nonzero filesize
	echo "## ERROR: ANNOVAR output empty"
	exit
fi

echo "##### (STEP 1) Running ANNOVAR ... done
##### see -- $ANNOUT"
#===============================================================================

#===============================================================================
##
## STEP (2) Fix new line issue
##
echo "
##### (STEP 2) Fixing newline ... "

FIX_NL=$BASE1'.fixnewline.txt'

awk -F '\t' '{if(NF!=0) print $0}' $NOFILT > $FIX_NL

echo "##### (STEP 2) Fixing newline ... done
##### see -- $FIX_NL
"
#===============================================================================

#===============================================================================
##
## STEP (3) Filter for rare (MAF<1e-4) exonic SNVs
##
echo "##### (STEP 3) Filtering for rare exonic SNVs ... "

FILT=$BASE0'.filt.txt'

python $SDIR/filter_denovo.py $FIX_NL $ANNOUT > $FILT

echo "##### (STEP 3) Filtering for rare exonic SNVs ... done
##### see -- $FILT
"
#===============================================================================

#===============================================================================
##
## STEP (4) Calculate cohort allele frequency using tab-separated list of sites:carriers
##
echo "
##### (STEP 4) Calculating cohort allele frequency ... "

CAF=$BASE0'.filt.cohortAF.txt'
CAF_REF='/home/local/ARCS/alh2194/mosaic/code_update_11Jan2019/pcgc_vcf_cohort_ct.txt'
COHORTSIZE='2530'

python $SDIR/check_cohortAF.py $CAF_REF $FILT $COHORTSIZE > $CAF

echo "##### (STEP 4) Cohort AF ... done
##### see -- $CAF
"
#===============================================================================

#===============================================================================
##
## STEP (5) Parse only samples found in Jin 2017 publication
##
echo "
##### (STEP 5) Parsing Jin 2017 samples ... "

JINF=$BASE0'.filt.cohortAF.jin.txt'
JIN_IDS='/home/local/ARCS/alh2194/mosaic/pcgc_pipeline/raw_denovo/jin_2530_ids.txt'

python $SDIR/parse_jin.py $CAF $JIN_IDS > $JINF

echo "##### (STEP 5) Parsing Jin 2017 samples ... done
##### see -- $JINF
"
#===============================================================================

#===============================================================================
##
## STEP (6) Remove variant clusters within $DIST bp
##
echo "
##### (STEP 6) Removing variant clusters ... "

DIST='10'
NOCLUSTF=$BASE0'.filt.cohortAF.jin.no_clust_'$DIST'.txt'

python $SDIR/remove_vclust.py $JINF $DIST > $NOCLUSTF

echo "##### (STEP 6) Removing variant clusters ... done
##### see -- $NOCLUSTF
"
#===============================================================================

#===============================================================================
##
## STEP (6) Filter out variants located in repeat regions
##
echo "
##### (STEP 7) Filter out variants in repeat regions ... "

LCR='/home/local/ARCS/alh2194/reference_files/LCR.addchr.bed'
MAP='/home/local/ARCS/alh2194/reference_files/mappability_lt1_300mer.bed'
SEGDUP='/home/local/ARCS/alh2194/reference_files/hg19_segdup.bed'
FILT2F=$BASE0'.filt.cohortAF.jin.no_clust_'$DIST'.filtered.txt'

if [[ $CALLER == 'samtools' ]]; then
	#python $SDIR/runfilters_sam.sh $NOCLUSTF $LCR $MAP $SEGDUP > $FILT2F
	# run bedtools intersect
	echo "	(Running bedtools intersect)"
	awk -F '\t' '{if($1!="id") print "chr"$3"\t"$4"\t"$4}' $NOCLUSTF | sort -k1,1 -k2,2n > 'rfilt.'$BASE0'.bed'
	intersectBed -wa -wb -a 'rfilt.'$BASE0'.bed' -b $LCR $MAP $SEGDUP -filenames > 'rfilt.'$BASE0'.txt'
	# run strand bias filter
	echo "	(Calculating strand bias)"
	Rscript $SDIR/strand_bias.R $NOCLUSTF 'rfilt.'$BASE0
	# parse filter outputs
	echo "	(Parsing filtering results)"
	python $SDIR/parsefilters.py $NOCLUSTF 'rfilt.'$BASE0'.txt' 'rfilt.'$BASE0'.sb.txt' > $FILT2F
	#rm -rf 'rfilt.'$BASE0'*'
	echo "	(Filtering results)
	"
	python $SDIR/getfiltct.py $FILT2F
fi

if [[ $CALLER == 'gatk' ]]; then
	python $SDIR/runfilters_gatk.sh $NOCLUSTF > $FILT2F
	## INCOMPLETE - NEED TO UPDATE
fi

if [[ !($CALLER == 'samtools') && !($CALLER == 'gatk') ]]; then
	echo "Please check caller argument - only 'samtools' and 'gatk' accepted"
fi

echo "
##### (STEP 7) Filter out variants in repeat regions ... done
##### see -- $FILT2F
"
#===============================================================================

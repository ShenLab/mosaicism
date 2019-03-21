#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -q shenlab.q
#$ -N SamtoolsVariantCalling 
#$ -l h_rt=128:00:00
#$ -l h_vmem=4G
#$ -cwd


# This script takes a list of gVCF files generated by the HaplotypeCaller (filename must end ".list") and performs the multi-sample joint aggregation step and merges the records together.
# The job should be called as an array. The the variants calling will then be split into X jobs, where X is the number of jobs in the array. This should be larger for bigger jobs (more samples).
#    InpFil - (required) - List of gVCF files. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file- only required if calling pipeline
#    VcfNam - (optional) - A name for the analysis - to be used for naming output files. Will be derived from input filename if not provided
#    LogFil - (optional) - File for logging progress
#    JobNum - (optional) - The number of jobs in array; normally this is derived automatically from the -t argument of qsub command. However, if specific sections of an array job fail this can be used to rerun them, ensuring that the target file is split in the same way.
#    Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#    Flag - B - BadET - prevent GATK from phoning home
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $DBSNP - dbSNP vcf from GATK
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# GATK <https://www.broadinstitute.org/gatk/> <https://www.broadinstitute.org/gatk/download>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

## modified from Jiayao Wang's original script

###############################################################

#set default arguments
usage="
-t 1-<X> ExmVC.7b.SamtoolsCall.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

-i (required) - List of BAM files. List file name must end \".list\"
-r (required) - shell file containing variables with locations of reference files and resource directories
-t (required) - Exome capture kit targets or other genomic intervals bed file. 
-n (optional) - Analysis/output VCF name - will be derived from input filename if not provided; only used if calling pipeline
-l (optional) - Log file
-j (optional) - number of jobs in array; normally this is derived automatically from the -t argument of qsub command; this can be used to rerun individual sections that failed
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -X (flag) - Do not run Variant Quality Score Recalibration - only if calling pipeline
     -B (flag) - Prevent GATK from phoning home
     -H (flag) - echo this message and exit
"

AllowMisencoded="false"
PipeLine="false"
NoRecal="false"
BadET="false"

while getopts i:r:a:t:n:l:j:PXBH opt; do
    case "$opt" in

		i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";;
        a) ArrNum="$OPTARG";; 
        t) TgtBed="$OPTARG";;
        n) VcfNam="$OPTARG";;
        j) JobNum="$OPTARG";;
        l) LogFil="$OPTARG";;
        P) PipeLine="true";;
        X) NoRecal="true";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]] || [[ -z "$TgtBed" ]]; then
 echo "Missing/Incorrect required arguments"
 echo "provided arguments: -i $InpFil -r $RefFil -t $TgtBed"
 echo "usage: $usage"
 exit
fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"


#Set local Variables
funcGetTargetFile #If the target file has been specified using a code, get the full path from the exported variable

GatkLog=$VcfNam.VCF.samtoolslog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$VcfNam.VCF.temp.log #temporary log file
TmpDir=$VcfNam.VCF.tempdir; mkdir -p $TmpDir #temporary directory
TgtFil=$TgtBed
VcfFil=$VcfNam.vcf
VcfAnnFil=$VcfNam.ann.vcf
VcfLeftAlnFil=$VcfNam.LA.vcf

##Run Joint Variant Calling
StepName="Joint call Variants with Samtools/bcftools"
StepCmd="samtools mpileup -C50 -u -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -q 20 -d 10000 -f $REF -l $TgtFil  -b $InpFil 2>$GatkLog | bcftools call - -c -v -p 1.1 2>$GatkLog > $VcfFil" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

##Write completion log
touch $PrgDir/$PrgFil

#End Log
funcWriteEndLog

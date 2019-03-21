#!/bin/bash
VcfFil=$1
if [[ ! "$VcfFil" ]] ; then
 echo "Missing/Incorrect required arguments"
 echo "Required: <variants file>"
 exit
fi

## scripts and arguments
TABLE_ANNOVAR="/home/local/users/jw/software_packages/annovar/table_annovar.pl"
ANNHDB="/home/local/users/jw/resources/ANNOVAR_DB/"
BUILD="hg19"

## convert variants to .avinput format
AVINPUT=`basename $VcfFil | sed 's/.txt/.avinput/g'`
awk -F '\t' '{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$3}' $VcfFil > $AVINPUT

## run table_annovar
CMD="$TABLE_ANNOVAR $AVINPUT $ANNHDB --buildver $BUILD --remove -protocol refGene,gnomad_exome,gnomad_genome,1000g2015aug_all,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_afr,1000g2015aug_sas,exac03,dbnsfp33a,cadd13gt10,cosmic70,genomicSuperDups,mcap,revel,avsnp147,mvp -operation g,f,f,f,f,f,f,f,f,f,f,f,f,r,f,f,f,f -otherinfo  -nastring . "

$CMD 

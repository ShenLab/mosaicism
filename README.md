# Mosaicism in Congenital Heart Disease

The aims of this project are to (1) develop a method to detect mosaic (post-zygotic) SNVs from exome sequencing data and (2) estimate the contribution of mosaicism to Congenital Heart Disease.  

directories
```
variant_calling/ - contains scripts used to call de novo variants from BAM files
detection_pipeline/ - contains scripts used to QC de novo variants and call mosaic SNVs
analysis/ - contains scripts used to analyze detected mosaics
```

### Prerequisites

SAMtools (http://samtools.sourceforge.net/)
GATK (https://software.broadinstitute.org/gatk/download/)
ANNOVAR (http://annovar.openbioinformatics.org/en/latest/)

## Authors

* **Alexander Hsieh** - [alexanderhsieh](https://github.com/alexanderhsieh)

## Acknowledgments

Thanks to Jiayao and Hongjian for help with the variant calling
* **Jiayao Wang** - [explorerwjy](https://github.com/explorerwjy)
* **Hongjian Qi** - [7lagrange](https://github.com/7lagrange)

README.md template from PurpleBooth (https://gist.github.com/PurpleBooth)

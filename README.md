# Mosaicism in Congenital Heart Disease

The aims of this project:
1. to develop a method to detect mosaic (post-zygotic) SNVs from exome sequencing data 
2. to estimate the contribution of mosaicism to Congenital Heart Disease.  

Directories
* variant_calling/ - contains scripts used to call de novo variants from BAM files
* detection_pipeline/ - contains scripts used to QC de novo variants and call mosaic SNVs
* analysis/ - contains scripts used to analyze detected mosaics

Key Scripts
* detection_pipeline/generate_candidates.POST.R - detects mosaic SNVs from a list of de novo SNVs
* analysis/power_analysis.R - estimates mosaic detection power as a function of sample average depth

### Prerequisites

- SAMtools (http://samtools.sourceforge.net/)
- GATK (https://software.broadinstitute.org/gatk/download/)
- ANNOVAR (http://annovar.openbioinformatics.org/en/latest/)

## Authors

* **Alexander Hsieh** - [alexanderhsieh](https://github.com/alexanderhsieh)

## Acknowledgments

Thanks to Jiayao and Hongjian for help with the variant calling and IGV automation
* **Jiayao Wang** - [explorerwjy](https://github.com/explorerwjy)
* **Hongjian Qi** - [7lagrange](https://github.com/7lagrange)

README.md template from PurpleBooth (https://gist.github.com/PurpleBooth)

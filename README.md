# Mosaicism in Congenital Heart Disease

### Project Aims
1. to develop a method to detect mosaic (post-zygotic) SNVs from exome sequencing data 
2. to estimate the contribution of mosaicism to Congenital Heart Disease.  

### Directories
* variant_calling/ - contains scripts used to call de novo variants from BAM files
* detection_pipeline/ - contains scripts used to QC de novo variants and call mosaic SNVs
* analysis/ - contains scripts used to analyze detected mosaics

### Test Data
* ADfile.example_minimal.txt - example de novo callset for testing (contains only columns required for mosaic detection: id, chr, pos, ref, alt, refdp, altdp)
* ADfile.example_full_annotation.txt - contains additional columns and annotations used to QC/filter variants

### Key Scripts and Usage
* detection_pipeline/generate_candidates.POST.R - detects candidate mosaic SNVs by calculating posterior odds for each de novo SNV following annotation and QC
```
# Call mosaics from input de novo callset (ADfile.example_minimal.txt) with output prefix 'test', posterior odds cutoff '10', and cohort size '2400'

Rscript generate_candidates.POST.R ADfile.example_minimal.txt test 10 2400

# Output Files
* test.candidates.txt = contains mosaic variants with posterior odds > cutoff
* test.denovo.txt = contains all de novo variants in input passing filters, annotated with posterior odds score, etc
* test.all_denovos.txt = contains all de novo variants (including sites failing filters), annotated with exclusion criteria
* test.outlier_samples.txt = contains IDs of all outlier samples (based on cohort size)

# Output Plots
* test.EM.pdf = histogram of variant allele fraction of input variants, with EM-estimated mosaic fraction
* test.dp_vs_vaf.pdf = scatterplot of DP vs. VAF, colored by germline/mosaic
* test.vaf_vs_post.df = scatterplot of VAF vs. log-scaled posterior odds, colored by germline/mosiac
* test.fdr_min_nalt.pdf = scatterplot of DP vs. Nalt, with line indicating FDR-based minimum Nalt value
* test.overdispersion.pdf = overdispersion plot
* test.QQ.pdf = QQ plot
```
* analysis/power_analysis.R - estimates mosaic detection power as a function of sample average depth
```
# Plot detection power as a function of sample average depth with model parameters determined from generate_candidates.POST.R 
# Estimate the true frequency of mosaics with VAF>0.1 using 'test.denovo.txt' (with parameters LR cutoff = 41, Theta estimate = 59, cohortsize = 2400, sample avg depth = 80)

Rscript power_analysis.R test.denovo.txt test 41 59 2400 80

# Output Files
* test.pwsite.log.txt = data used for estimating detection power given variant site DP
* test.pwsamp.log.txt = data used to plot detection power as a function of sample avg DP
* test.vaf_pw_adj.log.txt = data used in adjusting mosaic counts

# Output Plots
* test.power_sample_dp.pdf = plot of detection power curves as a function of sample avg DP
* test.vaf_pw_adj.pdf = histogram of adjusted mosaic counts, raw mosaic rate, adjusted mosaic rate
```


### Prerequisites

- SAMtools (http://samtools.sourceforge.net/)
- GATK (https://software.broadinstitute.org/gatk/download/)
- ANNOVAR (http://annovar.openbioinformatics.org/en/latest/)
- IGV (https://software.broadinstitute.org/software/igv/)

## Authors

* **Alexander Hsieh** - [alexanderhsieh](https://github.com/alexanderhsieh)

## Acknowledgments

Thanks to Jiayao and Hongjian for help with the variant calling and IGV automation
* **Jiayao Wang** - [explorerwjy](https://github.com/explorerwjy)
* **Hongjian Qi** - [7lagrange](https://github.com/7lagrange)

README.md template from PurpleBooth (https://gist.github.com/PurpleBooth)

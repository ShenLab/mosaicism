## Usage: Rscript VAF_vs_vclass.PCGC.R <annotated denovos file> <output prefix>
## Purpose: For all mosaics, HHE mosaics, and non-HHE mosaics, plot VAF distribution of damaging sites vs. likely benign sites and run Mann-Whitney U Test
library(ggridges)
library(ggplot2)

## handle arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2){
  stop("Please provide <annotated variants file> <outfile prefix>")
}
fname <- toString(args[1])
# setwd('~/Desktop/pe/mosaic/pcgc_pipeline/freeze_26Feb2019/post2')
#fname <- 'union_post2.clean.anno.txt'
# fname <-'post2.denovo.anno.txt'
op <- toString(args[2])
print(paste("Input File", fname, sep=': '))
print(paste("Outfile prefix", op, sep=': '))

f <- read.table(fname, sep='\t', header=T, quote="")
if("VAF" %in% names(f)){
  f$vaf <- f$VAF
}
if(!('col' %in% names(f))){
  f$col <- 'red'
  f$vaf <- f$altdp/(f$refdp+f$altdp)
}


f.chd <- f[f$CHD_risk=='yes',]
f.nchd <- f[f$CHD_risk=='no',]
# syndromic
f.syn.chd <- f.chd[f.chd$NDD_diagnosis=="Yes" | f.chd$extracardiac=="Yes",]
f.syn.nchd <- f.nchd[f.nchd$NDD_diagnosis=="Yes" | f.nchd$extracardiac=="Yes",]
# isolated
f.non.chd <- f.chd[(f.chd$NDD_diagnosis=="No" & f.chd$extracardiac=="No"),]
f.non.nchd <- f.nchd[(f.nchd$NDD_diagnosis=="No" & f.nchd$extracardiac=="No"),]
# unknown
f.unk.chd <- f.chd[(f.chd$NDD_diagnosis == "Unknown" & f.chd$extracardiac=="No") | (f.chd$NDD_diagnosis=="No" & f.chd$extracardiac=="Unknown") | (f.chd$NDD_diagnosis=="Unknown" & f.chd$extracardiac=="Unknown"),]
f.unk.nchd <- f.nchd[(f.nchd$NDD_diagnosis == "Unknown" & f.nchd$extracardiac=="No") | (f.nchd$NDD_diagnosis=="No" & f.nchd$extracardiac=="Unknown") | (f.nchd$NDD_diagnosis=="Unknown" & f.nchd$extracardiac=="Unknown"),]

mo <- f[f$col=="red",]
mo.chd <- f.chd[f.chd$col=="red",]
mo.nchd <- f.nchd[f.nchd$col=="red",]

table(mo.chd$vclass)
table(mo.nchd$vclass)

## for union set
#f$vaf <- f$VAF
#mo <- f
#mo.chd <- mo[mo$CHD_risk=='yes',]
#mo.nchd <- mo[mo$CHD_risk=='no',]



## USE GGRIDGES TO PLOT DENSITY WITH ADDITIONAL INFORMATION
# https://www.r-bloggers.com/density-plot-with-ggplot/
## Function to take a set of variants and plot the VAF distribution of damaging sites against the VAF of likely benign sites
## Also run Mann-Whitney U Test to test significance of difference
plot_vaf_ridge <- function(tmp, set){
  tmp.mo <- tmp
  tmp.mo <- tmp.mo[tmp.mo$vclass %in% c('LGD', 'Dmis', 'Bmis', 'synonymous'),] 
  tmp.mo$grp <- rep('Likely Benign', nrow(tmp.mo))
  tmp.mo$grpcol <- rep('blue', nrow(tmp.mo))
  tmp.mo[tmp.mo$vclass %in% c('LGD', 'Dmis'),]$grp <- 'Damaging'
  tmp.mo[tmp.mo$vclass %in% c('LGD', 'Dmis'),]$grpcol <- 'red'
  wilcox.p <- wilcox.test(tmp.mo[tmp.mo$grp=="Damaging",]$vaf, tmp.mo[tmp.mo$grp=="Likely Benign",]$vaf)$p.value
  ks.p <- ks.test(tmp.mo[tmp.mo$grp=="Damaging",]$vaf, tmp.mo[tmp.mo$grp=="Likely Benign",]$vaf)$p.value
  title <- paste(set, paste(nrow(tmp), 'sites', sep=' '), sep='\n')
  p <- ggplot(tmp.mo, aes(x=vaf, y=grp, fill=grpcol, col=grpcol)) + 
    geom_density_ridges(jittered_points=TRUE, position=position_points_jitter(width=0.05, height=0), point_shape='|', point_size=5, alpha=0.4) + 
    annotate("text", x=0.08, y=0.75, label=paste('p-value', round(wilcox.p, 3), sep='='), size=6) +
    ggtitle(title) +
    xlab('Variant Allele Fraction (VAF)') +
    xlim(c(0, 0.4)) +
    theme_bw() +
    theme(plot.title = element_text(size=24), axis.title = element_text(size=20), axis.title.y = element_blank(), axis.text = element_text(size=20))
  p <- p+guides(fill=FALSE, col=FALSE)
  p <- p+scale_fill_manual(breaks = c("Likely Benign", "Damaging"), values=c("blue", "red")) + scale_color_manual(breaks = c("Likely Benign", "Damaging"), values=c("blue", "red"))
  return(p)
}

## GGRIDGES PLOTS
outname_all_mosaic <- paste(op, 'vaf_all_mosaic.RIDGE.pdf', sep='.')
pdf(outname_all_mosaic)
plot_vaf_ridge(mo, 'mosaics')
dev.off()

outname_chd_mosaic <- paste(op, 'vaf_chd_mosaic.RIDGE.pdf', sep='.')
pdf(outname_chd_mosaic)
plot_vaf_ridge(mo.chd, 'CHD-related Genes')
dev.off()

outname_nchd_mosaic <- paste(op, 'vaf_nchd_mosaic.RIDGE.pdf', sep='.')
pdf(outname_nchd_mosaic)
plot_vaf_ridge(mo.nchd, 'Other Genes')
dev.off()

data.frame(table(mo$vclass))
data.frame(table(mo.chd$vclass))
data.frame(table(mo.nchd$vclass))

print("Done !")
print(paste("output", c(outname_all_mosaic, outname_chd_mosaic, outname_nchd_mosaic), sep=' : '))


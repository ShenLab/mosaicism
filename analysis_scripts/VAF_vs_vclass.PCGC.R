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
op <- toString(args[2])
print(paste("Input File", fname, sep=': '))
print(paste("Outfile prefix", op, sep=': '))

f <- read.table(fname, sep='\t', header=T, quote="")

## temporarily remove Isolated CHD from analysis
f <- f[!(f$NDD_diagnosis=="No" & f$extracardiac=="No"),]

tmp.hhe <- f[!is.na(f$hexp),]
f.hhe <- tmp.hhe[tmp.hhe$hexp>=75,]
f.nhhe <- tmp.hhe[!(tmp.hhe$hexp>=75),]
# syndromic
f.syn.hhe <- f.hhe[f.hhe$NDD_diagnosis=="Yes" | f.hhe$extracardiac=="Yes",]
f.syn.nhhe <- f.nhhe[f.nhhe$NDD_diagnosis=="Yes" | f.nhhe$extracardiac=="Yes",]
# isolated
f.non.hhe <- f.hhe[(f.hhe$NDD_diagnosis=="No" & f.hhe$extracardiac=="No"),]
f.non.nhhe <- f.nhhe[(f.nhhe$NDD_diagnosis=="No" & f.nhhe$extracardiac=="No"),]

mo <- f[f$col=="red",]
mo.hhe <- f.hhe[f.hhe$col=="red",]
mo.nhhe <- f.nhhe[f.nhhe$col=="red",]

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

outname_hhe_mosaic <- paste(op, 'vaf_hhe_mosaic.RIDGE.pdf', sep='.')
pdf(outname_hhe_mosaic)
plot_vaf_ridge(mo.hhe, 'High Heart Expression')
dev.off()

outname_nhhe_mosaic <- paste(op, 'vaf_nhhe_mosaic.RIDGE.pdf', sep='.')
pdf(outname_nhhe_mosaic)
plot_vaf_ridge(mo.nhhe, 'non-High Heart Expression')
dev.off()

data.frame(table(mo$vclass))
data.frame(table(mo.hhe$vclass))
data.frame(table(mo.nhhe$vclass))

print("Done !")
print(paste("output", c(outname_all_mosaic, outname_hhe_mosaic, outname_nhhe_mosaic), sep=' : '))


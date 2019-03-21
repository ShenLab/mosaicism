## Usage: Rscript pcgc_enrichment_analysis.HHE_only.R <annotated denovos file> <output prefix>
## Purpose: Run Fisher's Exact Test using counts of Damaging Mosaics, Non-Damaging Mosaics, Damaging Germline, and Non-Damaging Germline sites
## Calculates enrichment for All samples, Syndromic CHD samples, Isolated CHD samples, and Unknown CHD samples
## Calculates enrichment for All genes and HHE genes
library(ggplot2)
library(gridExtra)

## handle arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2){
  stop("Please provide <annotated variants file> <outfile prefix>")
}
fname <- toString(args[1])
op <- toString(args[2])
print(paste("Input File", fname, sep=': '))
print(paste("Outfile prefix", op, sep=': '))

f <- read.table(fname, sep='\t', header=T, quote='"')
if("VAF" %in% names(f)){
  f$vaf <- f$VAF
}

f.syn <- f[f$NDD_diagnosis=="Yes" | f$extracardiac=="Yes",]
f.non <- f[(f$NDD_diagnosis=="No" & f$extracardiac=="No"),]
f.unk <- f[(f$NDD_diagnosis == "Unknown" & f$extracardiac=="No") | (f$NDD_diagnosis=="No" & f$extracardiac=="Unknown") | (f$NDD_diagnosis=="Unknown" & f$extracardiac=="Unknown"),]
print( c( paste('Total Variants', nrow(f), sep=' : '), paste('Syndromic', nrow(f.syn), sep=' : '), paste('Non-syndromic', nrow(f.non), sep=' : '), paste('NDD/CA Unknown', nrow(f.unk), sep=' : ') ) )
print( c( paste('Total Samples', length(unique(f$id)), sep=' : '), paste('Syndromic Samples', length(unique(f.syn$id)), sep=' : '), paste('Non-syndromic Samples', length(unique(f.non$id)), sep=' : '), paste('NDD/CA Unknown', length(unique(f.unk$id)), sep=' : ') ) )

summarize <- function(all, set){
  all.chdrisk <- all[all$CHD_risk == 'yes',]
  summarydf <- NULL
  cts <- c(nrow(all), nrow(all.chdrisk))
  summarydf <- rbind.data.frame(summarydf, cts)
  mocts <- c(nrow(all[all$col=="red",]), nrow(all.chdrisk[all.chdrisk$col=="red",]))
  summarydf <- rbind.data.frame(summarydf, mocts)
  germcts <- c(nrow(all[all$col=="black",]), nrow(all.chdrisk[all.chdrisk$col=="black",]))
  summarydf <- rbind.data.frame(summarydf, germcts)
  names(summarydf) <- c("All","CHD_risk")
  summarydf$src <- c('Total', 'Mosaic', 'Germline')
  summarydf$set <- c(set, set, set)
  return(summarydf)
}
all.summary <- summarize(f, 'all')
syn.summary <- summarize(f.syn, 'syndromic')
non.summary <- summarize(f.non, 'non-syndromic')
unk.summary <- summarize(f.unk, 'unknown')

print(all.summary)
print(syn.summary)
print(non.summary)
print(unk.summary)

## Parse variants into different gene sets
# all samples
f.chd <- f[f$CHD_risk=='yes',]
# syndromic
f.syn.chd <- f.chd[f.chd$NDD_diagnosis=="Yes" | f.chd$extracardiac=="Yes",]
# isolated
f.non.chd <- f.chd[f.chd$NDD_diagnosis=="No" & f.chd$extracardiac=="No",]
# unknown
f.unk.chd <- f.chd[(f.chd$NDD_diagnosis == "Unknown" & f.chd$extracardiac=="No") | (f.chd$NDD_diagnosis=="No" & f.chd$extracardiac=="Unknown") | (f.chd$NDD_diagnosis=="Unknown" & f.chd$extracardiac=="Unknown"),]

### RESUME HERE ###
### NEED TO COMBINE DFs for ALL, SYN, ISO
m_v_d <- function(v){
  mo <- v[(v$col=="red"),]
  no <- v[(v$col=="black"),]
  v.m.d <- mo[mo$vclass=="LGD" | mo$vclass=="Dmis",]
  v.m.nd <- mo[!(mo$vclass=="LGD" | mo$vclass=="Dmis"),]
  v.n.d <- no[no$vclass=="LGD" | no$vclass=="Dmis",]
  v.n.nd <- no[!(no$vclass=="LGD" | no$vclass=="Dmis"),]
  cont <- matrix(c(nrow(v.m.d), nrow(v.n.d), nrow(v.m.nd), nrow(v.n.nd)), nrow=2)
  fish <- fisher.test(cont)
  out <- c(nrow(no), nrow(mo), nrow(v.m.d), nrow(v.n.d), nrow(v.m.nd), nrow(v.n.nd), unname(fish$estimate), fish$p.value, fish$conf.int) # LR cutoff, OR, pval, low CI, high CI
  return(out)
}
df <- data.frame()
df <- rbind.data.frame(df, m_v_d(f))
df <- rbind.data.frame(df, m_v_d(f.chd))
df <- rbind.data.frame(df, m_v_d(f.syn))
df <- rbind.data.frame(df, m_v_d(f.syn.chd))
df <- rbind.data.frame(df, m_v_d(f.non))
df <- rbind.data.frame(df, m_v_d(f.non.chd))
df <- rbind.data.frame(df, m_v_d(f.unk))
df <- rbind.data.frame(df, m_v_d(f.unk.chd))
names(df) <- c('dn.ct', 'mo.ct', 'v.m.d', 'v.n.d', 'v.m.nd', 'v.n.nd', 'or', 'p', 'ci1', 'ci2')
src <- c('all', 'all_chd',  'syn', 'syn_chd',  'iso', 'iso_chd', 'unk', 'unk_chd')
df$src <- src
outtablename <- paste(op, 'enrichment.txt', sep='.')
write.table(df,outtablename, sep='\t', quote=F, row.names=F )

# For VAF bins 0.0-0.5 by0.1, are there bins for which mosaics contribute more to disease?
vaf_bin_analysis <- function(v){
  outdf <- NULL
  mo <- v[v$col=="red",]
  no <- v[v$col=="black",]
  c1 <- c(0.0, 0.125)
  c2 <- c(0.125, 0.4)
  for(i in 1:length(c1)){
    bin <- mo[mo$vaf>=c1[i] & mo$vaf<c2[i],]
    v.m.d <- bin[bin$vclass=="LGD" | bin$vclass=="Dmis",]
    v.m.nd <- bin[!(bin$vclass=="LGD" | bin$vclass=="Dmis"),]
    v.n.d <- no[no$vclass=="LGD" | no$vclass=="Dmis",]
    v.n.nd <- no[!(no$vclass=="LGD" | no$vclass=="Dmis"),]
    cont <- matrix(c(nrow(v.m.d), nrow(v.n.d), nrow(v.m.nd), nrow(v.n.nd)), nrow=2)
    fish <- fisher.test(cont)
    out <- c(c1[i], c2[i], nrow(no), nrow(mo), nrow(v.m.d), nrow(v.n.d), nrow(v.m.nd), nrow(v.n.nd), unname(fish$estimate), fish$p.value, fish$conf.int) # LR cutoff, OR, pval, low CI, high CI
    outdf <- rbind.data.frame(outdf, out)
  }
  names(outdf) <- c('c1', 'c2', 'germ.ct', 'mo.ct', 'v.m.d', 'v.n.d', 'v.m.nd', 'v.n.nd', 'or', 'p', 'ci1', 'ci2')
  return(outdf)
}
plot_vaf_bin <- function(df, name){
  title <- paste("Mosaic vs. Damaging", name, sep='\n')
  df$bincts <- df$v.m.d+df$v.m.nd
  labs <- paste(round(df$p, 2), paste(paste(df$v.m.d, df$v.m.nd, sep='/'), paste(df$v.n.d, df$v.n.nd, sep='/'), sep='\n'), sep='\n')
  fp <- ggplot(df, aes(df$c2, df$or)) +  
    geom_text(aes(label=labs), vjust=-0.05, size=8) +
    geom_bar(stat='identity') +
    geom_hline(yintercept=1, linetype="longdash") + 
    labs(x="VAF bin upper bound", y="Enrichment" ) +
    ylim(c(0, max(round(df$or, 0))+3)) +
    ggtitle(title) +
    theme(plot.title = element_text(size=20), axis.title = element_text(size=18), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) 
  return(fp)
}

## All samples
totdf <- NULL
labs <- c('all', 'all', 'all.chd', 'all.chd', 'syn', 'syn', 'syn.chd', 'syn.chd', 'iso', 'iso', 'iso.chd', 'iso.chd', 'unk', 'unk', 'unk.chd', 'unk.chd')
# all genes
vaf.df <- vaf_bin_analysis(f)
totdf <- rbind.data.frame(totdf, vaf.df)
# chd genes
vaf.chd.df <- vaf_bin_analysis(f.chd)
totdf <- rbind.data.frame(totdf, vaf.chd.df)
# plots
p.vaf.df <- plot_vaf_bin(vaf.df, 'All') + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
p.vaf.chd.df <- plot_vaf_bin(vaf.chd.df, 'chd')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotallname = paste(op, 'enrichment_vaf.ALL.pdf', sep='.')
pdf(outplotallname, width=14, height=7)
grid.arrange(p.vaf.df, p.vaf.chd.df, ncol=2)
dev.off()

## Syndromic Samples
# all genes
vaf.syn.df <- vaf_bin_analysis(f.syn)
totdf <- rbind.data.frame(totdf, vaf.syn.df)
# chd genes
vaf.syn.chd.df <- vaf_bin_analysis(f.syn.chd)
totdf <- rbind.data.frame(totdf, vaf.syn.chd.df)
# plots
p.vaf.syn.df <- plot_vaf_bin(vaf.syn.df, 'All') + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
p.vaf.syn.chd.df <- plot_vaf_bin(vaf.syn.chd.df, 'chd')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotsynname = paste(op, 'enrichment_vaf.SYN.pdf', sep='.')
pdf(outplotsynname, width=14, height=7)
grid.arrange(p.vaf.syn.df, p.vaf.syn.chd.df, ncol=2)
dev.off()

## Isolated Samples
# all genes
vaf.iso.df <- vaf_bin_analysis(f.non)
totdf <- rbind.data.frame(totdf, vaf.iso.df)
# chd genes
vaf.iso.chd.df <- vaf_bin_analysis(f.non.chd)
totdf <- rbind.data.frame(totdf, vaf.iso.chd.df)
# plots
p.vaf.iso.df <- plot_vaf_bin(vaf.iso.df, 'All') + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
p.vaf.iso.chd.df <- plot_vaf_bin(vaf.iso.chd.df, 'chd')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotisoname = paste(op, 'enrichment_vaf.ISO.pdf', sep='.')
pdf(outplotisoname, width=14, height=7)
grid.arrange(p.vaf.iso.df, p.vaf.iso.chd.df, ncol=2)
dev.off()

## Unknown Samples
# all genes
vaf.unk.df <- vaf_bin_analysis(f.unk)
totdf <- rbind.data.frame(totdf, vaf.unk.df)
# chd genes
vaf.unk.chd.df <- vaf_bin_analysis(f.unk.chd)
totdf <- rbind.data.frame(totdf, vaf.unk.chd.df)
# plots
p.vaf.unk.df <- plot_vaf_bin(vaf.unk.df, 'All') + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
p.vaf.unk.chd.df <- plot_vaf_bin(vaf.unk.chd.df, 'chd')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotunkname = paste(op, 'enrichment_vaf.UNK.pdf', sep='.')
pdf(outplotunkname, width=14, height=7)
grid.arrange(p.vaf.unk.df, p.vaf.unk.chd.df, ncol=2)
dev.off()

totdf$src <- labs
outVAFtablename <- paste(op, 'enrichment.VAFbins.txt')
write.table(totdf,outVAFtablename, sep='\t', quote=F, row.names=F )

print("Done !")
print(paste("output", c(outtablename, outVAFtablename, outplotallname, outplotsynname, outplotisoname), sep=' : '))

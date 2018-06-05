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

f <- read.table(fname, sep='\t', header=T, quote="")

f.syn <- f[f$NDD_diagnosis=="Yes" | f$extracardiac=="Yes",]
f.non <- f[(f$NDD_diagnosis=="No" & f$extracardiac=="No"),]
f.unk <- f[(f$NDD_diagnosis == "Unknown" & f$extracardiac=="No") | (f$NDD_diagnosis=="No" & f$extracardiac=="Unknown") | (f$NDD_diagnosis=="Unknown" & f$extracardiac=="Unknown"),]
print( c( paste('Total Variants', nrow(f), sep=' : '), paste('Syndromic', nrow(f.syn), sep=' : '), paste('Non-syndromic', nrow(f.non), sep=' : '), paste('NDD/CA Unknown', nrow(f.unk), sep=' : ') ) )
print( c( paste('Total Samples', length(unique(f$id)), sep=' : '), paste('Syndromic Samples', length(unique(f.syn$id)), sep=' : '), paste('Non-syndromic Samples', length(unique(f.non$id)), sep=' : '), paste('NDD/CA Unknown', length(unique(f.unk$id)), sep=' : ') ) )

summarize <- function(all, set){
  hhe <- all[!is.na(all$hexp),]
  all.hhe <- hhe[hhe$hexp>=75,]
  summarydf <- NULL
  cts <- c(nrow(all), nrow(all.hhe))
  summarydf <- rbind.data.frame(summarydf, cts)
  mocts <- c(nrow(all[all$col=="red",]), nrow(all.hhe[all.hhe$col=="red",]))
  summarydf <- rbind.data.frame(summarydf, mocts)
  germcts <- c(nrow(all[all$col=="black",]), nrow(all.hhe[all.hhe$col=="black",]))
  summarydf <- rbind.data.frame(summarydf, germcts)
  names(summarydf) <- c("All","HHE")
  summarydf$src <- c('Total', 'Mosaic', 'Germline')
  summarydf$set <- c(set, set, set)
  return(summarydf)
}
all.summary <- summarize(f, 'all')
syn.summary <- summarize(f.syn, 'syndromic')
non.summary <- summarize(f.non, 'non-syndromic')
unk.summary <- summarize(f.unk, 'unknown')

## Parse variants into different gene sets
# all samples
tmp.hhe <- f[!is.na(f$hexp),]
f.hhe <- tmp.hhe[tmp.hhe$hexp>=75,]
# syndromic
f.syn.hhe <- f.hhe[f.hhe$NDD_diagnosis=="Yes" | f.hhe$extracardiac=="Yes",]
# isolated
f.non.hhe <- f.hhe[!(f.hhe$NDD_diagnosis=="Yes" | f.hhe$extracardiac=="Yes"),]
# unknown
f.unk.hhe <- f.hhe[(f.hhe$NDD_diagnosis == "Unknown" & f.hhe$extracardiac=="No") | (f.hhe$NDD_diagnosis=="No" & f.hhe$extracardiac=="Unknown") | (f.hhe$NDD_diagnosis=="Unknown" & f.hhe$extracardiac=="Unknown"),]

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
df <- rbind.data.frame(df, m_v_d(f.hhe))
df <- rbind.data.frame(df, m_v_d(f.syn))
df <- rbind.data.frame(df, m_v_d(f.syn.hhe))
df <- rbind.data.frame(df, m_v_d(f.non))
df <- rbind.data.frame(df, m_v_d(f.non.hhe))
df <- rbind.data.frame(df, m_v_d(f.unk))
df <- rbind.data.frame(df, m_v_d(f.unk.hhe))
names(df) <- c('dn.ct', 'mo.ct', 'v.m.d', 'v.n.d', 'v.m.nd', 'v.n.nd', 'or', 'p', 'ci1', 'ci2')
src <- c('all', 'all_hhe',  'syn', 'syn_hhe',  'iso', 'iso_hhe', 'unk', 'unk_hhe')
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
labs <- c('all', 'all', 'all.hhe', 'all.hhe', 'syn', 'syn', 'syn.hhe', 'syn.hhe', 'iso', 'iso', 'iso.hhe', 'iso.hhe', 'unk', 'unk', 'unk.hhe', 'unk.hhe')
# all genes
vaf.df <- vaf_bin_analysis(f)
totdf <- rbind.data.frame(totdf, vaf.df)
# HHE genes
vaf.hhe.df <- vaf_bin_analysis(f.hhe)
totdf <- rbind.data.frame(totdf, vaf.hhe.df)
# plots
p.vaf.df <- plot_vaf_bin(vaf.df, 'All') + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
p.vaf.hhe.df <- plot_vaf_bin(vaf.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotallname = paste(op, 'enrichment_vaf.ALL.pdf', sep='.')
pdf(outplotallname, width=14, height=7)
grid.arrange(p.vaf.df, p.vaf.hhe.df, ncol=2)
dev.off()

## Syndromic Samples
# all genes
vaf.syn.df <- vaf_bin_analysis(f.syn)
totdf <- rbind.data.frame(totdf, vaf.syn.df)
# HHE genes
vaf.syn.hhe.df <- vaf_bin_analysis(f.syn.hhe)
totdf <- rbind.data.frame(totdf, vaf.syn.hhe.df)
# plots
p.vaf.syn.df <- plot_vaf_bin(vaf.syn.df, 'All') + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
p.vaf.syn.hhe.df <- plot_vaf_bin(vaf.syn.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotsynname = paste(op, 'enrichment_vaf.SYN.pdf', sep='.')
pdf(outplotsynname, width=14, height=7)
grid.arrange(p.vaf.syn.df, p.vaf.syn.hhe.df, ncol=2)
dev.off()

## Isolated Samples
# all genes
vaf.iso.df <- vaf_bin_analysis(f.non)
totdf <- rbind.data.frame(totdf, vaf.iso.df)
# HHE genes
vaf.iso.hhe.df <- vaf_bin_analysis(f.non.hhe)
totdf <- rbind.data.frame(totdf, vaf.iso.hhe.df)
# plots
p.vaf.iso.df <- plot_vaf_bin(vaf.iso.df, 'All') + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
p.vaf.iso.hhe.df <- plot_vaf_bin(vaf.iso.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotisoname = paste(op, 'enrichment_vaf.ISO.pdf', sep='.')
pdf(outplotisoname, width=14, height=7)
grid.arrange(p.vaf.iso.df, p.vaf.iso.hhe.df, ncol=2)
dev.off()

## Unknown Samples
# all genes
vaf.unk.df <- vaf_bin_analysis(f.unk)
totdf <- rbind.data.frame(totdf, vaf.unk.df)
# HHE genes
vaf.unk.hhe.df <- vaf_bin_analysis(f.unk.hhe)
totdf <- rbind.data.frame(totdf, vaf.unk.hhe.df)
# plots
p.vaf.unk.df <- plot_vaf_bin(vaf.unk.df, 'All') + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
p.vaf.unk.hhe.df <- plot_vaf_bin(vaf.unk.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotunkname = paste(op, 'enrichment_vaf.UNK.pdf', sep='.')
pdf(outplotunkname, width=14, height=7)
grid.arrange(p.vaf.unk.df, p.vaf.unk.hhe.df, ncol=2)
dev.off()

totdf$src <- labs
outVAFtablename <- paste(op, 'enrichment.VAFbins.txt')
write.table(totdf,outVAFtablename, sep='\t', quote=F, row.names=F )

print("Done !")
print(paste("output", c(outtablename, outVAFtablename, outplotallname, outplotsynname, outplotisoname), sep=' : '))

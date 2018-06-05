## Usage: Rscript pcgc_case_only_analysis.HHE_only.R <annotated denovos file> <output prefix>
## Purpose: Run Case-only analysis comparing LGD+Dmis(+Bmis) vs. Synonymous mosaic and germline sites to calculate enrichment and to plot trend
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

## Parse variants into different gene sets
# all samples
tmp.hhe <- f[!is.na(f$hexp),]
f.hhe <- tmp.hhe[tmp.hhe$hexp>=75,]
f.nhhe <- tmp.hhe[!(tmp.hhe$hexp>=75),]
# syndromic
f.syn.hhe <- f.hhe[f.hhe$NDD_diagnosis=="Yes" | f.hhe$extracardiac=="Yes",]
f.syn.nhhe <- f.nhhe[f.nhhe$NDD_diagnosis=="Yes" | f.nhhe$extracardiac=="Yes",]
# isolated
f.non.hhe <- f.hhe[(f.hhe$NDD_diagnosis=="No" & f.hhe$extracardiac=="No"),]
f.non.nhhe <- f.nhhe[(f.nhhe$NDD_diagnosis=="No" & f.nhhe$extracardiac=="No"),]
# unknown
f.unk.hhe <- f.hhe[(f.hhe$NDD_diagnosis == "Unknown" & f.hhe$extracardiac=="No") | (f.hhe$NDD_diagnosis=="No" & f.hhe$extracardiac=="Unknown") | (f.hhe$NDD_diagnosis=="Unknown" & f.hhe$extracardiac=="Unknown"),]
f.unk.nhhe <- f.nhhe[(f.nhhe$NDD_diagnosis == "Unknown" & f.nhhe$extracardiac=="No") | (f.nhhe$NDD_diagnosis=="No" & f.nhhe$extracardiac=="Unknown") | (f.nhhe$NDD_diagnosis=="Unknown" & f.nhhe$extracardiac=="Unknown"),]

## Print summary counts
# All Samples
mo.summarydf <- NULL
germ.summarydf <- NULL
# mosaic
all.mo.summarydf <- NULL
all.mo.summarydf <- rbind.data.frame(all.mo.summarydf, unname(table(f.hhe[f.hhe$col=="red",]$vclass)))
all.mo.summarydf <- rbind.data.frame(all.mo.summarydf, unname(table(f.nhhe[f.nhhe$col=="red",]$vclass)))
colnames(all.mo.summarydf) <- names(table(f$vclass))
rownames(all.mo.summarydf) <- c('HHE', 'nHHE')
print('## ALL: mosaic summary')
print(all.mo.summarydf)
# germline
all.germ.summarydf <- NULL
all.germ.summarydf <- rbind.data.frame(all.germ.summarydf, unname(table(f.hhe[f.hhe$col=="black",]$vclass)))
all.germ.summarydf <- rbind.data.frame(all.germ.summarydf, unname(table(f.nhhe[f.nhhe$col=="black",]$vclass)))
colnames(all.germ.summarydf) <- names(table(f$vclass))
rownames(all.germ.summarydf) <- c('HHE', 'nHHE')
print('## ALL: germline summary')
print(all.germ.summarydf)

mo.summarydf <- rbind.data.frame(mo.summarydf, all.mo.summarydf)
germ.summarydf <- rbind.data.frame(germ.summarydf, all.germ.summarydf)

# Syndromic Samples
# mosaic
syn.mo.summarydf <- NULL
syn.mo.summarydf <- rbind.data.frame(syn.mo.summarydf, unname(table(f.syn.hhe[f.syn.hhe$col=="red",]$vclass)))
syn.mo.summarydf <- rbind.data.frame(syn.mo.summarydf, unname(table(f.syn.nhhe[f.syn.nhhe$col=="red",]$vclass)))
colnames(syn.mo.summarydf) <- names(table(f$vclass))
rownames(syn.mo.summarydf) <- c('HHE', 'nHHE')
print('## SYNDROMIC: mosaic summary')
print(syn.mo.summarydf)
# germline
syn.germ.summarydf <- NULL
syn.germ.summarydf <- rbind.data.frame(syn.germ.summarydf, unname(table(f.syn.hhe[f.syn.hhe$col=="black",]$vclass)))
syn.germ.summarydf <- rbind.data.frame(syn.germ.summarydf, unname(table(f.syn.nhhe[f.syn.nhhe$col=="black",]$vclass)))
colnames(syn.germ.summarydf) <- names(table(f$vclass))
rownames(syn.germ.summarydf) <- c('HHE', 'nHHE')
print('## SYNDROMIC: germline summary')
print(syn.germ.summarydf)

mo.summarydf <- rbind.data.frame(mo.summarydf, syn.mo.summarydf)
germ.summarydf <- rbind.data.frame(germ.summarydf, syn.germ.summarydf)

# Isolated Samples
# mosaic
non.mo.summarydf <- NULL
non.mo.summarydf <- rbind.data.frame(non.mo.summarydf, unname(table(f.non.hhe[f.non.hhe$col=="red",]$vclass)))
non.mo.summarydf <- rbind.data.frame(non.mo.summarydf, unname(table(f.non.nhhe[f.non.nhhe$col=="red",]$vclass)))
colnames(non.mo.summarydf) <- names(table(f$vclass))
rownames(non.mo.summarydf) <- c('HHE', 'nHHE')
print('## ISOLATED: mosaic summary')
print(non.mo.summarydf)
# germline
non.germ.summarydf <- NULL
non.germ.summarydf <- rbind.data.frame(non.germ.summarydf, unname(table(f.non.hhe[f.non.hhe$col=="black",]$vclass)))
non.germ.summarydf <- rbind.data.frame(non.germ.summarydf, unname(table(f.non.nhhe[f.non.nhhe$col=="black",]$vclass)))
colnames(non.germ.summarydf) <- names(table(f$vclass))
rownames(non.germ.summarydf) <- c('HHE', 'nHHE')
print('## ISOLATED: germline summary')
print(non.germ.summarydf)

mo.summarydf <- rbind.data.frame(mo.summarydf, non.mo.summarydf)
germ.summarydf <- rbind.data.frame(germ.summarydf, non.germ.summarydf)

# Unknown Samples
# mosaic
unk.mo.summarydf <- NULL
unk.mo.summarydf <- rbind.data.frame(unk.mo.summarydf, unname(table(f.unk.hhe[f.unk.hhe$col=="red",]$vclass)))
unk.mo.summarydf <- rbind.data.frame(unk.mo.summarydf, unname(table(f.unk.nhhe[f.unk.nhhe$col=="red",]$vclass)))
colnames(unk.mo.summarydf) <- names(table(f$vclass))
rownames(unk.mo.summarydf) <- c('HHE', 'nHHE')
print('## UNKNOWN: mosaic summary')
print(unk.mo.summarydf)
# germline
unk.germ.summarydf <- NULL
unk.germ.summarydf <- rbind.data.frame(unk.germ.summarydf, unname(table(f.unk.hhe[f.unk.hhe$col=="black",]$vclass)))
unk.germ.summarydf <- rbind.data.frame(unk.germ.summarydf, unname(table(f.unk.nhhe[f.unk.nhhe$col=="black",]$vclass)))
colnames(unk.germ.summarydf) <- names(table(f$vclass))
rownames(unk.germ.summarydf) <- c('HHE', 'nHHE')
print('## UNKNOWN: germline summary')
print(unk.germ.summarydf)

mo.summarydf <- rbind.data.frame(mo.summarydf, unk.mo.summarydf)
germ.summarydf <- rbind.data.frame(germ.summarydf, unk.germ.summarydf)

## write out summary tables
labs <- c('all_HHE', 'all_nHHE', 'syn_HHE', 'syn_nHHE', 'iso_HHE', 'iso_nHHE', 'unk_HHE', 'unk_nHHE')
mo.summarydf$src <- labs
germ.summarydf$src <- labs
write.table(mo.summarydf, paste(op, 'mosaic_summary.txt', sep='.'), sep='\t', quote=F, row.names=F)
write.table(germ.summarydf, paste(op, 'germline_summary.txt', sep='.'), sep='\t', quote=F, row.names=F)

## For specified sets of genes, is there an enrichment of LGD+Dmis mutations vs. Synonymous, compared to other genes?
naive_enrichment <- function(v1, v2){
  mo1 <- v1[v1$col=="red",]
  mo2 <- v2[v2$col=="red",]
  mo1.dam <- mo1[(mo1$vclass=="LGD" | mo1$vclass=="Dmis"),]
  mo1.syn <- mo1[mo1$vclass=="synonymous",]
  mo2.dam <- mo2[(mo2$vclass=="LGD" | mo2$vclass=="Dmis"),]
  mo2.syn <- mo2[mo2$vclass=="synonymous",]
  cont <- matrix(c(nrow(mo1.dam),  nrow(mo2.dam), nrow(mo1.syn),nrow(mo2.syn)), nrow=2)
  fish <- fisher.test(cont)
  out <- c(nrow(mo1), nrow(mo2), nrow(mo1.dam), nrow(mo2.dam), nrow(mo1.syn), nrow(mo2.syn), unname(fish$estimate), fish$p.value, fish$conf.int) # LR cutoff, OR, pval, low CI, high CI
  return(out)  
}
mo.df <- NULL
# all samples
mo.df <- rbind.data.frame(mo.df, naive_enrichment(f.hhe, f.nhhe))
# syndromic 
mo.df <- rbind.data.frame(mo.df, naive_enrichment(f.syn.hhe, f.syn.nhhe))
# isolated
mo.df <- rbind.data.frame(mo.df, naive_enrichment(f.non.hhe, f.non.nhhe))
# unknown
mo.df <- rbind.data.frame(mo.df, naive_enrichment(f.unk.hhe, f.unk.nhhe))
names(mo.df) <- c('in.ct', 'out.ct', 'in.LGD_Dmis', 'out.LGD_Dmis', 'in.Syn', 'out.Syn', 'or', 'p', 'ci1', 'ci2')
src <- c('all_hhe', 'syn_hhe', 'iso_hhe', 'unk_hhe')
mo.df$src <- src
mo_lgd_dmis_name <- paste(op, 'LGD_Dmis_vs_Syn.MOSAIC.txt', sep='.')
write.table(mo.df, mo_lgd_dmis_name, sep='\t', quote=F, row.names=F)

naive_enrichment_germ <- function(v1, v2){
  mo1 <- v1[v1$col=="black",]
  mo2 <- v2[v2$col=="black",]
  mo1.dam <- mo1[(mo1$vclass=="LGD" | mo1$vclass=="Dmis"),]
  mo1.syn <- mo1[mo1$vclass=="synonymous",]
  mo2.dam <- mo2[(mo2$vclass=="LGD" | mo2$vclass=="Dmis"),]
  mo2.syn <- mo2[mo2$vclass=="synonymous",]
  cont <- matrix(c(nrow(mo1.dam),  nrow(mo2.dam), nrow(mo1.syn),nrow(mo2.syn)), nrow=2)
  fish <- fisher.test(cont)
  out <- c(nrow(mo1), nrow(mo2), nrow(mo1.dam), nrow(mo2.dam), nrow(mo1.syn), nrow(mo2.syn), unname(fish$estimate), fish$p.value, fish$conf.int) # LR cutoff, OR, pval, low CI, high CI
  return(out)  
}
germ.df <- NULL
# all samples
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ(f.hhe, f.nhhe))
# syndromic 
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ(f.syn.hhe, f.syn.nhhe))
# isolated
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ(f.non.hhe, f.non.nhhe))
# unknown
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ(f.unk.hhe, f.unk.nhhe))
names(germ.df) <- c('in.ct', 'out.ct', 'in.LGD_Dmis', 'out.LGD_Dmis', 'in.Syn', 'out.Syn', 'or', 'p', 'ci1', 'ci2')
src <- c('all_hhe', 'syn_hhe', 'iso_hhe', 'unk_hhe')
germ.df$src <- src
germ_lgd_dmis_name <- paste(op, 'LGD_Dmis_vs_Syn.GERMLINE.txt', sep='.')
write.table(germ.df, germ_lgd_dmis_name, sep='\t', quote=F, row.names=F)


vaf_bin_analysis <- function(v1, v2){
  outdf <- NULL
  mo1 <- v1[v1$col=="red",]
  mo2 <- v2[v2$col=="red",]
  c1 <- c(0.0, 0.125)
  c2 <- c(0.125, 0.4)
  #c1 <- c(0.0, 0.0625, 0.125, 0.25)
  #c2 <- c(0.0625, 0.125, 0.25, 0.4)
  for(i in 1:length(c1)){
    inbin <- mo1[mo1$vaf>=c1[i] & mo1$vaf<c2[i],]
    outbin <- mo2[mo2$vaf>=c1[i] & mo2$vaf<c2[i],]
    v.in.d <- inbin[(inbin$vclass=="LGD" | inbin$vclass=="Dmis"),]
    v.in.nd <- inbin[inbin$vclass=="synonymous",]
    v.out.d <- outbin[(outbin$vclass=="LGD" | outbin$vclass=="Dmis"),]
    v.out.nd <- outbin[outbin$vclass=="synonymous",]
    cont <- matrix(c(nrow(v.in.d), nrow(v.out.d), nrow(v.in.nd), nrow(v.out.nd)), nrow=2)
    fish <- fisher.test(cont)
    out <- c(c1[i], c2[i], nrow(inbin), nrow(outbin), nrow(v.in.d), nrow(v.out.d), nrow(v.in.nd), nrow(v.out.nd), unname(fish$estimate), fish$p.value, fish$conf.int) # LR cutoff, OR, pval, low CI, high CI
    outdf <- rbind.data.frame(outdf, out)
  }
  names(outdf) <- c('c1', 'c2', 'in.ct', 'out.ct', 'v.in.d', 'v.out.d', 'v.in.nd', 'v.out.nd', 'or', 'p', 'ci1', 'ci2')
  return(outdf)
}
plot_vaf_bin <- function(df, name){
  title <- paste("LGD+Dmis vs. Syn", name, sep='\n')
  df$bincts <- df$v.in.d+df$v.in.nd
  df$bin_mean <- (df$c1+df$c2)/2
  labs <- paste(round(df$p, 2), paste(paste(df$v.in.d, df$v.in.nd, sep='/'), paste(df$v.out.d, df$v.out.nd, sep='/'), sep='\n'), sep='\n')
  fp <- ggplot(df, aes(df$c2, df$or)) +  
    geom_text(aes(label=labs), vjust=-0.25, size=8) +
    geom_bar(stat='identity') +
    geom_hline(yintercept=1, linetype="longdash") + 
    labs(x="VAF bin upper bound", y="Enrichment" ) +
    ylim(c(0, min(20, max(round(df$or, 0))+3))) +
    ggtitle(title) +
    theme(plot.title = element_text(size=20), axis.title = element_text(size=18), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) 
  return(fp)
}

## All samples
totdf <- NULL
src <- c('all.hhe', 'all.hhe', 'syn.hhe', 'syn.hhe', 'iso.hhe', 'iso.hhe', 'unk.hhe', 'unk.hhe')
# HHE genes
vaf.hhe.df <- vaf_bin_analysis(f.hhe, f.nhhe)
totdf <- rbind.data.frame(totdf, vaf.hhe.df)
# plots
p.vaf.hhe.df <- plot_vaf_bin(vaf.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotall_dmis_name = paste(op, 'LGD_Dmis_vs_Syn_vaf.ALL.pdf', sep='.')
pdf(outplotall_dmis_name, width=7, height=7)
p.vaf.hhe.df
dev.off()

## Syndromic Samples
# HHE genes
vaf.syn.hhe.df <- vaf_bin_analysis(f.syn.hhe, f.syn.nhhe)
totdf <- rbind.data.frame(totdf, vaf.syn.hhe.df)
# plots
p.vaf.syn.hhe.df <- plot_vaf_bin(vaf.syn.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotsyn_dmis_name = paste(op, 'LGD_Dmis_vs_Syn_vaf.SYN.pdf', sep='.')
pdf(outplotsyn_dmis_name, width=7, height=7)
p.vaf.syn.hhe.df
dev.off()

## Isolated Samples
# HHE genes
vaf.iso.hhe.df <- vaf_bin_analysis(f.non.hhe, f.non.nhhe)
totdf <- rbind.data.frame(totdf, vaf.iso.hhe.df)
# plots
p.vaf.iso.hhe.df <- plot_vaf_bin(vaf.iso.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotiso_dmis_name = paste(op, 'LGD_Dmis_vs_Syn_vaf.ISO.pdf', sep='.')
pdf(outplotiso_dmis_name, width=7, height=7)
p.vaf.iso.hhe.df
dev.off()

## Unknown Samples
# HHE genes
vaf.unk.hhe.df <- vaf_bin_analysis(f.unk.hhe, f.unk.nhhe)
totdf <- rbind.data.frame(totdf, vaf.unk.hhe.df)
# plots
p.vaf.unk.hhe.df <- plot_vaf_bin(vaf.unk.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotunk_dmis_name = paste(op, 'LGD_Dmis_vs_Syn_vaf.UNK.pdf', sep='.')
pdf(outplotunk_dmis_name, width=7, height=7)
p.vaf.unk.hhe.df
dev.off()

totdf$src <- src
lgd_dmis_vaf_name <- paste(op, 'LGD_Dmis_vs_Syn_vaf.txt', sep='.')
write.table(totdf, lgd_dmis_vaf_name, sep='\t', quote=F, row.names=F)


## LGD + Dmis + Bmis vs. Syn
naive_enrichment_mis <- function(v1, v2){
  mo1 <- v1[v1$col=="red",]
  mo2 <- v2[v2$col=="red",]
  mo1.dam <- mo1[(mo1$vclass=="LGD" | mo1$vclass=="Dmis" | mo1$vclass=="Bmis"),]
  mo1.syn <- mo1[mo1$vclass=="synonymous",]
  mo2.dam <- mo2[(mo2$vclass=="LGD" | mo2$vclass=="Dmis"| mo2$vclass=="Bmis"),]
  mo2.syn <- mo2[mo2$vclass=="synonymous",]
  cont <- matrix(c(nrow(mo1.dam),  nrow(mo2.dam), nrow(mo1.syn),nrow(mo2.syn)), nrow=2)
  fish <- fisher.test(cont)
  out <- c(nrow(mo1), nrow(mo2), nrow(mo1.dam), nrow(mo2.dam), nrow(mo1.syn), nrow(mo2.syn), unname(fish$estimate), fish$p.value, fish$conf.int) # LR cutoff, OR, pval, low CI, high CI
  return(out)  
}
mo.df <- NULL
# all samples
mo.df <- rbind.data.frame(mo.df, naive_enrichment_mis(f.hhe, f.nhhe))
# syndromic 
mo.df <- rbind.data.frame(mo.df, naive_enrichment_mis(f.syn.hhe, f.syn.nhhe))
# isolated
mo.df <- rbind.data.frame(mo.df, naive_enrichment_mis(f.non.hhe, f.non.nhhe))
# unknown
mo.df <- rbind.data.frame(mo.df, naive_enrichment_mis(f.unk.hhe, f.unk.nhhe))
names(mo.df) <- c('in.ct', 'out.ct', 'in.LGD_mis', 'out.LGD_mis', 'in.Syn', 'out.Syn', 'or', 'p', 'ci1', 'ci2')
src <- c('all_hhe', 'syn_hhe', 'iso_hhe', 'unk_hhe')
mo.df$src <- src
mo_lgd_mis_name <- paste(op, 'LGD_MIS_vs_Syn.MOSAIC.txt', sep='.')
write.table(mo.df, mo_lgd_mis_name, sep='\t', quote=F, row.names=F)

naive_enrichment_germ_mis <- function(v1, v2){
  mo1 <- v1[v1$col=="black",]
  mo2 <- v2[v2$col=="black",]
  mo1.dam <- mo1[(mo1$vclass=="LGD" | mo1$vclass=="Dmis" | mo1$vclass=="Bmis"),]
  mo1.syn <- mo1[mo1$vclass=="synonymous",]
  mo2.dam <- mo2[(mo2$vclass=="LGD" | mo2$vclass=="Dmis" | mo2$vclass=="Bmis"),]
  mo2.syn <- mo2[mo2$vclass=="synonymous",]
  cont <- matrix(c(nrow(mo1.dam),  nrow(mo2.dam), nrow(mo1.syn),nrow(mo2.syn)), nrow=2)
  fish <- fisher.test(cont)
  out <- c(nrow(mo1), nrow(mo2), nrow(mo1.dam), nrow(mo2.dam), nrow(mo1.syn), nrow(mo2.syn), unname(fish$estimate), fish$p.value, fish$conf.int) # LR cutoff, OR, pval, low CI, high CI
  return(out)  
}
germ.df <- NULL
# all samples
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ_mis(f.hhe, f.nhhe))
# syndromic 
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ_mis(f.syn.hhe, f.syn.nhhe))
# isolated
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ_mis(f.non.hhe, f.non.nhhe))
# isolated
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ_mis(f.unk.hhe, f.unk.nhhe))
names(germ.df) <- c('in.ct', 'out.ct', 'in.LGD_mis', 'out.LGD_mis', 'in.Syn', 'out.Syn', 'or', 'p', 'ci1', 'ci2')
src <- c('all_hhe', 'syn_hhe', 'iso_hhe', 'unk_hhe')
germ.df$src <- src
#germ.df
germ_lgd_mis_name <- paste(op, 'LGD_MIS_vs_Syn.GERMLINE.txt', sep='.')
write.table(germ.df, germ_lgd_mis_name, sep='\t', quote=F, row.names=F)


vaf_bin_analysis_mis <- function(v1, v2){
  outdf <- NULL
  mo1 <- v1[v1$col=="red",]
  mo2 <- v2[v2$col=="red",]
  c1 <- c(0.0, 0.125)
  c2 <- c(0.125, 0.4)
  #c1 <- c(0.0, 0.0625, 0.125, 0.25)
  #c2 <- c(0.0625, 0.125, 0.25, 0.4)
  for(i in 1:length(c1)){
    inbin <- mo1[mo1$vaf>=c1[i] & mo1$vaf<c2[i],]
    outbin <- mo2[mo2$vaf>=c1[i] & mo2$vaf<c2[i],]
    v.in.d <- inbin[(inbin$vclass=="LGD" | inbin$vclass=="Dmis" | inbin$vclass=="Bmis"),]
    v.in.nd <- inbin[inbin$vclass=="synonymous",]
    v.out.d <- outbin[(outbin$vclass=="LGD" | outbin$vclass=="Dmis" | outbin$vclass=="Bmis"),]
    v.out.nd <- outbin[outbin$vclass=="synonymous",]
    cont <- matrix(c(nrow(v.in.d), nrow(v.out.d), nrow(v.in.nd), nrow(v.out.nd)), nrow=2)
    fish <- fisher.test(cont)
    out <- c(c1[i], c2[i], nrow(inbin), nrow(outbin), nrow(v.in.d), nrow(v.out.d), nrow(v.in.nd), nrow(v.out.nd), unname(fish$estimate), fish$p.value, fish$conf.int) # LR cutoff, OR, pval, low CI, high CI
    outdf <- rbind.data.frame(outdf, out)
  }
  names(outdf) <- c('c1', 'c2', 'in.ct', 'out.ct', 'v.in.d', 'v.out.d', 'v.in.nd', 'v.out.nd', 'or', 'p', 'ci1', 'ci2')
  return(outdf)
}
plot_vaf_bin_mis <- function(df, name){
  title <- paste("LGD+MIS vs. Syn", name, sep='\n')
  df$bincts <- df$v.in.d+df$v.in.nd
  df$bin_mean <- (df$c1+df$c2)/2
  labs <- paste(round(df$p, 2), paste(paste(df$v.in.d, df$v.in.nd, sep='/'), paste(df$v.out.d, df$v.out.nd, sep='/'), sep='\n'), sep='\n')
  fp <- ggplot(df, aes(df$c2, df$or)) +  
    geom_text(aes(label=labs), vjust=-0.25, size=8) +
    geom_bar(stat='identity') +
    geom_hline(yintercept=1, linetype="longdash") + 
    labs(x="VAF bin upper bound", y="Enrichment" ) +
    ylim(c(0, min(20, max(round(df$or, 0))+3))) +
    ggtitle(title) +
    theme(plot.title = element_text(size=20), axis.title = element_text(size=18), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) 
  return(fp)
}

## All samples
totdf <- NULL
src <- c('all.hhe', 'all.hhe', 'syn.hhe', 'syn.hhe', 'iso.hhe', 'iso.hhe', 'unk.hhe', 'unk.hhe')
# HHE genes
vaf.hhe.df <- vaf_bin_analysis_mis(f.hhe, f.nhhe)
totdf <- rbind.data.frame(totdf, vaf.hhe.df)
# plots
p.vaf.hhe.df <- plot_vaf_bin_mis(vaf.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotall_mis_name = paste(op, 'LGD_MIS_vs_Syn_vaf.ALL.pdf', sep='.')
pdf(outplotall_mis_name, width=7, height=7)
p.vaf.hhe.df
dev.off()

## Syndromic Samples
# HHE genes
vaf.syn.hhe.df <- vaf_bin_analysis_mis(f.syn.hhe, f.syn.nhhe)
totdf <- rbind.data.frame(totdf, vaf.syn.hhe.df)
# plots
p.vaf.syn.hhe.df <- plot_vaf_bin_mis(vaf.syn.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotsyn_mis_name = paste(op, 'LGD_MIS_vs_Syn_vaf.SYN.pdf', sep='.')
pdf(outplotsyn_mis_name, width=7, height=7)
p.vaf.syn.hhe.df
dev.off()

## Isolated Samples
# HHE genes
vaf.iso.hhe.df <- vaf_bin_analysis_mis(f.non.hhe, f.non.nhhe)
totdf <- rbind.data.frame(totdf, vaf.iso.hhe.df)
# plots
p.vaf.iso.hhe.df <- plot_vaf_bin_mis(vaf.iso.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotiso_mis_name = paste(op, 'LGD_MIS_vs_Syn_vaf.ISO.pdf', sep='.')
pdf(outplotiso_mis_name, width=7, height=7)
p.vaf.iso.hhe.df
dev.off()

## Unknown Samples
# HHE genes
vaf.unk.hhe.df <- vaf_bin_analysis_mis(f.unk.hhe, f.unk.nhhe)
totdf <- rbind.data.frame(totdf, vaf.unk.hhe.df)
# plots
p.vaf.unk.hhe.df <- plot_vaf_bin_mis(vaf.unk.hhe.df, 'HHE')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotunk_mis_name = paste(op, 'LGD_MIS_vs_Syn_vaf.UNK.pdf', sep='.')
pdf(outplotunk_mis_name, width=7, height=7)
p.vaf.unk.hhe.df
dev.off()

totdf$src <- src
lgd_mis_vaf_name <- paste(op, 'LGD_MIS_vs_Syn_vaf.txt', sep='.')
write.table(totdf, lgd_mis_vaf_name, sep='\t', quote=F, row.names=F)

print("Done !")
print(paste("output", c(paste(op, 'mosaic_summary.txt', sep='.'), paste(op, 'germline_summary.txt', sep='.'), lgd_dmis_vaf_name, outplotall_dmis_name, outplotsyn_dmis_name, outplotiso_dmis_name, outplotunk_dmis_name, lgd_mis_vaf_name, outplotall_mis_name, outplotsyn_mis_name, outplotiso_mis_name, outplotunk_mis_name), sep=' : '))

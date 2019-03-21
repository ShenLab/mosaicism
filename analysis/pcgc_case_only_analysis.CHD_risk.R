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
if("VAF" %in% names(f)){
  f$vaf <- f$VAF
}
f.syn <- f[f$NDD_diagnosis=="Yes" | f$extracardiac=="Yes",]
f.non <- f[(f$NDD_diagnosis=="No" & f$extracardiac=="No"),]
f.unk <- f[(f$NDD_diagnosis == "Unknown" & f$extracardiac=="No") | (f$NDD_diagnosis=="No" & f$extracardiac=="Unknown") | (f$NDD_diagnosis=="Unknown" & f$extracardiac=="Unknown"),]
print( c( paste('Total Variants', nrow(f), sep=' : '), paste('Syndromic', nrow(f.syn), sep=' : '), paste('Non-syndromic', nrow(f.non), sep=' : '), paste('NDD/CA Unknown', nrow(f.unk), sep=' : ') ) )
print( c( paste('Total Samples', length(unique(f$id)), sep=' : '), paste('Syndromic Samples', length(unique(f.syn$id)), sep=' : '), paste('Non-syndromic Samples', length(unique(f.non$id)), sep=' : '), paste('NDD/CA Unknown', length(unique(f.unk$id)), sep=' : ') ) )

## Parse variants into different gene sets
# all samples
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

## Print summary counts
# All Samples
mo.summarydf <- NULL
germ.summarydf <- NULL
# mosaic
all.mo.summarydf <- NULL
all.mo.summarydf <- rbind.data.frame(all.mo.summarydf, unname(table(f.chd[f.chd$col=="red",]$vclass)))
all.mo.summarydf <- rbind.data.frame(all.mo.summarydf, unname(table(f.nchd[f.nchd$col=="red",]$vclass)))
colnames(all.mo.summarydf) <- names(table(f$vclass))
rownames(all.mo.summarydf) <- c('CHD', 'nCHD')
print('## ALL: mosaic summary')
print(all.mo.summarydf)
# germline
all.germ.summarydf <- NULL
all.germ.summarydf <- rbind.data.frame(all.germ.summarydf, unname(table(f.chd[f.chd$col=="black",]$vclass)))
all.germ.summarydf <- rbind.data.frame(all.germ.summarydf, unname(table(f.nchd[f.nchd$col=="black",]$vclass)))
colnames(all.germ.summarydf) <- names(table(f$vclass))
rownames(all.germ.summarydf) <- c('CHD', 'nCHD')
print('## ALL: germline summary')
print(all.germ.summarydf)

mo.summarydf <- rbind.data.frame(mo.summarydf, all.mo.summarydf)
germ.summarydf <- rbind.data.frame(germ.summarydf, all.germ.summarydf)

# Syndromic Samples
# mosaic
syn.mo.summarydf <- NULL
syn.mo.summarydf <- rbind.data.frame(syn.mo.summarydf, unname(table(f.syn.chd[f.syn.chd$col=="red",]$vclass)))
syn.mo.summarydf <- rbind.data.frame(syn.mo.summarydf, unname(table(f.syn.nchd[f.syn.nchd$col=="red",]$vclass)))
colnames(syn.mo.summarydf) <- names(table(f$vclass))
rownames(syn.mo.summarydf) <- c('chd', 'nchd')
print('## SYNDROMIC: mosaic summary')
print(syn.mo.summarydf)
# germline
syn.germ.summarydf <- NULL
syn.germ.summarydf <- rbind.data.frame(syn.germ.summarydf, unname(table(f.syn.chd[f.syn.chd$col=="black",]$vclass)))
syn.germ.summarydf <- rbind.data.frame(syn.germ.summarydf, unname(table(f.syn.nchd[f.syn.nchd$col=="black",]$vclass)))
colnames(syn.germ.summarydf) <- names(table(f$vclass))
rownames(syn.germ.summarydf) <- c('chd', 'nchd')
print('## SYNDROMIC: germline summary')
print(syn.germ.summarydf)

mo.summarydf <- rbind.data.frame(mo.summarydf, syn.mo.summarydf)
germ.summarydf <- rbind.data.frame(germ.summarydf, syn.germ.summarydf)

# Isolated Samples
# mosaic
non.mo.summarydf <- NULL
non.mo.summarydf <- rbind.data.frame(non.mo.summarydf, unname(table(f.non.chd[f.non.chd$col=="red",]$vclass)))
non.mo.summarydf <- rbind.data.frame(non.mo.summarydf, unname(table(f.non.nchd[f.non.nchd$col=="red",]$vclass)))
colnames(non.mo.summarydf) <- names(table(f$vclass))
rownames(non.mo.summarydf) <- c('chd', 'nchd')
print('## ISOLATED: mosaic summary')
print(non.mo.summarydf)
# germline
non.germ.summarydf <- NULL
non.germ.summarydf <- rbind.data.frame(non.germ.summarydf, unname(table(f.non.chd[f.non.chd$col=="black",]$vclass)))
non.germ.summarydf <- rbind.data.frame(non.germ.summarydf, unname(table(f.non.nchd[f.non.nchd$col=="black",]$vclass)))
colnames(non.germ.summarydf) <- names(table(f$vclass))
rownames(non.germ.summarydf) <- c('chd', 'nchd')
print('## ISOLATED: germline summary')
print(non.germ.summarydf)

mo.summarydf <- rbind.data.frame(mo.summarydf, non.mo.summarydf)
germ.summarydf <- rbind.data.frame(germ.summarydf, non.germ.summarydf)

# Unknown Samples
# mosaic
unk.mo.summarydf <- NULL
unk.mo.summarydf <- rbind.data.frame(unk.mo.summarydf, unname(table(f.unk.chd[f.unk.chd$col=="red",]$vclass)))
unk.mo.summarydf <- rbind.data.frame(unk.mo.summarydf, unname(table(f.unk.nchd[f.unk.nchd$col=="red",]$vclass)))
colnames(unk.mo.summarydf) <- names(table(f$vclass))
rownames(unk.mo.summarydf) <- c('chd', 'nchd')
print('## UNKNOWN: mosaic summary')
print(unk.mo.summarydf)
# germline
unk.germ.summarydf <- NULL
unk.germ.summarydf <- rbind.data.frame(unk.germ.summarydf, unname(table(f.unk.chd[f.unk.chd$col=="black",]$vclass)))
unk.germ.summarydf <- rbind.data.frame(unk.germ.summarydf, unname(table(f.unk.nchd[f.unk.nchd$col=="black",]$vclass)))
colnames(unk.germ.summarydf) <- names(table(f$vclass))
rownames(unk.germ.summarydf) <- c('chd', 'nchd')
print('## UNKNOWN: germline summary')
print(unk.germ.summarydf)

mo.summarydf <- rbind.data.frame(mo.summarydf, unk.mo.summarydf)
germ.summarydf <- rbind.data.frame(germ.summarydf, unk.germ.summarydf)

## write out summary tables
labs <- c('all_chd', 'all_nchd', 'syn_chd', 'syn_nchd', 'iso_chd', 'iso_nchd', 'unk_chd', 'unk_nchd')
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
mo.df <- rbind.data.frame(mo.df, naive_enrichment(f.chd, f.nchd))
# syndromic 
mo.df <- rbind.data.frame(mo.df, naive_enrichment(f.syn.chd, f.syn.nchd))
# isolated
mo.df <- rbind.data.frame(mo.df, naive_enrichment(f.non.chd, f.non.nchd))
# unknown
mo.df <- rbind.data.frame(mo.df, naive_enrichment(f.unk.chd, f.unk.nchd))
names(mo.df) <- c('in.ct', 'out.ct', 'in.LGD_Dmis', 'out.LGD_Dmis', 'in.Syn', 'out.Syn', 'or', 'p', 'ci1', 'ci2')
src <- c('all_chd', 'syn_chd', 'iso_chd', 'unk_chd')
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
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ(f.chd, f.nchd))
# syndromic 
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ(f.syn.chd, f.syn.nchd))
# isolated
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ(f.non.chd, f.non.nchd))
# unknown
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ(f.unk.chd, f.unk.nchd))
names(germ.df) <- c('in.ct', 'out.ct', 'in.LGD_Dmis', 'out.LGD_Dmis', 'in.Syn', 'out.Syn', 'or', 'p', 'ci1', 'ci2')
src <- c('all_chd', 'syn_chd', 'iso_chd', 'unk_chd')
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
  #labs <- paste(round(df$p, 2), paste(paste(df$v.in.d, df$v.in.nd, sep='/'), paste(df$v.out.d, df$v.out.nd, sep='/'), sep='\n'), sep='\n')
  labs <- paste(paste(paste(df$v.in.d, df$v.in.nd, sep='/'), paste(df$v.out.d, df$v.out.nd, sep='/'), sep='\n'), sep='\n')
  fp <- ggplot(df, aes(df$c2, df$or)) +  
    geom_text(aes(label=labs), vjust=-0.25, size=8) +
    geom_bar(stat='identity') +
    geom_hline(yintercept=1, linetype="longdash") + 
    labs(x="VAF bin upper bound", y="Enrichment" ) +
    #ylim(c(0, min(20, max(round(df$or, 0))+3))) +
    ylim(c(0, min(20, max(round(df$or, 0))+1))) +
    ggtitle(title) +
    theme(plot.title = element_text(size=24), axis.title = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
  return(fp)
}

## All samples
totdf <- NULL
src <- c('all.chd', 'all.chd', 'syn.chd', 'syn.chd', 'iso.chd', 'iso.chd', 'unk.chd', 'unk.chd')
# chd genes
vaf.chd.df <- vaf_bin_analysis(f.chd, f.nchd)
totdf <- rbind.data.frame(totdf, vaf.chd.df)
# plots
#p.vaf.chd.df <- plot_vaf_bin(vaf.chd.df, 'Plausible CHD Risk Genes')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
p.vaf.chd.df <- plot_vaf_bin(vaf.chd.df, 'CHD-related Genes')+ theme(axis.text.x = element_text(size=26), axis.text.y = element_text(size=26)) 
outplotall_dmis_name = paste(op, 'LGD_Dmis_vs_Syn_vaf.ALL.pdf', sep='.')
pdf(outplotall_dmis_name, width=7, height=7)
p.vaf.chd.df
dev.off()

## All samples EXCLUDING ISOLATED
totdf <- NULL
src <- c('all.chd', 'all.chd', 'syn.chd', 'syn.chd', 'iso.chd', 'iso.chd', 'unk.chd', 'unk.chd')
# chd genes
#vaf.chd.df <- vaf_bin_analysis(f.chd[f.chd$CHD_subtype!='isolated',], f.nchd[f.nchd$CHD_subtype!='isolated',])
vaf.chd.df <- vaf_bin_analysis(f.chd, f.nchd)
totdf <- rbind.data.frame(totdf, vaf.chd.df)
# plots
p.vaf.chd.df <- plot_vaf_bin(vaf.chd.df, 'CHD-related Genes')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotall_dmis_name = paste(op, 'LGD_Dmis_vs_Syn_vaf.ALL_noISO.pdf', sep='.')
pdf(outplotall_dmis_name, width=7, height=7)
p.vaf.chd.df
dev.off()

## Syndromic Samples
# chd genes
vaf.syn.chd.df <- vaf_bin_analysis(f.syn.chd, f.syn.nchd)
totdf <- rbind.data.frame(totdf, vaf.syn.chd.df)
# plots
p.vaf.syn.chd.df <- plot_vaf_bin(vaf.syn.chd.df, 'CHD-related Genes')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotsyn_dmis_name = paste(op, 'LGD_Dmis_vs_Syn_vaf.SYN.pdf', sep='.')
pdf(outplotsyn_dmis_name, width=7, height=7)
p.vaf.syn.chd.df
dev.off()

## Isolated Samples
# chd genes
vaf.iso.chd.df <- vaf_bin_analysis(f.non.chd, f.non.nchd)
totdf <- rbind.data.frame(totdf, vaf.iso.chd.df)
# plots
p.vaf.iso.chd.df <- plot_vaf_bin(vaf.iso.chd.df, 'CHD-related Genes')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotiso_dmis_name = paste(op, 'LGD_Dmis_vs_Syn_vaf.ISO.pdf', sep='.')
pdf(outplotiso_dmis_name, width=7, height=7)
p.vaf.iso.chd.df
dev.off()

## Unknown Samples
# chd genes
vaf.unk.chd.df <- vaf_bin_analysis(f.unk.chd, f.unk.nchd)
totdf <- rbind.data.frame(totdf, vaf.unk.chd.df)
# plots
p.vaf.unk.chd.df <- plot_vaf_bin(vaf.unk.chd.df, 'CHD-related Genes')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotunk_dmis_name = paste(op, 'LGD_Dmis_vs_Syn_vaf.UNK.pdf', sep='.')
pdf(outplotunk_dmis_name, width=7, height=7)
p.vaf.unk.chd.df
dev.off()

pdf(paste(op, 'LGD_Dmis_vs_Syn_vaf.combined.pdf', sep='.'), width=28, height=7)
grid.arrange(p.vaf.chd.df, p.vaf.syn.chd.df, p.vaf.iso.chd.df, p.vaf.unk.chd.df, ncol=4)
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
mo.df <- rbind.data.frame(mo.df, naive_enrichment_mis(f.chd, f.nchd))
# syndromic 
mo.df <- rbind.data.frame(mo.df, naive_enrichment_mis(f.syn.chd, f.syn.nchd))
# isolated
mo.df <- rbind.data.frame(mo.df, naive_enrichment_mis(f.non.chd, f.non.nchd))
# unknown
mo.df <- rbind.data.frame(mo.df, naive_enrichment_mis(f.unk.chd, f.unk.nchd))
names(mo.df) <- c('in.ct', 'out.ct', 'in.LGD_mis', 'out.LGD_mis', 'in.Syn', 'out.Syn', 'or', 'p', 'ci1', 'ci2')
src <- c('all_chd', 'syn_chd', 'iso_chd', 'unk_chd')
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
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ_mis(f.chd, f.nchd))
# syndromic 
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ_mis(f.syn.chd, f.syn.nchd))
# isolated
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ_mis(f.non.chd, f.non.nchd))
# isolated
germ.df <- rbind.data.frame(germ.df, naive_enrichment_germ_mis(f.unk.chd, f.unk.nchd))
names(germ.df) <- c('in.ct', 'out.ct', 'in.LGD_mis', 'out.LGD_mis', 'in.Syn', 'out.Syn', 'or', 'p', 'ci1', 'ci2')
src <- c('all_chd', 'syn_chd', 'iso_chd', 'unk_chd')
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
src <- c('all.chd', 'all.chd', 'syn.chd', 'syn.chd', 'iso.chd', 'iso.chd', 'unk.chd', 'unk.chd')
# chd genes
vaf.chd.df <- vaf_bin_analysis_mis(f.chd, f.nchd)
totdf <- rbind.data.frame(totdf, vaf.chd.df)
# plots
p.vaf.chd.df <- plot_vaf_bin_mis(vaf.chd.df, 'CHD-related Genes')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotall_mis_name = paste(op, 'LGD_MIS_vs_Syn_vaf.ALL.pdf', sep='.')
pdf(outplotall_mis_name, width=7, height=7)
p.vaf.chd.df
dev.off()

## Syndromic Samples
# chd genes
vaf.syn.chd.df <- vaf_bin_analysis_mis(f.syn.chd, f.syn.nchd)
totdf <- rbind.data.frame(totdf, vaf.syn.chd.df)
# plots
p.vaf.syn.chd.df <- plot_vaf_bin_mis(vaf.syn.chd.df, 'CHD-related Genes')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotsyn_mis_name = paste(op, 'LGD_MIS_vs_Syn_vaf.SYN.pdf', sep='.')
pdf(outplotsyn_mis_name, width=7, height=7)
p.vaf.syn.chd.df
dev.off()

## Isolated Samples
# chd genes
vaf.iso.chd.df <- vaf_bin_analysis_mis(f.non.chd, f.non.nchd)
totdf <- rbind.data.frame(totdf, vaf.iso.chd.df)
# plots
p.vaf.iso.chd.df <- plot_vaf_bin_mis(vaf.iso.chd.df, 'CHD-related Genes')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotiso_mis_name = paste(op, 'LGD_MIS_vs_Syn_vaf.ISO.pdf', sep='.')
pdf(outplotiso_mis_name, width=7, height=7)
p.vaf.iso.chd.df
dev.off()

## Unknown Samples
# chd genes
vaf.unk.chd.df <- vaf_bin_analysis_mis(f.unk.chd, f.unk.nchd)
totdf <- rbind.data.frame(totdf, vaf.unk.chd.df)
# plots
p.vaf.unk.chd.df <- plot_vaf_bin_mis(vaf.unk.chd.df, 'CHD-related Genes')+ theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) 
outplotunk_mis_name = paste(op, 'LGD_MIS_vs_Syn_vaf.UNK.pdf', sep='.')
pdf(outplotunk_mis_name, width=7, height=7)
p.vaf.unk.chd.df
dev.off()

pdf(paste(op, 'LGD_MIS_vs_Syn_vaf.combined.pdf', sep='.'), width=28, height=7)
grid.arrange(p.vaf.chd.df, p.vaf.syn.chd.df, p.vaf.iso.chd.df, p.vaf.unk.chd.df, ncol=4)
dev.off()

totdf$src <- src
lgd_mis_vaf_name <- paste(op, 'LGD_MIS_vs_Syn_vaf.txt', sep='.')
write.table(totdf, lgd_mis_vaf_name, sep='\t', quote=F, row.names=F)

print("Done !")
print(paste("output", c(paste(op, 'mosaic_summary.txt', sep='.'), paste(op, 'germline_summary.txt', sep='.'), lgd_dmis_vaf_name, outplotall_dmis_name, outplotsyn_dmis_name, outplotiso_dmis_name, outplotunk_dmis_name, lgd_mis_vaf_name, outplotall_mis_name, outplotsyn_mis_name, outplotiso_mis_name, outplotunk_mis_name), sep=' : '))

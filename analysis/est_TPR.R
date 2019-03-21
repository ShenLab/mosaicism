f <- read.table('~/Dropbox (Personal)/Alex/mosaic_manuscript/pcgc1115/pcgc1115.denovo.anno.txt', sep='\t', header=T, quote='"')

tmp <- f[f$post>1 & f$post<10 & f$vaf<0.5,]

pdf('TPR_sites_below_10.pdf')
p.title <- paste(nrow(tmp), "sites 1<posterior odds<10", sep=' ')
hist(tmp$post, breaks=seq(0, 10,  by=1), main=paste("Posterior Odds", p.title, sep='\n'), xlab="posterior odds", cex.axis=1.5, cex.lab=1.5)
dev.off()

# estimated true fraction for posterior odds < 10
est_tp <- function(v){
  c1 <- seq(1, 9, by=1)
  c2 <- seq(2, 10, by=1)
  mofrac <- 0.1212 # from EM estimate 08/01/2018 after excluding sites failing IGV review
  df <- NULL
  for(i in 1:length(c1)){
    #print(c1[i])
    #ratio <- c1[i]/mofrac # post_odds/prior = LR
    #est.tprate <- ratio/(ratio+1) # est TP = LR/(1+LR), FP = 1-TP
    est.tprate <- c1[i]/(1+c1[i])
    bin.ct <- nrow(v[v$post > c1[i] & v$post <= c2[i],])
    est.tp <- est.tprate*bin.ct
    est.fp <- (1-est.tprate)*bin.ct
    out <- c(c1[i], c2[i], est.tprate, bin.ct, est.tp, est.fp)
    df <- rbind.data.frame(df, out)
  }
  names(df) <- c('c1', 'c2', 'est.TPrate', 'bin.ct', 'est.TP', 'est.FP')
  #print(paste('Est.TP',round(sum(df$est.TP), 1), sep=' : '))
  #print(paste('Est.FP',round(sum(df$est.FP), 1), sep=' : '))
  return(df)
}
est.df <- est_tp(tmp)
est.df

tmp$est_tpr <- tmp$post/(1+tmp$post)
pdf('TPR_est_FN.pdf')
tmp.est_prop <- paste(paste('Est. TP', round(sum(est.df$est.TP), 1), sep=':'), paste('Est. FP', round(nrow(tmp)-sum(est.df$est.TP), 1), sep=':'), sep='\n')
hist(tmp$est_tpr, main=paste(p.title ,tmp.est_prop, sep='\n'), xlab="TPR", cex.axis=1.5, cex.lab=1.5)
dev.off()

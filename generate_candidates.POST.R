library("bbmle")
library("emdbook") 
library("ggplot2") 


#############################################
## FUNCTIONS
#############################################
## Function to run Expectation Maximization to decompose VAF distribution into Mosaic and Germline
## Output: Mosaic Fraction, Mean Mosaic VAF, Mean Germline VAF
mofracEM <- function(x, thetahat, op, printopt) {
  idx=1
  d = x$altdp
  n = x$refdp+x$altdp
  ind = rep(0,nrow(x))
  maxi = 1e6
  delta=0.0001  # convergence condition for r estimate
  pp2 = 0.6 # prior germline mean VAF
  pp1 = 0.05 # prior mosaic mean VAF
  pr = 0.5 # prior mosaic fraction
  while (idx < maxi) {
    # E step
    for (i  in 1:nrow(x)) {
      j = d[i]
      k = n[i]
      mo.afs = seq(0.05, 0.4, by=0.05)
      l0 = dbetabinom(j, p=pp2, size=k, theta=thetahat) * (1-pr)
      tmp.l1 = rep(0, length(mo.afs))
      for (z in 1:length(mo.afs)){
        tmp.l1[z] = dbetabinom(j, p=mo.afs[z], size=k, theta=thetahat)
      }
      l1 = mean(tmp.l1) * pr
      if (l1 > l0) {
        ind[i] = 1	
      } else {
        ind[i] = 2
      }
    }
    # M step
    nr = length(ind[ind==1])/nrow(x)
    np1 = mean(d[ind==1]/n[ind==1])
    np2 = max(mean(d[ind==2]/n[ind==2]), 0.49)
    if (abs(nr-pr) < delta) # converge
    {	
      break
      
    } else
    {
      pr = nr
      pp1 = np1
      pp2 = np2
    }
    idx = idx + 1  # avoid infinite loops
  }
  print(c(length(ind[ind==1]), nrow(x), nr, np1, np2))
  # print EM plot only for 2nd run (after excluding likely mosaics)
  if(printopt == 'yes'){
    em.p = paste(op, 'EM.pdf', sep='.')
    pdf(em.p)
    p.title = paste('Variant Allele Fraction', paste('Est. Mosaic Fraction', round(nr, 4), sep=' = '), sep='\n')
    hist(d/n, br=seq(0.0, 1.0, by=0.05), freq=F, ylim=c(0, max(density(d[ind==2]/n[ind==2])$y)+1), main=p.title, xlab="Variant Allele Fraction (VAF)", cex.axis=1.5, cex.lab=1.5)
    mo.dens <- density(d[ind==1]/n[ind==1])$y * (length(d[ind==1])/length(d)) # adjust density for mosaics by mosaic fraction
    germ.dens <- density(d[ind==2]/n[ind==2])$y * (length(d[ind==2])/length(d)) # adjust density for germline by germline fraction
    lines(density(d[ind==1]/n[ind==1])$x, mo.dens,  col='red', lwd=2)
    lines(density(d[ind==2]/n[ind==2])$x, germ.dens,  col='blue', lwd=2)
    legend("topright", c('germline', 'mosaic'), col=c('blue', 'red'), lty=1, cex=1)
    dev.off()
  }
  return(c(nr, np1, np2))
}

## Function to find minimum number of alternate read support, given DP and Expected FP 0.01 and Exome size 3e7
find_min_alt <- function(N){
  out.min = 6
  for(i in 1:(N/2)){
    exp.fp <- (1-ppois(i, N*(0.005/3)))*3e7
    if(exp.fp <= 0.01){
      out.min = i - 1
      break
    }
  }
  return(out.min)
}

## Function to remove outlier samples with abnormally high denovo counts
## 1. calculate expected number of samples with 0, 1, 2, ... denovo counts given a mean(# denovos/sample) and cohortsize
## 2. cutoff = denovo count value for which expected number of samples is <1
## 3. remove all samples with denovo counts > cutoff as outliers
remove_outliers <- function(x, cohortsize, op){
  samples <- names(table(x$id))
  counts <- unname(table(x$id))
  df <- cbind.data.frame(samples, counts)
  exp.counts <- (1-ppois(0:max(max(counts), 20), mean(counts)))*cohortsize
  names(exp.counts) <- c(seq(0, max(max(counts), 20), by=1))
  cutoff <- as.numeric(names(exp.counts[exp.counts<1])[1])
  df.outliers <- df[df$Freq>cutoff,]
  outliers <- c(as.character(df.outliers$samples))
  outlier.frac <- paste(length(outliers), cohortsize, sep='/')
  print(paste("Mean denovo count", mean(counts), sep=' : '))
  print(paste("Max denovo count cutoff", cutoff, sep=' : '))
  print(paste("Outlier samples", outlier.frac, sep=' : '))
  outfile <- paste(op, 'outlier_samples.txt', sep='.')
  write.table(df.outliers[,c(1,3)], outfile ,quote=F, row.names=F, sep='\t')
  print(paste("# outlier samples written to", outfile, sep=' : '))
  return(outliers)
}

## Function to calculate Likelihood Ratio
llratio <- function(alt, dp, th, mp){
  phat <- alt/dp
  px <- dbetabinom(alt, prob=phat, size=dp, theta=th) # P(D|Mx)
  p1 <- dbetabinom(alt, prob=mp, size=dp, theta=th) # P(D|M1)
  p0 <- dbetabinom(alt, prob=0, size=dp, theta=th) # P(D|M0)
  x_1 <- px/p1
  x_0 <- px/p0 # irrelevant, since p0 = 0
  return(x_1)
}

## Function to estimate Theta given N and Nalt
betaBinomEst <- function(counts, total) {
  mtmp <- function(prob,size,theta) { -sum(dbetabinom(counts,prob,size,theta,log=TRUE)) } 
  m0 <- mle2(mtmp,start=list(prob=0.5,theta=10),data=list(size=total))
  # MLE of theta
  t = coef(m0)["theta"]
  return(t)
}

## Function to estimate theta for entire variant callset
estimateTheta <- function(x) {
  mn = min(max(x$N), 500)
  results = matrix(0, nrow=length(levels(factor(x$N))) , ncol=6)
  results = as.data.frame(results)
  colnames(results) = c("N", "m", "p", "var", "theta","nalt")
  i = 1
  for (t in 1:mn) {
    v = x[x$N == t, ]$altdp
    if (length(v) > 3 ) { 
      theta = betaBinomEst(v, t)[1]
      results[i,1] = t
      results[i,2] = length(v)
      results[i,3] = mean(v)/t
      results[i,4] = var(v)
      results[i,5] = theta  
      results[i,6] = mean(v) # mean of nalts of variants used to support this
      
      i = i + 1
    }
  }
  return(results)
}

## Function to calculate Beta-Binomial p-value 
testBetaBinomial <- function(x, mp , th) {
  ### now given p and theta, test each variant against a null. The alternative hypothesis is that the variant is a mosaic
  pvalues = rep(0, nrow(x))
  for ( i in 1:nrow(x)) {
    p = sum(  dbetabinom(1:x[i,]$altdp,   prob= mp, size = x[i,]$N,  theta = th ) )
    pvalues[i] = p
  }
  return(pvalues) 
}

## Main Function
fitReadCounts <- function(a, op) {
  a$refdp <- as.numeric(a$refdp)
  a$altdp <- as.numeric(a$altdp) 
  N = a$refdp + a$altdp 
  x = cbind.data.frame(a, N)
  x$vaf <- x$altdp/x$N
  ## save record of all denovos filtered due to outliers, fdr.min.nalt, etc
  all <- x
  ## save record of all denovos filtered due to outliers, fdr.min.nalt, etc
  all$filter2 <- rep('.', nrow(all))
  
  x<-x[!is.na(x$N),]
  x <- x[x$refdp>0&x$altdp>0,] # ignore 0 ref or alt DP entries
  
  if("filter" %in% colnames(x)){
    x <- x[x$filter==".",] # remove variants failing QC
  }
  
  if("IGV" %in% colnames(x)){
    x <- x[x$IGV %in% c('ok', '.'),] # remove variants failing IGV inspection
  }

  
  print(paste("# denovos passing filters & nonzero AD", paste(nrow(x), nrow(a), sep='/'), sep=':'))
  
  tmp.fdr.min.alt <- mapply(find_min_alt, all$N)
  all = cbind.data.frame(all, tmp.fdr.min.alt)
  
  fdr.min.alt <- mapply(find_min_alt, x$N)
  x = cbind.data.frame(x, fdr.min.alt)
  old.total <- nrow(x)
  
  if(nrow(all[all$altdp <= all$tmp.fdr.min.alt,]) > 1){
    all[all$altdp <= all$tmp.fdr.min.alt,]$filter2 <- 'fdr.min.nalt'
  }


    
  fdrname <- paste(op, 'fdr_min_nalt.pdf', sep='.')
  pdf(fdrname)
  df <- NULL
  for(d in seq(10, 500, by=1)){
    out <- c(d, find_min_alt(d))
    df <- rbind.data.frame(df, out)
  }
  names(df) <- c('dp', 'fdr.min.nalt')
  plot(x$N, x$altdp, col='black', type='p', cex=0.35, main="FDR-based min. Nalt", xlab="DP", ylab="Nalt", xlim=c(0, 500), ylim=c(0, 60), cex.axis=1.5, cex.lab=1.5)
  lines(df$dp, df$fdr.min.nalt, type='s', col='red' )
  legend("topright", c("FDR-based min. Nalt"), col=c('red'), lty=1, cex=0.8)
  dev.off()
  
  x <- x[x$altdp>x$fdr.min.alt,]
  
  print(paste("# sites surviving FDR-based min Nalt cutoff", paste(nrow(x), old.total, sep='/'), sep=':'))
  
  ## remove bad samples 
  badsamples <- remove_outliers(x, cohortsize, op)
  x <- x[!(x$id %in% badsamples),]
  
  #all[all$id %in% badsamples,]$filter2 <- 'outlier'
  if(nrow(x[(x$id %in% badsamples),]) > 1){
    all[all$id %in% badsamples,]$filter2 <- 'outlier'
  }
  
  print(paste("# denovos surviving bad sample removal", paste(nrow(x), nrow(a), sep='/'), sep=':'))
  
  ## write out filter log for all denovos
  allname <- paste(op, '.all_denovos.txt', sep='')
  write.table(all, allname, sep='\t', row.names=F, quote=F)

  ## generate depth table
  results = estimateTheta(x)  
  results = results[results$m  > 0, ]
  
  ## theta and p estimations 
  if(unname(table(results$m>20))[1]<=0.8*nrow(results)){ # if more than 20% of DP values in results have >20 variants supporting, use m>20 (only data from DP values with enough support)
    thetahat = sum(results[results$N > 12 & results$m > 20, ]$m*results[results$N > 12 & results$m > 20, ]$theta)/sum(results[results$N > 12 & results$m > 20, ]$m)
  }else { # if less than 20% of DP values in results have >20 variants supporting, use m>0 (all data available)
    thetahat = sum(results$m*results$theta)/sum(results$m)
  }
  
  ## initial EM estimation of mosaic fraction
  EMout = mofracEM(x, thetahat, op, 'no')
  mofrac = EMout[1]
  movaf = EMout[2]
  germvaf = EMout[3]
  mp = 0.49
  lrcut = postcut/mofrac # use LR cutoff corresponding to posterior odds of postcut, given EM mosaic fraction

  ## calculate p values for candidates
  pvalues <- testBetaBinomial(x, mp, thetahat)
  padj <- p.adjust(pvalues, method="BH")
  x = cbind.data.frame(x, pvalues)
  x = cbind.data.frame(x, padj)
  
  ##  generate candidates
  x$lr <- mapply(llratio, x$altdp, x$N, thetahat, mp)
  x$post <- x$lr * mofrac # posterior odds = LR * prior = LR * mosaic fraction
  x$col = 'black' # for later plots
  x[x$post>postcut & x$altdp/x$N<=0.5,]$col='red'
  z = x[x$post>postcut & x$altdp/x$N<=0.5,] 
  
  ## Exclude likely mosaic sites and re-estimate parameters
  x.non <- x[x$col=="black",]
  results.non <- estimateTheta(x.non)
  results.non = results.non[results.non$m  > 0, ]
  if(unname(table(results.non$m>20))[1]<=0.8*nrow(results.non)){
    thetahat.non = sum(results.non[results.non$N > 12 & results.non$m > 20, ]$m*results.non[results.non$N > 12 & results.non$m > 20, ]$theta)/sum(results.non[results.non$N > 12 & results.non$m > 20, ]$m)
  } else {
    thetahat.non = sum(results.non$m*results.non$theta)/sum(results.non$m)
  }
  EMout2 = mofracEM(x, thetahat.non, op, 'yes')
  mofrac = EMout2[1]
  movaf = EMout2[2]
  germvaf = EMout2[3]
  lrcut = postcut/mofrac # use LR cutoff corresponding to posterior odds of 5, given EM mosaic fraction
  mp.non = 0.49
  print(c(mofrac, movaf, germvaf))
  print(c(paste("# mp", mp.non, sep=": "), paste("thetahat", thetahat.non, sep=": ")))
  print(paste("# LR cutoff", lrcut, sep=': '))

  pvalues <- testBetaBinomial(x, mp.non, thetahat.non) 
  padj <- p.adjust(pvalues, method="BH") # add column for adjusted p.value
  
  ##  generate updated set of candidates
  x$pvalues = pvalues
  x$padj = padj
  x$lr <- mapply(llratio, x$altdp, x$N, thetahat.non, mp.non)
  x$post <- x$lr * mofrac # posterior odds = LR * p(0/1*)/p(0/1)
  x[x$post>postcut & x$altdp/x$N<=0.5,]$col='red'
  z = x[x$post>postcut & x$altdp/x$N<=0.5,]
  print(paste('# Mosaic candidates', nrow(z), sep=': '))
  print(paste('# Estimated mosaic fraction', nrow(z)/nrow(x), sep=': '))

  cname = paste(op, 'candidates.txt', sep='.')
  c.pname = paste(op, 'p_candidates.txt', sep='.')
  dnname = paste(op, 'denovo.txt', sep='.')
  write.table(z, cname, quote=F,row.names=F, sep="\t")
  print(paste('####### OUTPUT candidates posterior odds > cutoff -- see', cname, sep=': '))
  write.table(x, dnname, quote=F,row.names=F, sep="\t")
  print(paste('####### OUTPUT annotated denovos -- see', dnname, sep=': '))

  #############################################
  ## OUTPUT PLOTS
  ############################################# 
  ## VAF vs. LR
  p1name = paste(op, 'vaf_vs_post.pdf', sep='.')
  pdf(p1name)
  plot(x$altdp/x$N, log10(x$post), xlab="VAF", ylab="log10(posterior odds)",  cex=0.5, xlim=c(0,1), cex.axis=2, cex.lab=2)
  title(main="VAF vs. posterior odds", cex.main=2)
  lines(z$altdp/z$N, log10(z$post), type="p", col="red", cex=0.5)
  abline(h=log10(postcut), col="red")
  abline(v=0.35, col="red")
  legend("topright", c("all variants", "candidate mosaics"), pch=1, col=c("black", "red"), cex=1.25)
  dev.off()
  
  ## DP vs. VAF
  p2name = paste(op, 'dp_vs_vaf.pdf', sep='.')
  pdf(p2name)
  plot(x$refdp + x$altdp, x$altdp/ (x$altdp+x$refdp),  xlab = "DP", ylab="VAF", pch=16, cex=0.8, ylim=range(c(0.05, 0.95)), xlim=c(0, 500), col=adjustcolor(x$col, alpha=0.2), cex.axis=2, cex.lab=2)
  title(main='DP vs. VAF', cex.main=2)
  mn = min(max(x$N), 500)
  ci = matrix(0, nrow = mn, ncol=3)
  for (i in 10:mn) {
    lci = 0
    hci = mn 
    d = 0
    for (j in 1:i) {
      d = sum(dbetabinom(0:j, prob = mp, size = i, theta = thetahat.non  ) )
      if (d > 0.025) {
        lci = j - 1
        break
      } 
    }   
    for (j in i:1) {
      d = sum(dbetabinom(i:j, prob = mp, size = i, theta = thetahat.non  ) )
      if (d > 0.025) {
        hci = j + 1
        break
      } 
    }
    ci[i, 1] = lci
    ci[i, 2] = hci
    ci[i, 3] = i
  }
  lines(1:mn, ci[,1]/ci[,3], col='red')
  lines(1:mn, ci[,2]/ci[,3], col='red')
  abline(h=mp, col='blue')
  legend("topright", c("germline de novo", "candidate mosaic"), col=c("black", "red"), pch=16, cex=1.25)
  dev.off()
  
  ## overdispersion
  p3name = paste(op, 'overdispersion.pdf', sep='.')
  pdf(p3name)
  plot(results.non$N, results.non$var,  xlab = "DP", ylab = "Var(Nalt)", log="y", cex=0.5, main="DP vs. Var(Nalt) \n Binomial vs. Beta-Binomial", cex.axis=1.5, cex.lab=1.5)
  lines(results.non$N, results.non$N * results.non$p * (1 - results.non$p), col='blue')
  lines(results.non$N, results.non$N * mp * (1 - mp) * (results.non$N + thetahat.non) / (thetahat.non + 1), col='red')
  legend("bottomright", c("Binomial", "Beta-Binomial"), col=c("blue", "red"), lty=1)
  dev.off()
  
  ## QQ plot
  qq.p = paste(op, 'QQ.pdf', sep='.')
  pdf(qq.p)
  obs.p <- x$pvalues
  n = nrow(x)
  exp.p <- seq(from=1/(n+1), to=n/(n+1), by= 1/(n+1) ) # expected: n/n+1 .... 1/n+1
  obs.p.log <- -log10(obs.p)
  exp.p.log <- -log10(exp.p)
  plot(sort(exp.p.log, decreasing=TRUE), sort(obs.p.log, decreasing=TRUE), main="QQ plot", xlab="exp -log10(p)", ylab="obs -log10(p)", xlim=c(0, 1), ylim=c(0, 1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
  lines(sort(exp.p.log, decreasing=TRUE), sort(exp.p.log, decreasing=TRUE), col='red')
  dev.off()
  
  pnamelist = c(fdrname, p1name, p2name, p3name, qq.p, paste(op, "EM.png", sep='.'))
  print(paste('####### OUTPUT PLOT', pnamelist, sep=': '))
  
  return(results)
}

#############################################
## HANDLE ARGUMENTS
#############################################
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4){
  stop("Please provide 4 arguments: <ADfile> <outfile prefix> <posterior odds cutoff> <cohortsize>")
}
fname <- toString(args[1])
outprefix <- toString(args[2])
postcut <- as.integer(args[3])
cohortsize <- as.integer(args[4])
print(c(paste("Input File", fname, sep=': ') , paste("Outfile prefix",outprefix, sep=': '), paste("Posterior Odds Cutoff",postcut, sep=': ')))

#############################################
## LOAD DATA AND RUN MAIN FUNCTION
#############################################
a <- read.table(fname, sep='\t', header=T, quote='"')
results = fitReadCounts(a, outprefix)
print("################ DONE !")

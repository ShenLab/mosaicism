library("bbmle")
library("emdbook") 
library("ggplot2") 
library("fitdistrplus")


#############################################
## FUNCTIONS
#############################################
## Function to run Expectation Maximization to decompose VAF distribution into Mosaic and Germline
## Input: dataframe of variants, thetahat estimate, output plot file prefix, flag for yes/no printing plot
## Output: Dataframe with Mosaic Fraction, Mean Mosaic VAF, Mean Germline VAF, min mosaic VAF, max mosaic VAF, all VAF, index (1=mosaic, 2=germline)
mofracEM <- function(x, thetahat, op, printopt) {
  idx <- 1 # loop index 
  maxi <- 1e6 # max index
  delta <- 0.0001  # convergence condition for r estimate
  
  x$ind <- 0 # Z; stores latent label; 1 for mosaic, 2 for germline
  
  pp2 = 0.6 # prior germline mean VAF
  pp1 = 0.01 # prior mosaic mean VAF
  pr = 0.5 # prior mosaic fraction
  
  while (idx < maxi) {
    # E step
    
    
    ## version 2: calculate posterior odds and compare
    x$l0 <- dbetabinom(x$altdp, p=pp2, size=x$N, theta=thetahat) 
    x$l1 <- dbetabinom(x$altdp, p=pp1, size=x$N, theta=thetahat) 
    x$lr <- x$l1/x$l0
    x$post <- x$lr * pr/(1-pr)
    
    x$pp <- x$post/(1+x$post)
    
    x[x$post>1,]$ind <- 1
    x[x$post<=1,]$ind <- 2
    
    
    # M step
    nr = nrow(x[x$ind==1,])/nrow(x)
    np1 = mean(x[x$ind==1,]$altdp/x[x$ind==1,]$N)
    #np2 = max(mean(x[x$ind==2,]$altdp/x[x$ind==2,]$N), 0.49)
    np2 = mean(x[x$ind==2,]$altdp/x[x$ind==2,]$N)
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
  options(scipen=999)
  
  print('## EM results:')
  print(paste('## mosaic count', nrow(x[x$ind==1,]), sep=': '))
  print(paste('## all denovos count', nrow(x), sep=': '))
  print(paste('## prior mosaic fraction', nr, sep=': '))
  print(paste('## mosaic mean VAF', np1, sep=': '))
  print(paste('## germline mean VAF', np2, sep=': '))
  print(paste('## sum of posterior prob approach', sum(x$pp)/nrow(x), sep=': '))
  print(paste('## true fraction', nrow(x[x$src=='mosaic',])/nrow(x), sep=': '))
  
  ## get limits of integration for BF calculation
  minvaf <- min(x[x$ind==1,]$vaf)
  maxvaf <- max(x[x$ind==1,]$vaf)
  
  # print EM plot only for 2nd run (after excluding likely mosaics)
  if(printopt == 'yes'){
    em.p = paste(op, 'EM.pdf', sep='.')
    pdf(em.p)
    true_frac <- nrow(x[x$src=='mosaic',])/nrow(x)
    #p.title = paste('Variant Allele Fraction', paste('Est. Mosaic Fraction', round(nr, 4), sep=' = '),  paste('True Fraction', round(true_frac, 4), sep=' = '), sep='\n')
    p.title = paste('Variant Allele Fraction', paste('Est. Mosaic Fraction', round(sum(x$pp)/nrow(x), 4), sep=' = '),  paste('True Fraction', round(true_frac, 4), sep=' = '), sep='\n')
    hist(x$altdp/x$N, br=seq(0.0, 1.0, by=0.05), freq=F, ylim=c(0, max(density(x[x$ind==2,]$altdp/x[x$ind==2,]$N)$y)+1), main=p.title, xlab="Variant Allele Fraction (VAF)", cex.axis=1.5, cex.lab=1.5)
    mo.dens <- density(x[x$ind==1,]$altdp/x[x$ind==1,]$N)$y * (nrow(x[x$ind==1,])/nrow(x)) # adjust density for mosaics by mosaic fraction
    germ.dens <- density(x[x$ind==2,]$altdp/x[x$ind==2,]$N)$y * (nrow(x[x$ind==2,])/nrow(x)) # adjust density for germline by germline fraction
    lines(density(x[x$ind==1,]$altdp/x[x$ind==1,]$N)$x, mo.dens,  col='red', lwd=2)
    lines(density(x[x$ind==2,]$altdp/x[x$ind==2,]$N)$x, germ.dens,  col='blue', lwd=2)
    legend("topright", c('germline', 'mosaic'), col=c('blue', 'red'), lty=1, cex=1)
    dev.off()
  }
  
  outmovaf <- x[x$ind==1,]$vaf
  outgermvaf <- x[x$ind==2,]$vaf
  outdf <- data.frame(mofrac = rep(sum(x$pp)/nrow(x), nrow(x)), movaf = rep(np1, nrow(x)), germvaf = rep(np2, nrow(x)), minvaf = rep(minvaf,nrow(x)), maxvaf = rep(maxvaf, nrow(x)), allvaf = x$vaf, ind=x$ind)
  return(outdf) 
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
  print(paste("# Mean denovo count", mean(counts), sep=' : '))
  print(paste("# Max denovo count cutoff", cutoff, sep=' : '))
  print(paste("# Outlier samples", outlier.frac, sep=' : '))
  outfile <- paste(op, 'outlier_samples.txt', sep='.')
  write.table(df.outliers[,c(1,3)], outfile ,quote=F, row.names=F, sep='\t')
  print(paste("# outlier samples written to", outfile, sep=' : '))
  return(outliers)
}


## Function to calculate Likelihood Ratio
likelihood_ratio <- function(alt, dp, th, mp){
  phat <- alt/dp
  px <- dbetabinom(alt, prob=phat, size=dp, theta=th) # P(D|Mx)
  p1 <- dbetabinom(alt, prob=mp, size=dp, theta=th) # P(D|M1)
  p0 <- dbetabinom(alt, prob=0, size=dp, theta=th) # P(D|M0)
  x_1 <- px/p1
  x_0 <- px/p0 # irrelevant, since p0 = 0
  return(x_1)
}


## Function to calculate Bayes Factor -- version 2 (directly sampling from posterior mosaic VAF distribution)
## a, b correspond to the posterior Beta(a, b)
bayes_factor <- function(alt, dp, th, mp, a, b, n_sample, sample_avg_dp){
  
  if(alt/dp > 0.5){ ## BF is 0 for sites VAF>0.5
    px_w <- 0
  } else{
    ## sample VAFs from the posterior VAF distribution (Beta(alpha+nalt, beta+nref))
    #n_test <- 50
    set.seed(1)
    #vaftest <- rbeta(n_test, shape1=a+alt, shape2=b+(dp-alt))
    vaftest <- get_vaf_dist(a, b, n_sample, sample_avg_dp)
    
    
    ## avg by the number of tests
    px <- mapply(dbetabinom, alt, vaftest, dp, th) * 1/length(vaftest)
    px_w <- sum(px) # weighted sum of P(D|Mx) (essentially avg)
  }
  
  p1 <- dbetabinom(alt, prob=mp, size=dp, theta=th) # P(D|M1)
  
  bf <- px_w / p1
  
  return(bf)
  
}

## Function to estimate Theta given N and Nalt
betaBinomEst <- function(counts, total) {
  mtmp <- function(prob,size,theta){
    -sum(dbetabinom(counts,prob,size,theta,log=TRUE)) 
  } 
  m0 <- mle2(mtmp,start=list(prob=0.5,theta=10),data=list(size=total))
  # MLE of theta
  t = coef(m0)["theta"]
  return(t)
}

## Input: dataframe of variants and a single DP value
## Output: vector containing c(dp, value, # variants, mean VAF, var Nalt, theta, mean Nalt)
est_params <- function(df, dp){
  v <- df[df$N == dp,]$altdp
  
  out  <- c(dp, 0, 0, 0, 0, 0) # default value
  
  if(length(v) > 3){
    theta <- betaBinomEst(v, dp)[1]
    
    nvar <- length(v)
    vaf_mean <- mean(v)/dp
    nalt_var <- var(v)
    th <- unname(theta)
    nalt_mean <- mean(v)
    
    out <- c(dp, nvar, vaf_mean, nalt_var, th, nalt_mean)
  }
  return(out)
}

## Function to estimate theta for entire variant callset -- calls est_params()
## Input: dataframe of variants
## Output: dataframe of results
estimateTheta <- function(df){
  mn <- min(max(df$N), 500)
  dps <- seq(1, mn, by=1)
  results <- data.frame(t(sapply(dps, est_params, df=df))) ## run est_params for each value in dps
  names(results) <- c('N', 'm', 'p', 'var', 'theta','nalt')
  
  results <- results[results$m>3,] # only return results in cases where there were enough variants to estimate parameters
  return(results)
}

## Function to calculate Beta-Binomial p-value
## Input: dataframe of variants, VAF of null germline model, theta estimate
## Output: vector of pvalues, 1 for each variant in dataframe
testBetaBinomial <- function(df, mp, th){
  ## function to be applied to each row of dataframe
  get_p <- function(altdp, N, ...){
    extra.args <- list(...)
    mp <- extra.args[[1]]
    th <- extra.args[[2]]
    sum(dbetabinom(1:altdp, prob=mp, size=N, theta=th))
  }
  pvalues <- mapply(get_p, df$altdp, df$N, mp, th)
  return(pvalues)
}


## Function to fit Beta() distribution to observed mosaic VAF distribution
## output plot summarizing fit
## choose best parameter combination on the basis of K-S distance to observed VAF
fit_mosaic <- function(EMvaf, sample_avg_dp, op){
  alpha <- seq(0.1, 1, by=0.1)
  beta <- seq(5, 15, by=0.5)
  
  n <- 500
  #n <- 1000 # for low coverage datasets
  sample_avg_dp <- sample_avg_dp
  df <- data.frame()
  for(i in 1:length(alpha)){
    for(j in 1:length(beta)){
      set.seed(1)
      tmpvaf0 <- rbeta(n, shape1=alpha[i], shape2=beta[j])/2
      tmpdp <- rnbinom(n, size=4, mu=sample_avg_dp)
      tmpnalt <- round(tmpvaf0 * tmpdp) # round Nalt to whole numbers
      tmpvaf <- tmpnalt/tmpdp
      
      tmpdf <- cbind.data.frame(tmpvaf0, tmpvaf, tmpdp, tmpnalt)
      tmpdf <- tmpdf[tmpdf$tmpnalt >= 6,]
      
      tmp.kstest <- ks.test(tmpdf$tmpvaf, EMvaf)
      ksD <- unname(tmp.kstest$statistic)
      ksp <- tmp.kstest$p.value
      out <- c(alpha[i], beta[j], ksD, ksp)
      df <- rbind.data.frame(df, out)
    }
  }
  
  names(df) <- c('a', 'b', 'ks.d', 'ks.p')
  df[order(df$ks.d),][1:5,]
  
  
  
  ## print fitting results
  print('## Mosaic distribution fitting ')
  msg.fit <- paste( paste('alpha', df[order(df$ks.d),][1,]$a, sep='='), 
                    paste('beta', df[order(df$ks.d),][1,]$b, sep='='), 
                    paste('KS.dist', df[order(df$ks.d),][1,]$ks.d, sep='='),
                    paste('KS.pvalue', df[order(df$ks.d),][1,]$ks.p, sep='='), sep='   ')
  print(msg.fit)
  
  ## plot fit
  test <- get_vaf_dist(df[order(df$ks.d),][1,]$a, df[order(df$ks.d),][1,]$b, n, sample_avg_dp)
  pdf(paste(op, 'mosaic_fit.pdf', sep='.'))
  t.fit <- paste( paste(paste('alpha', df[order(df$ks.d),][1,]$a, sep='='), paste('beta', df[order(df$ks.d),][1,]$b, sep='='), sep=' '), 
                  paste('KS.dist', df[order(df$ks.d),][1,]$ks.d, sep='='),
                  paste('KS.pvalue', df[order(df$ks.d),][1,]$ks.p, sep='='), sep='\n')
  
  hist(EMvaf, main=t.fit, xlab='Variant Allele Fraction', freq=F)
  lines(density(test), col='red')
  legend('topright', c('fitted distribution'), col=c('red'), lty=1)
  dev.off()
  
  best.alpha <- df[order(df$ks.d),][1,]$a
  best.beta <- df[order(df$ks.d),][1,]$b
  
  return(c(best.alpha, best.beta))
}

## generate VAFs for plotting
get_vaf_dist <- function(a, b, n, sample_avg_dp){
  set.seed(1)
  tmpvaf0 <- rbeta(n, shape1=a, shape2=b)/2
  tmpdp <- rnbinom(n, size=4, mu=sample_avg_dp)
  tmpnalt <- round(tmpvaf0 * tmpdp) # round Nalt to whole numbers
  tmpvaf <- tmpnalt/tmpdp
  
  tmpdf <- cbind.data.frame(tmpvaf0, tmpvaf, tmpdp, tmpnalt)
  tmpdf <- tmpdf[tmpdf$tmpnalt >= 6,]
  
  return(tmpdf$tmpvaf)
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
  print("Running E-M (iteration 1)")
  EMout = mofracEM(x, thetahat, op, 'no')
  mofrac = EMout[1,]$mofrac
  movaf = EMout[1,]$movaf
  germvaf = EMout[1,]$germvaf
  mp = germvaf
  vmin = EMout[1,]$minvaf
  vmax = EMout[1,]$maxvaf
  tmp.mo.vaf = EMout[EMout$ind==1,]$allvaf
  tmp.germ.vaf = EMout[EMout$ind==2,]$allvaf
  
  ## fit mosaic distribution to estimate parameters
  #fit.mo <- fitdist(tmp.mo.vaf, 'beta') # shape1=6.655298, shape2=45.098363
  fit.mo <- fit_mosaic(tmp.mo.vaf, sample_avg_dp, outprefix)
  fit.mo.shape1 <- fit.mo[1]
  fit.mo.shape2 <- fit.mo[2]
  
  ## fit germline distribution to estimate parameters
  fit.germ <- fitdist(tmp.germ.vaf, 'beta') # shape1=16.03405, shape2=16.54083
  fit.germ.shape1 <- fit.germ$estimate[['shape1']]
  fit.germ.shape2 <- fit.germ$estimate[['shape2']]
  
  
  ## calculate p values for candidates
  pvalues <- testBetaBinomial(x, mp, thetahat)
  padj <- p.adjust(pvalues, method="BH")
  x = cbind.data.frame(x, pvalues)
  x = cbind.data.frame(x, padj)
  
  ##calculate likelihood ratio and bayes factor
  lrcut = postcut/mofrac # use LR cutoff corresponding to posterior odds of postcut, given EM mosaic fraction
  x$lr <- mapply(likelihood_ratio, x$altdp, x$N, thetahat, mp)
  
  x$bf <- mapply(bayes_factor, x$altdp, x$N, thetahat, mp, fit.mo.shape1, fit.mo.shape2, n_sample=200, sample_avg_dp) # note: use n_sample=1000 for low coverage
  
  
  ## calculate posterior odds
  x$post_lr <- x$lr * mofrac/(1-mofrac) 
  x$post <- x$bf * mofrac/(1-mofrac)
  
  
  x$col = 'black' # for later plots
  x[x$post>postcut & x$altdp/x$N<=0.5,]$col='red'
  z = x[x$post>postcut & x$altdp/x$N<=0.5,] 
  
  ## Exclude likely mosaic sites and re-estimate parameters
  x.non <- x[x$col=="black",]
  results.non <- estimateTheta(x.non)
  results.non = results.non[results.non$m  > 0, ]
  
  ## theta and p estimations
  if(unname(table(results.non$m>20))[1]<=0.8*nrow(results.non)){ # if more than 20% of DP values in results have >20 variants supporting, use m>20 (only data from DP values with enough support)
    thetahat.non = sum(results.non[results.non$N > 12 & results.non$m > 20, ]$m*results.non[results.non$N > 12 & results.non$m > 20, ]$theta)/sum(results.non[results.non$N > 12 & results.non$m > 20, ]$m)
  } else { # if less than 20% of DP values in results have >20 variants supporting, use m>0 (all data available)
    thetahat.non = sum(results.non$m*results.non$theta)/sum(results.non$m)
  }
  
  print("Running E-M (iteration 2)")
  EMout2 = mofracEM(x, thetahat.non, op, 'yes')
  mofrac = EMout2[1,]$mofrac
  movaf = EMout2[1,]$movaf
  germvaf = EMout2[1,]$germvaf
  mp.non = germvaf
  vmin = EMout2[1,]$minvaf
  vmax = EMout2[1,]$maxvaf
  tmp.mo.vaf = EMout2[EMout2$ind==1,]$allvaf
  tmp.germ.vaf = EMout2[EMout2$ind==2,]$allvaf

  ## fit mosaic distribution to estimate parameters
  #fit.mo <- fitdist(tmp.mo.vaf, 'beta') # shape1=6.655298, shape2=45.098363
  fit.mo <- fit_mosaic(tmp.mo.vaf, sample_avg_dp, outprefix)
  fit.mo.shape1 <- fit.mo[1]
  fit.mo.shape2 <- fit.mo[2]
  
  ## fit germline distribution to estimate parameters
  fit.germ <- fitdist(tmp.germ.vaf, 'beta') # shape1=16.03405, shape2=16.54083
  fit.germ.shape1 <- fit.germ$estimate[['shape1']]
  fit.germ.shape2 <- fit.germ$estimate[['shape2']]
  
  lrcut = postcut/mofrac # use LR cutoff corresponding to posterior odds of 5, given EM mosaic fraction
  mp.non = germvaf
  print(c(mofrac, movaf, germvaf))
  print(c(paste("# mp", mp.non, sep=": "), paste("thetahat", thetahat.non, sep=": ")))
  print(paste("# LR cutoff", lrcut, sep=': '))
  
  ## calculate p-values
  pvalues <- testBetaBinomial(x, mp.non, thetahat.non) 
  padj <- p.adjust(pvalues, method="BH") # add column for adjusted p.value
  
  ##  generate updated set of candidates
  x$pvalues = pvalues
  x$padj = padj
  
  ##  generate updated set of candidates
  ## calculate LR and BF
  x$lr <- mapply(likelihood_ratio, x$altdp, x$N, thetahat.non, mp.non)
  x$bf <- mapply(bayes_factor, x$altdp, x$N, thetahat.non, mp.non, fit.mo.shape1, fit.mo.shape2, n_sample=200, sample_avg_dp) # note: use n_sample=1000 for low coverage
  
  
  ## calculate posterior odds
  x$post_lr <- x$lr * mofrac/(1-mofrac) 
  x$post <- x$bf * mofrac/(1-mofrac)
  
  ## identify candidate mosaics
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
if(length(args)!=5){
  stop("Please provide 5 arguments: <ADfile> <outfile prefix> <posterior odds cutoff> <cohortsize> <sample avg dp>")
}
fname <- toString(args[1])
outprefix <- toString(args[2])
op <- outprefix
postcut <- as.integer(args[3])
cohortsize <- as.integer(args[4])
sample_avg_dp <- as.integer(args[5])
print(c(paste("Input File", fname, sep=': ') , paste("Outfile prefix",outprefix, sep=': '), paste("Posterior Odds Cutoff",postcut, sep=': '), paste("Cohortsize",cohortsize, sep=': '), paste("Dataset Sample Avg DP",sample_avg_dp, sep=': ')))


#############################################
## LOAD DATA AND RUN MAIN FUNCTION
#############################################
a <- read.table(fname, sep='\t', header=T, quote='"')
results = fitReadCounts(a, outprefix)
print("################ DONE !")

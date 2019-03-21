## Usage: Rscript power_analysis.R <denovos file> <output prefix> <LR cutoff> <Thetahat> <Cohort Size> <Sample Avg DP>
## Purpose: For a given pipeline run, generate:
#           1. Curves for estimated Mosaic Detection Power as a function of Sample Avg DP
#           2. Power-adjusted VAF histogram + projected distribution for mosaics VAF>0.1
#           3. Logfile containing a table of counts and adjustment multipliers
library("bbmle")
library("emdbook") 
library("ggplot2")
library("MASS")

#############################################
## HANDLE ARGUMENTS
#############################################
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=6){
  stop("Please provide 6 arguments: <denovos file> <outprefix> <LR cutoff> <theta> <cohort size> <sample avg depth>")
}
fname <- toString(args[1])
op <- toString(args[2])
lrcutoff <- as.integer(args[3])
thetahat <- as.integer(args[4])
cohortsize <- as.integer(args[5])
sample_avg_depth <- as.integer(args[6]) # sample avg dp value to be used in adjusting
print(c(paste("Input File", fname, sep=': ') ,paste("Outfile Prefix",op, sep=': '), paste("LRcut",lrcutoff, sep=': '), paste("Theta",thetahat, sep=': '), paste("Cohort Size",cohortsize, sep=': ')))

#############################################
## LOAD DATA
#############################################
x <- read.table(fname, sep='\t', header=T, quote="")
dn <- x[x$col=="black",]
mo <- x[x$col=="red",]

#############################################
## FUNCTIONS
#############################################
## calculate LR (given Nalt, DP, theta)
llratio <- function(alt, dp, th){
  phat <- alt/dp
  px <- dbetabinom(alt, prob=phat, size=dp, theta=th) # P(D|Mx)
  p1 <- dbetabinom(alt, prob=0.49, size=dp, theta=th) # P(D|M1)
  p0 <- dbetabinom(alt, prob=0, size=dp, theta=th) # P(D|M0)
  x_1 <- px/p1
  x_0 <- px/p0 # irrelevant, since p0 = 0
  return(x_1)
}
# function to find minimum number of alternate read support, given DP and Expected FP 0.01
find_min_alt <- function(N){
  out.min = 6
  for(i in 1:(N/2)){
    exp.fp <- (1-ppois(i, N*(0.005/3)))*3e7
    if(exp.fp <= 0.01){
      out.min = i
      break
    }
  }
  return(out.min)
}
## Find maxBAF for which LR>lcut (given theta and DP) and set minBAF = 6/N
findBAFcuts <- function(N,  lcut, theta) {
  maxBAF = 1e-5 # initial value
  #minBAF = 6/N # min_AD = 6 so minBAF = 6/N
  minBAF = find_min_alt(N)/N
  for (af in c(1:100) / 200 ) { # get maxBAF starting from 1e-5 and moving up
    lr = 0
    ad = round(N * af)
    if(ad>0){
      lr = llratio(ad, N, theta)
    }
    if (lr >= lcut) { # if LR>lcut
      if (af>maxBAF){
        maxBAF = af # update maxBAF
      }
    }
  }
  minBAF = round(minBAF, 3)
  maxBAF = round(maxBAF, 3)
  out = c(minBAF, maxBAF)
  return(out)
}


#############################################
## MOSAIC DETECTION VS VARIANT DP
## Pr(LR>cutoff | DPsite, VAF, LRcutoff) 
#############################################
print("Starting power calculation ...")
dp = seq(1, 1500,by=1)
baf = seq(0.0, 0.5, by=0.01)
LRcut= lrcutoff
theta = thetahat 
prior = 10/LRcut
pwsite = data.frame()
for (i in 1:length(dp)) { ## for each site DP value between 1:1500...
  ## progress message
  if(dp[i] %% 100 == 0){
    print(dp[i])
  }
  bafcuts = findBAFcuts(dp[i], LRcut, theta) # get minBAF and maxBAF for which LR>LRcut given theta
  minBAF = bafcuts[1]
  maxBAF = bafcuts[2]
  min = minBAF
  max = maxBAF
  if(min < max) { 
    for (af in seq(0.0, 0.5, by=0.01)) { # for a given DPsite and theta, for each VAF 0.01:0.50, get the CDF for the mosaic distribution centered at VAF for the range minBAF:maxBAF
      maxcdf = sum(dbetabinom(1:round(max * dp[i]), dp[i], p=af, theta = theta))
      mincdf = sum(dbetabinom(1:round(min * dp[i]), dp[i], p=af, theta = theta))
      cdf = maxcdf - mincdf
      ## adjust for strand bias
      p_no_sb = 1-dbinom(0, round(af*dp[i]), 0.5) # adjust by strand_bias probability
      adj_cdf = cdf*p_no_sb # P(detect mosaic) * P(not strand bias at given dp and af)
      out = c(dp[i], af, min, max, mincdf, maxcdf, cdf, p_no_sb, adj_cdf)
      pwsite = rbind.data.frame(pwsite, out) # store cdf value for fixed LR cutoff and DP
    }
    
  }
  
}
names(pwsite) <- c('dp','af', 'minBAF', 'maxBAF','mincdf', 'maxcdf', 'cdf', 'p_no_sb', 'adj_cdf')

pwsite$altdp <- round(pwsite$dp * pwsite$af, 0)
pwsite$fdr_min <-  mapply(find_min_alt, pwsite$dp)
pwsite$true_min <- mapply(max, 6, pwsite$fdr_min)
table( pwsite$adj_cdf>0 & (pwsite$altdp < pwsite$true_min) )

tmp_pwsite <- pwsite # backup for troubleshooting

pwsite[pwsite$adj_cdf>0 & (pwsite$altdp < pwsite$true_min),]$adj_cdf <- 0.0 ## for conditions with Nalt<max(6, fdr.min.nalt), power should be 0.0

pwsitelogname <- paste(op, 'pwsite.log.txt', sep='.')
write.table(pwsite, pwsitelogname, quote=F, row.names=F, sep='\t')


#############################################
## MOSAIC DETECTION VS SAMPLE AVG DP 
#############################################
sample_power <- function(pw, sdp){
  x <- data.frame()
  baf = seq(0.0, 0.5, by=0.01)
  theta2 = 4 # overdispersion parameter from Sampson 2011
  for(i in 1:length(baf)){ # for each baf value
    tmp <- pwsite[pwsite$af==baf[i],] # get all rows from DPsite dataframe with specific baf value
    m <- nrow(tmp)
    p_dpsite <- dnbinom(tmp$dp, mu=sdp, size=theta2) # for each DPsite value with this particular baf, get the Pr(DPsite | DPsample)
    cdfs <- tmp$adj_cdf * p_dpsite # multiply CDF for each DPsite by Pr(DPsite | DPsample)
    pw <- sum(cdfs) # DPsample power total = sum of adjusted DPsite cdfs
    out <- c(sdp, baf[i], m, pw)
    x <- rbind.data.frame(x, out)
  }
  names(x) <- c('samp_dp', 'af', 'm', 'pw')
  return(x)
}
pwsamp = data.frame()
p_dpsite = data.frame()
for(d in c(40, 60, 80, 100, 120, 150, 200, 300, 400, 500)){
  out <- sample_power(pwsite, d)
  pwsamp <- rbind.data.frame(pwsamp, out)
}
names(pwsamp) <- c('samp_dp', 'af', 'm', 'pw')

pwsamplogname <- paste(op, 'pwsamp.log.txt', sep='.')
write.table(pwsamp, pwsamplogname, quote=F, row.names=F, sep='\t')

title <- 'Mosaic Detection Power vs. Sample Avg. Depth'
p <- ggplot(pwsamp, aes(x=pwsamp$af, y=pwsamp$pw, colour=factor(pwsamp$samp_dp))) + 
  geom_line() +
  xlab(c("VAF")) +
  ylab(c("power")) +
  ggtitle(title) +
  scale_x_continuous(breaks=seq(0.0, 0.5, by=0.05)) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0.0, 1.0, by=0.05)) +
  scale_colour_discrete(name='Sample DP') +
  theme(axis.title = element_text(size=14), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12), legend.justification=c(1,1), legend.position=c(1,1), axis.text = element_text(size=12))+
  guides(colour=guide_legend(ncol=1))+
  theme_bw()
#p

sampdp.fname <- paste(op, 'power_sample_dp.pdf', sep='.')
ggsave(sampdp.fname, p)

#############################################
## POWER ADJUSTED VAF HISTOGRAM   
#############################################
afs <- seq(0.00, 1, by=0.05)
h.mo <- hist(mo[mo$vaf>0.1,]$vaf, breaks=afs, plot=F) ## NEED TO TRY USING LARGER BIN SIZES LIKE 0.05!!!!!!!
cts <- c(0, h.mo$counts)

#names(cts) <- afs
pcgc.sampdp <- sample_avg_depth
pcgc.pwsamp <- pwsamp[pwsamp$samp_dp==pcgc.sampdp,]
cts.pwadj <- rep(0, length(cts))
logdf <- NULL
for(i in 1:length(afs)){
  if(afs[i] > 0.1 & afs[i]<=0.5){   ## ONLY ADJUST BY POWER FOR VAF BINS > 0.1 -- otherwise will hugely inflate low VAF counts
    tmp_pw <- pcgc.pwsamp[pcgc.pwsamp$af==toString(afs[i]), ]$pw
  }
  else{
    tmp_pw <- 1
  }
  cts.pwadj[i] <- cts[i]/tmp_pw
  logdf <- rbind.data.frame(logdf, c(afs[i], tmp_pw, cts[i], cts.pwadj[i]))
}
names(logdf) <-  c('VAF', 'cts', 'cts.pwadj')

cts.pwadj <- cts.pwadj[!is.na(cts.pwadj)]
h.mo$counts <- cts.pwadj
adj.mo.frac <- sum(cts.pwadj)/(nrow(dn)+sum(cts.pwadj))
raw.mo.rate <- sum(cts)/cohortsize
adj.mo.rate <- sum(cts.pwadj)/cohortsize
print(paste('### Raw mosaic count', nrow(mo), sep=' : '))
print(paste('### Raw mosaic fraction', nrow(mo)/nrow(x), sep=' : '))
print(paste('### Raw mosaic rate', raw.mo.rate, sep=' : '))
print(paste('### Adj mosaic count', sum(cts.pwadj), sep=' : '))
print(paste('### Adj mosaic fraction', adj.mo.frac, sep=' : '))
print(paste('### Adj mosaic rate', adj.mo.rate, sep=' : '))
adjvafout <- paste(op, 'vaf_pw_adj.pdf', sep='.')
pdf(adjvafout)
raw.title <- paste('Raw (VAF>10%)', paste(sum(cts), cohortsize, sep='/'), paste(round(raw.mo.rate, 3), 'exome', sep='/') , sep=' = ')
adj.title <- paste('Adjusted (VAF>10%)', paste( round(sum(cts.pwadj),0), cohortsize, sep='/'), paste(round(adj.mo.rate, 3), 'exome', sep='/'), sep=' = ')
title.tf <- paste(raw.title, adj.title, sep='\n')
hist(x$vaf, breaks=seq(0, 1, by=0.05), col=adjustcolor('grey', alpha=0.4), xlab='Variant Allele Fraction (VAF)', main=title.tf, cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
hist(mo$vaf, breaks=seq(0, 1, by=0.05), col=adjustcolor('red', alpha=0.4), add=T)
countdf <- data.frame(afs, cts.pwadj)
lines(countdf[countdf$afs>0.1 & countdf$afs<=0.5,]$afs-0.025, countdf[countdf$afs>0.1 & countdf$afs<=0.5,]$cts.pwadj, col='red')
lines(countdf[countdf$afs>0.1 & countdf$afs<=0.5,]$afs-0.025, countdf[countdf$afs>0.1 & countdf$afs<=0.5,]$cts.pwadj, col='red', type='p')
abline(v=0.1, col='red', lty=2)
legend('topright', c('germline', 'mosaic'), fill=c(adjustcolor('grey', alpha=0.4), adjustcolor('red', alpha=0.4)), cex=1)
dev.off()


## write out log file
logdf <- data.frame(afs, cts, cts.pwadj)
names(logdf) <- c('VAF', 'cts', 'cts.pwadj')
logdf$pw <- c(pcgc.pwsamp[pcgc.pwsamp$af %in% seq(0.00, 1, by=0.05),]$pw, rep(0.0, 12))
logdf <- logdf[logdf$pw>0,]
logname <- paste(op, 'vaf_pw_adj.log.txt', sep='.')
write.table(logdf, logname, quote=F, row.names=F, sep='\t')

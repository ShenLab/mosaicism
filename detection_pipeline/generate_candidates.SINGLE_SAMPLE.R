#############################################
## Script for calculating posterior odds for a given set of variants
## Note: thetahat, mp, mofrac parameters are hard-coded in the main function
##        parameters estimated from cohort of n=2530 WES blood samples (60x)
## Input: ADfile
## Output: annotated ADfile, with columns post=posterior odds score, flag=germline or mosaic
## Usage: Rscript generate_candidates.SINGLE_SAMPLE.R <ADfile> <output file prefix> <posterior odds cutoff>
#############################################

library("bbmle")
library("emdbook") 
library("ggplot2") 


#############################################
## FUNCTIONS
#############################################

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
  
  ## apply FDR-based minimum Nalt filter
  fdr.min.alt <- mapply(find_min_alt, x$N)
  x = cbind.data.frame(x, fdr.min.alt)
  
  #############################################
  ## use parameter values from previous WES studies

  thetahat = 111 ## accounts for overdispersion -- higher = less variance in Nalt, lower = more variance in Nalt
  mp = 0.49 ## mean VAF of germline variants
  mofrac = 0.1205 ## estimated prior mosaic fraction among input denovos
  #############################################

  ##  score variants
  x$lr <- mapply(llratio, x$altdp, x$N, thetahat, mp)
  x$post <- x$lr * mofrac # posterior odds = LR * prior = LR * mosaic fraction
  
  ## define mosaics as any variants with posterior odds > cutoff
  x$flag <- 'germline'
  x[x$post > postcut,]$flag <- 'mosaic'

  ## write out file
  dnname = paste(op, 'denovos.scored.txt', sep='.')

  write.table(x, dnname, quote=F,row.names=F, sep="\t")
  print(paste('####### OUTPUT annotated denovos -- see', dnname, sep=': '))

}

#############################################
## HANDLE ARGUMENTS
#############################################
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3){
  stop("Please provide 3 arguments: <ADfile> <outfile prefix> <posterior odds cutoff>")
}
fname <- toString(args[1])
outprefix <- toString(args[2])
postcut <- as.integer(args[3])

print(c(paste("Input File", fname, sep=': ') , paste("Outfile prefix",outprefix, sep=': '), paste("Posterior Odds Cutoff",postcut, sep=': ')))

#############################################
## LOAD DATA AND RUN MAIN FUNCTION
#############################################
a <- read.table(fname, sep='\t', header=T, quote='"')
results = fitReadCounts(a, outprefix)
print("################ DONE !")

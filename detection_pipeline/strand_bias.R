## Usage: Rscript strand_bias.R <DP4 annotated file> <OUTPUT PREFIX>
## Purpose: Taking DP4 (ADFREF, ADFALT, ADRREF, ADRALT) readcounts, perform Fisher's Exact Test to calculate enrichment and p-value
##          Flag sites as potential false positives if ( (OR<0.33 or OR>3) and p<1e-3 )

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2){
  stop("Please provide 2 arguments: <denovos file> <outfile prefix>")
}
fname <- toString(args[1])
op <- toString(args[2])
print(c(paste("Input File", fname, sep=': '), paste("Outfile prefix", op, sep=': ')))
f <- read.table(fname, sep='\t', header=T)

strand_bias_test <- function(v){
  flag = FALSE
  if(min(v$adfref, v$adfalt, v$adrref, v$adralt) < 1){
    flag = TRUE
  }
  cont <- matrix(c(v$adfref, v$adfalt, v$adrref, v$adralt), nrow=2)
  fish <- fisher.test(cont)
  or <- unname(fish$estimate)
  p <- fish$p.value
  if(p<1e-3){
    if(or<0.33 | or>3){
      flag = TRUE
    }
  }
  return(c(flag, or, p))
}


flags <- c()
ors <- c()
ps <- c()
for(i in 1:nrow(f)){
  if(i %% 1000 == 0){
    print(i)
  }
  tmp <- strand_bias_test(f[i,])
  flag <- tmp[1]
  or <- tmp[2]
  p <- tmp[3]
  flags <- c(flags, flag)
  ors <- c(ors, or)
  ps <- c(ps, p)
}
f$strand_bias_flag <- flags
f$strand_bias_or <- ors
f$strand_bias_p <- ps

print(table(flags))
outname <- paste(op, 'sb.txt', sep='.')
write.table(f, outname, row.names=F, quote=F, sep='\t')
print(paste("DONE see", outname, sep=':'))
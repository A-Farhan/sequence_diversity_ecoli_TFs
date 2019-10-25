## generate random TRNs to get the probability of getting observed div diff b/w TFs & TGs
# using diversity estimates from above runs

## require
source('rm_tf_trn.R')
source('colblindpal.R')

## files
frn = 'trn.tab'
fdiv = 'snsdiv_asm.tab'
fp = 'ran_trn_div.svg'

## diversity
divs <- read.table(file = fdiv,header = T,sep = '\t',stringsAsFactors = F)
divs <- divs[divs$runs >= 100, ]

## trn
trn <- read.table(file = frn,header = F,sep = '\t',stringsAsFactors = F)
trn <- split( trn[,2], trn[,1])
# without tfs
trn <- rm_tf_trn(trn)
# tfs
fs <- names(trn)
# no. of tfs
nf = length(fs)
# targets
gs <- unique( unlist( sapply( trn, function(x) strsplit(x,',')), recursive = T))
gs <- gs[ !gs %in% fs]
# only keep targets for which div is available
gs <- gs[ gs %in% divs$gene]

# TF diversity
fd <- divs$nonsyn[match(fs,divs$gene)]
# mean diversity of correspondin tus
gd <- sapply( trn_wo_f, function(x) 
  mean( sapply( strsplit(x,','), function(y) 
    mean( divs$nonsyn[ match(y,divs$gene)], na.rm = T)), na.rm = T ))

# get statistic of Wilcoxon signed rank test
ob <- wilcox.test(x = fd, y = gd, alternative = 'l', paired = T)
vo <- as.numeric(ob$statistic)

## function to generate random trn and get div test statistic
random_trn_div <- function(trn,targets,divs,targets_div){
  # number of tfs
  nf <- length(trn)
  # pick transcription factors
  fs <- sample(targets,nf)
  ## generate random trn
  rtrn <- trn
  names(rtrn) <- fs
  # get diversity for random tfs
  fdiv <- divs$nonsyn[match(fs,divs$gene)]
  # get random statistic
  vr <- as.numeric(wilcox.test(x = fdiv,y = targets_div,alternative = 'l',paired = T)$statistic)
  return(vr)
}

drb <- replicate(n = 1000,
      expr = random_trn_div(trn = trn_wo_f,targets = gs,divs = divs,targets_div = gd))

## plotting
svg(fp,width = 5,height = 7)
par(mar=c(5,4.5,2,2))
X = sort(drb)
L = length(X)
Y = sapply(X, function(x) sum(X<=x)/L)
# random TRN curve
plot( X, Y, type='l',cex.lab = 1.5, cex.axis = 1.5, xlab = 'W statistic, Wilcoxon test', ylab = 'Cumulative Frequency',
     lwd=2, col = colorblind[3])
# real TRN line
segments(x0 = vo,y0 = 0,x1 = vo,y1 = 1,col = colorblind[2],lwd=2)
segments( x0 = min(drb)-100, y0 = 1, x1 =  max(drb)+100, y1 = 1, lty=2)
segments( x0 = min(drb)-100, y0 = 0, x1 =  max(drb)+100, y1 = 0, lty=2)

# significance threshold
segments(x0 = min(drb)-100,y0 = 0.05,x1 = max(drb)+100,y1 = 0.05,col = colorblind[1],lwd=1)
mtext(text = 0.05,side = 2,line = 0.5,at = 0.05,adj = 0)
legend(x = 3500,y = 0.2,legend = c('real TRN','random TRN'),lty=1,lwd=2,col = c(colorblind[2],colorblind[3]), bty='n', cex = 1.2, x.intersp = 0.2)
dev.off()

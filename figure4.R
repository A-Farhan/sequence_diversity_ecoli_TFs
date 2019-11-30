## conservation of E coli TFs across species

## require 
source('pvalue_stars.R')
source('colblindpal.R')
source('sci_notation.R')
library(aqfig)
library(RColorBrewer)

## arguments
# minimum no. of TUs for a TF to be GR
minu = 10

## files
f1 = 'mean_std_nsdiv_wgs.tab' ## mean standardised Non-Synonymous Diversity
f2 = 'gene_orthocount.tab' # number of orthologs per gene
frn = 'trn.tab' # TRN
fp1 = 'figure4.svg' # plot output

## TRN
trn <- read.table(file = frn,header = F,sep = '\t',stringsAsFactors = F)
trn <- split( trn[,2], trn[,1])
# tfs
fs <- names(trn)
# no. of tfs
nf = length(fs)
# identify non global regulators
ngr <- fs[ sapply( trn, function(x) length(x) < minu)]
# no. of non-GRs
nr = length(ngr)
# targets
gs <- unique( unlist( sapply( trn, function(x) strsplit(x,',')), recursive = T))
# no. of targets
ng <- length(gs)
# number of TUs
ntus = sapply(trn,length)

## load diversity
d1 = read.table(file = f1,header = T,sep = '\t',stringsAsFactors = F)
# remove NA rows
d1 <- d1[ !is.na(d1$div), , drop=F]
# TFs'
fdiv = d1[fs,'div']
names(fdiv) = fs
# means of mean diversity of TUs of TFs
gdiv = sapply( trn, function(x) 
  mean( sapply( strsplit( x,',',fixed=T),function(y) 
    mean( d1[ y,'div'], na.rm = T)), na.rm = T))

## load conservation
d2 = read.table(file = f2,header = T,sep = '\t',
                stringsAsFactors = F, row.names = 1)
# no. of species
N = max(d2$selected)
# column for conservation fraction
d2$fc = d2$selected/N

# conservation and diversity for all genes
d1 <- d1[ rownames(d1) %in% rownames(d2), , drop=F]
X <- d1$div
Y <- d2[ rownames(d1), 'fc']
umat <- cbind(X,Y)
rownames(umat) <- rownames(d1)
umat  <- umat[ order(umat[ , 1]), ]

# conservation fraction for different classes
cf = d2[ fs, 'fc']
cg = d2[ gs, 'fc']
cr = d2[ ngr, 'fc']
cuf = sapply( fs, function(x) mean( 
  sapply( strsplit(trn[[x]],','), function(y) 
    mean( d2[ y, 'fc'], na.rm = T)), na.rm=T) )
cur = cuf[ngr]

### plotting
svg(fp1)
par(mar = c(4,4.5,2,2))
par(mfrow=c(2,2))
cols <- colorblind

# plot1: Conservation against diversity for all genes
nb = 30 # no. of bins
X = umat[, 1] - min( umat[, 1]) + 1
x <- log(X)
y <- umat[,2]

scatterplot.density(x, y,num.bins = nb,col.regression.line = NULL,col = kristen.colors(12), 
                    col.one.to.one.line = NULL, xlab='Diversity within species', 
                    ylab='Conservation across species', main='',xaxt='n',
                    cex.lab = 1.2, cex.axis = 1.2, density.in.percent = F)
# x-axis ( Diversity)
axis(side = 1, at = seq(min(x),max(x),length.out = 5), 
     labels = round( seq(min(umat[,1]),max(umat[,1]),length.out = 5)), 
     cex.axis = 1.2)
# test
ob <- cor.test( umat[,1], umat[,2], alternative = 'l',method = 's')
d <- round(ob$estimate,2)
pv <- scinot(ob$p.value)
Lines <- bquote( paste( rho, ' = ',.(d),', ', 
                        italic(P)['spearman correlation'], ' = ', .(pv[1])^.(pv[2]) ))
mtext(text = Lines, side = 3, adj = 0, cex = 0.9)
# sample size
mtext( text = bquote( paste( italic('N'), ' = ', .(length(x)))), 
       side = 4, line = 1.8, adj = 1, cex = 0.9)

## plot2: conservation of TFs and TGs, paired (ALL TFs)
plot( x = cuf, y = cf, xlim = c(0,1), ylim = c(0,1), pch = 20, col = cols[3],
      xlab = 'TG Conservation', ylab = 'TF conservation',
      cex.axis = 1.2, cex.lab = 1.2)
segments( x0 = 0, y0 = 0, x1 = 1, y1 = 1)
# sample size
mtext( text = bquote( paste( italic('N'), ' = ', .(length(cf)))), 
       side = 4, line = 0.5, adj = 1, cex = 0.9)
# test
ob <- wilcox.test( x = cf, y = cuf, alternative = 'l', paired = T)
pv <- round(ob$p.value,3)
Lines <- bquote( paste( italic(P)['wilcoxon signed rank'], ' = ', .(pv)) )
mtext(text = Lines, side = 3, adj = 1, cex = 0.9)

## plot 3: TF conservation and nTUs
mat = cbind( ntus, cf)
plot( mat, pch=20, log='x', xlab='Number of regulated TUs', ylab = 'TF Conservation',
      cex.lab=1.2, cex.axis=1.2, xaxt = 'n', col = cols[3])
# x-axis
axis(side = 1, at = c(1,2,5,10,50,200), cex.axis = 1.2 )
# sample size
mtext( text = bquote( paste( italic('N'), ' = ', .(length(cf)))), 
       side = 4, line = 0.5, adj = 1, cex = 0.9)
# test
ob <- cor.test( x = mat[,1], y = mat[,2], 
                alternative = 'g', method = 's')
d <- round( ob$estimate, 2)
pv <- scinot(ob$p.value)
Lines <- bquote( paste( rho, ' = ',.(d),', ', 
                        italic(P)['spearman correlation'], ' = ', .(pv[1])^.(pv[2]) ))
mtext(text = Lines, side = 3, adj = 0, cex = 0.9)

## plot 4: conservation of TFs and TGs, paired (non-GR TFs)
plot( x = cur, y = cr, xlim = c(0,1), ylim = c(0,1), pch = 20, col = cols[3],
      xlab = 'TG Conservation', ylab = 'TF conservation',
      cex.axis = 1.2, cex.lab = 1.2)
segments( x0 = 0, y0 = 0, x1 = 1, y1 = 1)
# sample size
mtext( text = bquote( paste( italic('N'), ' = ', .(length(cr)))), 
       side = 4, line = 0.5, adj = 1, cex = 0.9)
# test
ob <- wilcox.test( x = cr, y = cur, alternative = 'l', paired = T)
pv <- round(ob$p.value,3)
Lines <- bquote( paste( italic(P)['wilcoxon signed rank'], ' = ', .(pv)) )
mtext(text = Lines, side = 3, adj = 1, cex = 0.9)

dev.off()
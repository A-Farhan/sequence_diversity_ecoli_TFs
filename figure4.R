## conservation of E coli TFs across species

## require 
source('pvalue_stars.R')
source('colblindpal.R')
library(aqfig)
library(RColorBrewer)
library(dichromat)

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
gs <- gs[ !gs %in% fs]
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
umat  <- umat[ order(umat[ , 1]), ]
rownames(umat) <- rownames(d1)

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

scatterplot.density(x, y, num.bins = nb, col.regression.line = NULL, col = kristen.colors(12), 
                    col.one.to.one.line = NULL, xlab='Diversity within species', 
                    ylab='Conservation across species', main='',xaxt='n',
                    cex.lab = 1.2, cex.axis = 1.2, density.in.percent = F)
# x-axis
axis(side = 1, at = seq(min(x),max(x),length.out = 5), 
     labels = round( seq(min(umat[,1]),max(umat[,1]),length.out = 5)), 
     cex.axis = 1.2)
# mention sample size
mtext( text = paste('N = ',length(x)), side = 4, line = 1.5, adj = 0, cex = 0.8)
# test
ob <- suppressWarnings( cor.test( umat[,1], umat[,2], alternative = 'l',method = 's'))
pv <- pv_stars(ob$p.value)
d <- round(ob$estimate,2)
Lines <- bquote( paste( 'Spearman correlation, ', 
                        rho, ' = ',.(d),.(pv_stars(pv))))
mtext(text = Lines, side = 3, adj = 1, cex = 0.9)

## plot2: TF conservation and nTUs
mat = cbind( ntus, cf)
plot( mat, pch=20, log='x', xlab='Number of regulated TUs', ylab = 'TF Conservation', col = cols[3],
      cex.lab=1.2, cex.axis=1.2, xaxt = 'n')
# x-axis
axis(side = 1, at = c(1,2,5,10,50,200), cex.axis = 1.2 )

# test
ob <- suppressWarnings( cor.test( x = mat[,1], y = mat[,2], 
                alternative = 'g', method = 's'))
pv <- pv_stars(ob$p.value)
d <- round( ob$estimate, 2)
Lines <- bquote( paste( 'Spearman correlation, ', 
                        rho, ' = ',.(d),.(pv_stars(pv))))
mtext(text = Lines, side = 3, adj = 1, cex = 0.9)

## plot3: conservation of TF, TG and non-GR
# rescale paired differences
df = cf-cuf
dr = cr-cur
mn = min(df,na.rm = T)
mx = max(df,na.rm = T)
ef = (df - mn) / (mx-mn)
er = (dr - mn) / (mx-mn)

boxplot( cf,cg,cr,ef,er,notch=T, col='grey', cex.lab = 1.2, border = 'grey', 
         cex.axis = 0.8, axes = F, pch = 20, cex = 0.6, boxwex = 0.5 )
box()
# x-axis 
# with ticks but without labels
axis(1, labels = F)
# labels
labels = c('TF','TG','LR','TF','LR')
# plot x labs at default x position
text(x =  seq_along(labels), y = par("usr")[3]-0.1, srt = 0,
     labels = labels, xpd = T, cex = 1.2)

# first y-axis
axis(side = 2, at = seq(0,1,0.2), cex.axis=1.2)
mtext(text = 'Conservation', side = 2, line = 3)
# second y-axis
axis(side = 4,at = seq(0,1,0.2), labels = round(seq(mn,mx,0.2),1),
     cex.axis=1.2)
mtext(text = 'Conservation [ TF-TG ]', side = 4, line = 2.2)
# sample sizes at the top
mtext( text = c(nf,ng,nr), side = 3, at = 1:3, cex = 1)
# vertical line dividing two sides
segments(x0 = 3.5,y0 = 0,x1 = 3.5,y1 = 1)
# line marking TG median for left side plot
md = median(cg,na.rm = T)
segments(x0 = 0,y0 = md,x1 = 3.5,y1 = md,lty = 1)
# line marking null for the right side plot
x = -mn/(mx-mn)
segments(x0 = 3.5,y0 = x,x1 = 5.5,y1 = x,lty = 1)
# test
ob <- suppressWarnings( wilcox.test( x = cf, y = cg, alternative = 'l'))
pv <- round(ob$p.value,3)
d <- round( median(cf,na.rm = T) - median(cg,na.rm = T), 2)

## plot4: TF TG conservation distribution
plot( density( cf, na.rm = T), col = cols['orange'], lwd=2,main='', xaxt = 'n', 
      xlab = 'Conservation', ylab = '', cex.lab=1.2, cex.axis=1.2)
lines( density( cg,na.rm = T), col = cols['blue'], lwd=2)
# x-axis
axis( side = 1, at = seq(0,1,0.1), labels = F)
mtext( text = c(0,0.5,1), side = 1, line = 1, at = c(0,0.5,1))
# vertical lines for second peaks
segments(x0 = 0.7, y0 = 0, x1 = 0.7, y1 = 0.6)
segments(x0 = 0.91, y0 = 0, x1 = 0.91, y1 = 0.55)
legend('topright',legend = c('TF','TG'), col = c(cols['orange'],cols['blue']), lty = 1, lwd = 2, cex = 1.2)
mtext(text = 'Density',side = 4,line = 0.5)
dev.off()

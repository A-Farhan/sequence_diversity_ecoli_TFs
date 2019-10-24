## plots for diversity-gene length correlation

# require
source('pvalue_stars.R')
source('colblindpal.R')
source('sci_notation.R')

## files
fdiv1 = 'nucdiv_asm.tab'
fdiv2 = 'snsdiv_asm.tab'
frn = 'trn.tab'
fp = 'figure1.svg'

## arguments
span = 0.5 # span for loess smoothing
trans = 0.4 # transparency in colors

## trn
trn <- read.table(file = frn,header = F,sep = '\t',stringsAsFactors = F)
trn <- split( trn[,2], trn[,1])
# tfs
fs <- names(trn)
# no. of tfs
nf = length(fs)
# targets
gs <- unique( unlist( sapply( trn, function(x) strsplit(x,',')), recursive = T))
# no. of targets
ng <- length(gs)

## load data on diversity
div = read.table(file = fdiv1, header = T, sep = '\t', 
                 stringsAsFactors = F, row.names = 1)

div2 = read.table(file = fdiv2,header = T,sep = '\t',stringsAsFactors = F)
# find indices of genes of first data in the second data
ixs = match( rownames(div), div2$gene)
div$syn = div2$syn[ixs]
div$nonsyn = div2$nonsyn[ixs]
# exclude all rows with NA
div <- div[ !apply(div,1,anyNA), ]
# sort by no. of sites
div <- div[ order(div$nuc_nogap), ]

### diversity differences b/w TFs and TGs
## nuc
nucf = div[ fs, 'nuc_div']
nucg = div[ gs, 'nuc_div']
nucu = sapply( trn, function(x) mean( 
  sapply( strsplit(x,',',fixed=T), function(y) mean( div[ y,'nuc_div'], na.rm = T)), na.rm = T))
# unpaired
nucwt = wilcox.test( x = nucf, y = nucg, alternative = 'l')
# paired
pnucwt = wilcox.test( x = nucf, y = nucu, alternative = 'l', paired = T)

## syn
synf = div[ fs, 'syn']
syng = div[ gs, 'syn']
synu = sapply( trn, function(x) mean( 
  sapply( strsplit(x,',',fixed=T), function(y) mean( div[ y,'syn'], na.rm = T)), na.rm = T))
# unpaired
synwt = wilcox.test( x = synf, y = syng, alternative = 'l')
# paired
psynwt = wilcox.test( x = synf, y = synu, alternative = 'l', paired = T)

## non-syn
nonsynf = div[ fs, 'nonsyn']
nonsyng = div[ gs, 'nonsyn']
nonsynu = sapply( trn, function(x) mean( 
  sapply( strsplit(x,',',fixed=T), function(y) mean( div[ y,'nonsyn'], na.rm = T)), na.rm = T))
# unpaired
nonsynwt = wilcox.test( x = nonsynf, y = nonsyng, alternative = 'l')
# paired
pnonsynwt = wilcox.test( x = nonsynf, y = nonsynu, alternative = 'l', paired = T)

### spearman correlation b/w paired div diff (TF-TG) and TF length
# do length based sorting first
mat = cbind( div[ fs, 'nuc_len'], nucf-nucu, synf-synu, nonsynf-nonsynu )
mat = mat[ order(mat[,1]), ]
# scale div values
mat[,2] = ( mat[,2] - min(mat[,2]) ) / ( max(mat[,2]) - min(mat[,2]) )
mat[,3] = ( mat[,3] - min(mat[,3]) ) / ( max(mat[,3]) - min(mat[,3]) )
mat[,4] = ( mat[,4] - min(mat[,4]) ) / ( max(mat[,4]) - min(mat[,4]) )

# nuc
nucsc <- cor.test( mat[,1], mat[,2], alternative = 'g', method = 'spearman')
# syn
synsc <- cor.test( mat[,1], mat[,3], alternative = 'g', method = 'spearman')
# non-syn
nonsynsc <- cor.test( mat[,1], mat[,4], alternative = 'g', method = 'spearman')

### linear regression
# nuc
nuclr <- lm( mat[,2]~mat[,1])
snuclr <- summary(nuclr)
# syn
synlr <- lm( mat[,3]~mat[,1])
ssynlr <- summary(synlr)
# non-syn
nonsynlr <- lm( mat[,4]~mat[,1])
snonsynlr <- summary(nonsynlr)

## function to print P-values
fun <- function(num){ 
  sci = scinot(num)
  out = bquote( paste( italic('P = '), .(sci[1])^.(sci[2]) ))
}

### plotting
svg(fp)
# color-blind friendly colors
cols1 = c( rgb(0,0.6,0.5,trans), rgb(0.9,0.6,0,trans), rgb(0.8,0.6,0.7,trans))
cols2 = colorblind[c(4,2,8)]

par(mar = c(4,4.1,2,2), oma = c(3,0,0,0), bty='n')
layout( matrix( c(1,1,1,2,2,2,3,3,3,4,4,5,5,6,6,7,7,7), nrow = 2, ncol = 9, byrow = T) )

## plot 1-3: boxplot for unpaired div diff b/w TF and TG

# plot 1: nucdiv
boxplot( nucf, nucg, col = cols1[1], notch = T, outline = F, pch = 20, names=c('TF','TG'), border = cols2[1], 
         ylab = '', cex.axis = 1.2, cex.lab = 1.2, boxwex = 0.5, yaxt = 'n', ylim = c(0, 0.04) )
axis( side = 2, at = seq(0,0.04,0.01), cex.axis = 1.2)
segments( x0 = 0.5, y0 = median(nucg), x1 = 2.5, y1 = median(nucg))
mtext( text = fun(nucwt$p.value), side = 3, line = -1, cex = 0.95)

# sample sizes
mtext( text = c(nf,ng), side = 1, at = 1:2, line = 2.7)

# plot 2: syndiv
boxplot( synf, syng, col = cols1[2], notch = T,outline = F, pch = 20, names=c('TF','TG'), border = cols2[2], 
         ylab = '', cex.axis = 1.2, cex.lab = 1.2, boxwex = 0.5, yaxt = 'n', ylim = c(0, 0.14) )
axis( side = 2, at = round( seq(0,0.14,length.out = 5), 3), cex.axis = 1.2)
segments( x0 = 0.5, y0 = median(syng), x1 = 2.5, y1 = median(syng))
mtext( text = fun(synwt$p.value), side = 3, line = -1, cex = 0.95)

# plot 3: nonsyndiv
boxplot( nonsynf, nonsyng, col = cols1[3], notch = T,outline = F, pch = 20, names=c('TF','TG'), border = cols2[3], 
         ylab = '', cex.axis = 1.2, cex.lab = 1.2, boxwex = 0.5)
segments( x0 = 0.5, y0 = median(nonsyng), x1 = 2.5, y1 = median(nonsyng))
mtext( text = fun(nonsynwt$p.value), side = 3, line = -1, cex = 0.95)

## plot 4-6: boxplot for paired div diff b/w TF and TG
par( mar = c(2,3,1,0))
# plot 4: nucdiv
boxplot( nucf-nucu, col = cols1[1], notch = T,outline = F, pch = 20, border = cols2[1], 
         ylab = '', cex.axis = 1.2, cex.lab = 1.2, boxwex = 0.5, ylim = c(-0.02,0.02))
mtext( text = 'TF - TG', side = 1, line = 0)
segments( x0 = 0.5, y0 = 0, x1 = 1.5, y1 = 0)
mtext( text = fun(pnucwt$p.value), side = 3, line = -1, cex = 0.95)

# plot 5: syndiv
boxplot( synf-synu, col = cols1[2], notch = T,outline = F, pch = 20, border = cols2[2], 
         ylab = '', cex.axis = 1.2, cex.lab = 1.2, boxwex = 0.5, ylim = c(-0.07,0.07), yaxt = 'n')
axis( side = 2, at = round( seq( -0.07, 0.07, length.out = 5), 3), cex.axis = 1.2)
mtext( text = 'TF - TG', side = 1, line = 0)
segments( x0 = 0.5, y0 = 0, x1 = 1.5, y1 = 0)
mtext( text = fun(psynwt$p.value), side = 3, line = -1, cex = 0.95)

# plot 6: nonsyndiv
boxplot( nonsynf-nonsynu, col = cols1[3], notch = T,outline = F, pch = 20, border = cols2[3], 
         ylab = '', cex.axis = 1.2, cex.lab = 1.2, boxwex = 0.5)
mtext( text = 'TF - TG', side = 1, line = 0)
segments( x0 = 0.5, y0 = 0, x1 = 1.5, y1 = 0)
mtext( text = fun(pnonsynwt$p.value), side = 3, line = -1, cex = 0.95)

## plot 7: correlation b/w gene length and nucleotide diversity #############
par( mar = c(4,2,1,1), bty = 'o')
# scatter
plot( mat[,1], mat[,2], pch = 20, col = cols1[1], cex.lab = 1.5, cex.axis = 1.2, xlim = c(200,2000), 
      ylim = c(0.2,0.7), xlab='Gene length', ylab = '')
points( mat[,1], mat[,3], pch = 20, col = cols1[2])
points( mat[,1], mat[,4], pch = 20, col = cols1[3])
# y-lab
# mtext( text = 'TF - TG scaled diversity', side = 2, line = 2.5)

## loess curves
Z1 <- predict( loess( mat[,2] ~ mat[,1], span=span))
lines( mat[,1], Z1, lwd=2,col=cols2[1])
Z2 <- predict( loess( mat[,3] ~ mat[,1], span=span))
lines( mat[,1], Z2, lwd=2,col=cols2[2])
Z3 <- predict( loess( mat[,4] ~ mat[,1], span=span))
lines( mat[,1], Z3, lwd=2,col=cols2[3])

# legend
labels = expression( 'nucleotide', paste(italic('S'),'ynonymous'), paste(italic('N'),'on-synonymous'))
legend(x = -4700, y = 0.1, legend = labels, cex = 1.5,
       fill = cols2, bty = 'n', xpd = NA, horiz = T, title = 'Diversity', title.adj = 0, 
       text.width = 1000, x.intersp = 0.5, border = NA)

# print correlation values
s1 = round( nucsc$estimate, 2)
s2 = round( synsc$estimate, 2)
s3 = round( nonsynsc$estimate, 2)

Lines <- bquote( paste( rho, ' = ',.(s1),'*, ',
                        rho[italic('S')], ' = ',.(s2),'*, ', 
                        rho[italic('N')], ' = ',.(s3) ) )
mtext( text = Lines, side = 3, adj = 1, cex = 0.95, line = 0.8)

# describe P
ln = expression( paste( '[ ', italic('P'), ': Wilcoxon test ]'))
mtext( text = ln, side = 1, line = 0, outer = T, at = 0.3)

# describe asterisk
ln = bquote( paste( rho, ': Spearman correlation', ', *: ', italic(P), '<', 0.05))
mtext( text = ln, side = 1, cex = 0.9, line = 1, outer = T, adj = 1)

dev.off()

# require
source('pvalue_stars.R')
source('colblindpal.R')
source('rm_tf_trn.R')

# files
fsn = '_snsdiv.tab'
frn = 'trn.tab'
frun = '_runs.txt'
fsp = 'specific_tfs.txt'
fp = 'figure3.svg'

# all div files
fls = Sys.glob(paste0('*',fsn))
# correspondin projects
prjs = gsub(pattern = fsn, replacement = '', fls)
# no. of projects
N = length(prjs)

# no. of TU threshold
t = 10
# minimum no. of runs
mr = 50

## TRN
trn = read.table(file = frn, header = F, sep = '\t', stringsAsFactors = F)
trn = split( trn[,2], trn[,1])
trn = rm_tf_trn(trn)
fs = names(trn)
nf = length(fs)
gs = unique( unlist( sapply( trn, function(x) strsplit(x,','))))
gs = gs[!gs %in% fs]
ng = length(gs)
# gene of interest
goi <- c(fs,gs)
ngoi <- length(goi)

## number of TUs per TF
ntu = sapply(trn,length)
# specific TFs
spf = readLines(fsp)
# general TFs
genf <- c('lrp','argR','crp','fruR','fnr','arcA','narL','fis','ihfA','purR')
genf <- genf[ genf %in% fs]

## initialize matrix to hold diversity values across projects
pdiv = matrix( data = NA, nrow = ngoi, ncol = N)

# for every project
for( i in 1:N){
  # load data
  dsn = read.table( paste0( prjs[i], fsn ), 
                    header = T, sep = '\t', stringsAsFactors = F, row.names = 1)
  # no. of runs
  nr = length( readLines( paste0( prjs[i], frun)) )
  
  # skip the project, if not sufficient runs
  if(nr < mr) next
  
  # scale data with median & MAD of TGs
  nsg = dsn[ gs, 'nonsyn']
  mg = median( nsg, na.rm = T)
  madg =  mad(x = nsg, constant = 1, na.rm = T)
  dsn$nonsyn = (dsn$nonsyn - mg)/madg
  pdiv[,i] = dsn[ goi, 'nonsyn']

}

# remove empty columns
pdiv = pdiv[ , !apply( pdiv, 2, function(x) all(is.na(x)))]
rownames(pdiv) <- goi
# take means across projects
mdiv = rowMeans( pdiv, na.rm = T)

## plotting
svg(fp)
layout( mat = matrix( c(1,2), nrow = 1, ncol = 2, byrow = F), 
        widths = c(1.8,1))
par( mar = c(4.1,2,1.2,2), oma = c(0,2.2,0,0))
colors <- colorblind

## plot1: TF Div and nTUs correlation
fd = mdiv[fs]
mat = cbind(ntu,fd)
mn = min( fd, na.rm = T)
mx = max( fd, na.rm = T)

plot( mat, log = 'x', ylim = c( mn, mx), cex.lab = 1.2, cex.axis = 1.2, pch = 19, 
      col = colors['sky_blue'], xaxt = 'n',
      xlab='Number of regulated TUs', ylab = '')
# x-axis
axis( side = 1, at = c(1,2,5,10,50,200), cex.axis = 1.2)
mtext(text = paste('N = ', nrow(mat), sep = ' '), side = 3, adj = 0)
# y-axis label
mtext( text = 'TF Diversity', side = 2, line = 1, outer = T, cex = 1.2)
# test
ob1 <- suppressWarnings( cor.test( x = mat[,1], y = mat[,2], alternative = 'l', method = 'spearman'))
pv <- ob1$p.value
cat( 'Correlation b/w TF diversity and its regulon size, P = ', pv, '\n')
ss <- pv_stars(pv)
d <- round( ob1$estimate, 3)
Lines <- bquote( paste( 'Spearman correlation, ', rho, ' = ',.(d),.(ss)))
mtext(text = Lines, side = 4, adj = 1, line = 0.4, cex = 1)

## plot2: boxplot of div of general and specific TFs
mdgf = mdiv[genf]
mdsf = mdiv[spf]

boxplot( mdgf, mdsf, outline = F, col='grey', xaxt = 'n', 
         ylab = '', cex.axis = 1.2, cex.lab = 1.2)
# x-axis
# with ticks but without labels
axis( side = 1, at = 1:2, labels = F)
labels = c('General','Specific')
# Plot x labs at default x position
text( x =  1:2, y = par("usr")[3] - 0.1, srt = 45, adj = 1,
     labels = labels, xpd = T)
# sample sizes on top
mtext(text = c( length(mdgf), length(mdsf)), side = 3, at = 1:2)
# test
ob2 <- wilcox.test( x = mdgf,y = mdsf, alternative = 'l')
pv <- ob2$p.value
cat( 'Diversity General < Specific, P = ', pv, '\n')
ss <- pv_stars(pv)
d <- round( median( mdgf, na.rm = T) - median( mdsf, na.rm = T), 2)
# math expression
ex = bquote(paste('M'['G']-'M'['S'],' = ', .(d), .(ss), ' Wilcoxon test' ))
mtext( text = ex, side = 4, adj = 1, cex = 1, line = 0.4)
# describe asterisks
ln = bquote( paste( '** :', italic(P), '<', 0.01))
mtext(text = ln, side = 4, line = 0.4, adj = 0, cex = 1)
dev.off()
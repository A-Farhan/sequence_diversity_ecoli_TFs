# require
source('pvalue_stars.R')

# files
fin1 = '_snsdiv.tab'
fin2 = '_runs.txt'
frn = 'trn.tab'
fp = 'figure2.svg'

# all nucdiv files
fls = Sys.glob(paste0('*',fin1))
# correspondin projects
prjs = gsub(pattern = fin1, replacement = '', fls)
# no. of projects
N = length(prjs)

# minimum no. of runs
mr = 50

## TRN
trn = read.table(file = frn,header = F,sep = '\t', stringsAsFactors = F)
trn = split(trn[,2],trn[,1])
# tfs
fs = names(trn)
# targets
gs = unique( unlist( sapply( trn, function(x) strsplit( x, ',', fixed=T))))

## p vector
pv <- NULL

# vector of number of runs
nruns <- rep(NA,N)
names(nruns) <- prjs

## initialize list to hold TFs standardised diversity distributions
out = list()

for( i in 1:N){
  # no. of runs
  nr = length( readLines( paste0( prjs[i], fin2)) )
  # capture no. of runs in the vector
  nruns[i] <- nr
  # skip the project, if not sufficient runs
  if(nr < mr) next

  # load data
  dsn = read.table( paste0( prjs[i], fin1 ), 
                    header = T, sep = '\t', stringsAsFactors = F, row.names = 1)
  
  # get Non-syn div of target genes and its median & MAD
  nsg = dsn[ gs, 'nonsyn']
  mg = median( nsg, na.rm = T)
  madg =  mad(x = nsg, constant = 1, na.rm = T)
  # scale whole data with above median and MAD
  dsn$nonsyn = (dsn$nonsyn - mg)/madg
  
  # non-syn div of tfs
  X = dsn[ fs, 'nonsyn']
  # means of mean diversity of TUs of TFs
  Y = sapply( trn, function(x) mean( 
    sapply( strsplit(x,',',fixed = T), function(y) mean(dsn[y,'nonsyn'], na.rm = T)), na.rm = T))
  
  # save the difference
  out[[ prjs[i] ]] = X-Y
  # test
  ob <- wilcox.test( x = X, y = Y, alternative = 'l', paired = T)
  # extend p vector
  pv <- c( pv, ob$p.value)
}

## adjust p-values for multiple testing
padj <- p.adjust( p = pv, method = 'holm')

## plotting
svg(fp, width = 8)
par( mar = c(5.5,4.2,1.5,2), bty = 'o')

boxplot( x = out, notch = T, outline = T, col = 'grey', pch = 20, cex = 0.5, 
         ylab = 'Non-synonymous diversity TF - TG', ylim = c(-5,5), border = 'grey40',
         cex.lab = 1.2, cex.axis = 1.2, names = F)

# line marking expectation
segments(x0 = 0, y0 = 0, x1 = length(out)+0.5, y1 = 0, lty = 1)

# significance symbol
ss <- sapply( padj, pv_stars)
ss[ss == 'ns'] = ''
mtext(text = ss, side = 3, at = 1:length(out), cex = 1.2, line = 0, col = 'grey20')

# labels
labels <- names(out)
text(x =  seq_along(labels), y = par("usr")[3] - 0.5, srt = 40, adj = 1, cex = 1.1,
     labels = labels, xpd = T)

# describe asterisks
ln = bquote( paste('** :', italic(P), '<', 0.01))
mtext(text = ln, side = 4, line = 0.5, adj = 0)

# sample size
ln = bquote( paste( italic(N), ' = ', .(length(fs))) )
mtext(text = ln, side = 4, line = 0.5, adj = 1)

dev.off()

## fraction of data within y-limits
t = unlist(out)
nt = length(t)
fna = sum(is.na(t))/nt
t = t[!is.na(t)]
f = sum(t >= -5 & t <= 5)/length(t)

cat('Fraction of NA: ',fna,'\n')
cat('Fraction of data in the plot: ',f,'\n')
cat('P-values were in the range: ', min(padj), '-', max(padj), '\n')
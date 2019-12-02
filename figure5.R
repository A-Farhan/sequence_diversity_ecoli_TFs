## test if mutations were more frequent in TFs compared to TGs in lab evolution experiments

# require
source('colblindpal.R')

## significance threshold
mp = 0.05

## plot file
fp = 'figure5.svg'

# file suffix
fsx = '_mfreq.tab'

# all files with above suffix
fls = Sys.glob( paste0('*',fsx) )

# types of projects
tps = unique( sapply( strsplit(x = fls, split = '_', fixed = T),'[',1))

## initialize a list to hold p-value vectors for each project type
out = list()

for(k in 1:length(tps)){
  fnames = grep( pattern = paste0(tps[k],'_'), x = fls, value = T)
  N = length(fnames)
  pvs = rep(NA,N)
  names(pvs) <- gsub(fsx,'',fnames)
  for(i in 1:N){
    f = fnames[i]
    mat = read.table(file = f,header = T,sep = '\t', row.names = 1, stringsAsFactors = F)
    pvs[i] = fisher.test(mat,alternative = 'g')$p.value
  }
  out[[k]] = pvs
}
names(out) = tps

## log transform
prjp = unlist(out)
prjp = round(-log10(prjp),3)

## filter data
sel = grep( pattern = 'le|ma|50000', x = names(prjp), value=T)
prjp = prjp[sel]
N = length(prjp)

## plotting
svg(fp)
par(mar=c(4.5,5,2,2))

# color for different categories
cols <- rep( NA, N)
for(i in 1:N){
  is_le <- grepl( pattern = 'le_P', x = names(prjp)[i], ignore.case = F)
  is_ab <- grepl( pattern = 'le_ab', x = names(prjp)[i], ignore.case = F)
  is_ltee <- grepl( pattern = 'ltee', x = names(prjp)[i], ignore.case = F)
  if(is_le)
    cols[i] <- colorblind[2]
  else if(is_ab) 
    cols[i] <- colorblind[3]
  else if(is_ltee)
    cols[i] <- colorblind[5]
  else
    cols[i] <- colorblind[8]
}

barplot( height = prjp, space = 0, xlab = 'Projects', 
         ylab = expression(paste('-log'[10],'P')), xaxt = 'n', cex.lab = 1.5, 
         col = cols )
# threshold line
cut = -log10(mp)
segments(x0 = -10,y0 = cut, x1 = 100, y1 = cut)
# legend
legend('topright', legend = c('LE','AR','LTEE','MA'), fill = colorblind[c(2,3,5,8)], bty = 'n', cex = 1.2)
dev.off()
## overall test b/w adaptive lab evolution and MA
ma = prjp[ grep(pattern = 'ma', names(prjp))]
ale = prjp[ !prjp %in% ma]
ob <- wilcox.test( x = ale, y = ma, alternative = 'g', exact = F)
cat('Adaptive lab evolution experiments show TF mutation enrichment more often than MA studies, P = ', ob$p.value,'\n')

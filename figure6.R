## make a panel of 2x3 figures, rows for non-mutator & mutators
# 3 columns for DAF, mutated site fraction, diversity

source('colblindpal.R')
source('rm_tf_trn.R')
source('pvalue_stars.R')

## arguments
# minimum no. of TUs for TF to be GR
mu = 10
# no. of classes of proteins
nclass = 3
# generation interval size
isize = 10000
# significance threshold
t1 = -log10(0.05)
# color transparency
trans = 0.3

# files
frn <- 'trn.tab'
fft <- 'reference.ft'
plotf <- 'figure6.svg'

## TRN
trn = read.table(file = frn, header = F, sep = '\t', stringsAsFactors = F)
trn = split( trn[,2], trn[,1])
trn = rm_tf_trn(trn)
# all TFs
tf = names(trn)
ntf = length(tf)
# targets
tg = unique( unlist( sapply( trn, function(x) strsplit(x,','))))
tg = tg[!tg %in% tf]
ntg = length(tg)
# GRs
gr = names(trn)[ sapply(trn, function(x) length(x) >= mu)]
ngr = length(gr)
# "local" regulators
lr = tf[ !tf %in% gr]
nlr = length(lr)

# reference feature table
rft <- read.table(file = fft, header = F, sep = '\t', quote = "", stringsAsFactors = F)
rft <- rft[ rft$V1 == 'CDS' & rft$V12 != '', ]
gene_len <- rft$V18
names(gene_len) <- rft$V15

# no. of sites of target, local regulators and GRs
nsites <- c( sum(gene_len[tg]), sum(gene_len[lr]), sum(gene_len[gr]))

# left and right end of generation intervals & total no.
lastg = 60000
glow = seq( 0, lastg-isize, isize)
ghigh = seq( isize, lastg, isize)
nint = length(ghigh)

# suffix for input files
sfx = '_ts.tab'
# populations
pops = c('m5','m6','p1','p2','p4','p5','p6','p3','m2','m4')
np = length(pops)
# non mutators
npops = pops[1:6]
# mutators
mpops = pops[7:np]

# get time points common to all populations
gens <- as.integer( Reduce( f = intersect, 
                            x = sapply( pops, function(x) 
                              unlist( strsplit( x = readLines( paste0(x,sfx))[1], split = '\t'))[-(1:3)])))
ngen = length(gens)
# corresponding column names
xgen <- paste0('X',gens) 

## function to get no. of new mutations 
# given a binary matrix of mutations presence over time intervals
# 1 row per mutation
# 1 col per time interval
n_dnm_inv <- function(bmat,nint){
  
  # list of rows indices corresponding to de novo mutations in the interval
  X = list()
  for(j in 1:nint){
    # for the first interval, simply keep all mutations that are present
    if(j==1)
      X[[j]] = which(bmat[,j]==1)
    # for others, exclude the indices already captured for previous intervals
    else{
      t = which(bmat[,j]==1)
      X[[j]] = t[!t %in% unlist(X)]
    }
  }
  # no. of de novo mutations for each interval 
  nm = sapply( X, length)
  return(nm)
}

## function to give max freq of each mutation over a time period, 
# given a freq matrix and list of columns indices corresponding to time points, 
maxf_inv <- function( fmat, col_inxs){
  # find maximum freq of each mutation over the time interval
  maxf <- sapply( col_inxs, function(X)
    apply( fmat[ , X], 1, max))
  if( length(dim(maxf)) != 0)
  # exclude zero freq mutations
  maxf <- apply( maxf, 2, function(x) x[ x != 0])
  return(maxf)
}

## function that takes as input a file of time series and returns
# 1. DAF 2. de novo mutation frequency 3. freq distribution over time intervals
fun <- function(infile,tp){
  # load data
  data = read.table(file = infile, header = T, sep = '\t', stringsAsFactors = F)
  # only keep mutations that are at zero freq at first time point
  data = data[ data$X0 == 0, ]
  # limit data to selected time points
  data = data[ , c( colnames(data)[1:3], tp)]
  
  # separate frequencies and mutation info
  fdata = data[ , -(1:3)]
  info = data[ , 1:3]
  
  # column names corresponding to generations
  gcols = colnames(fdata)
  # corresponding generations
  gens = as.integer( gsub(pattern = 'X', replacement = '', x = gcols))
 
  ## part 1: DAF
  
  # no. of time points
  ngen = length(gens)
  # split data by gene name
  gdata <- split( fdata, data$gene)
  # get gene names
  gnames = names(gdata)
  # no. of genes
  ngenes = length(gnames)
  # get daf for all genes
  daf <- t( sapply( gdata, function(X) apply( X, 2, function(x) sum(x))))

  # initialize separate matrices for TG, LR and GR
  mtg <-  matrix( 0, ntg, ngen)
  rownames(mtg) <- tg
  mlr <-  matrix( 0, nlr, ngen)
  rownames(mlr) <- lr
  mgr <-  matrix( 0, ngr, ngen)
  rownames(mgr) <- gr
  # fill above matrices
  tmp <- daf[ rownames(daf) %in% tg,]
  mtg[ rownames(tmp), ] <- tmp
  tmp <- daf[ rownames(daf) %in% lr,]
  mlr[ rownames(tmp), ] <- tmp
  tmp <- daf[ rownames(daf) %in% gr,]
  mgr[ rownames(tmp), ] <- tmp

  # take average per time point
  avtg <- colMeans(mtg, na.rm = T)
  avlr <- colMeans(mlr, na.rm = T)
  avgr <- colMeans(mgr, na.rm = T)
  out1 = rbind( avtg,avlr,avgr)
  colnames(out1) <- gens
  
  ## part 2. De novo mutation frequency
  # change class info to identify GRs and non-GRs
  x = info$class
  x[info$gene %in% gr] = 2
  info$class = x
  
  # for every generation interval, get vector of indices for columns corresponding to the interval
  inxs = lapply( 1:nint, function(x) which( gens > glow[x] & gens <= ghigh[x]))
  # for the list of indices above, sum frequencies over the indices
  sumf = sapply( inxs, function(x) rowSums( fdata[, x]) )
  # convert the above to a binary matrix
  # such that 1 represent non-zero freq of a mutation over an interval
  fbin = ifelse(test = sumf != 0, yes = 1, no = 0)
  # attach info back
  fbin = cbind( info, fbin)  
  
  # split the above by class
  cdata = split(x = fbin[,-(1:3)], f = fbin$class)
  # apply function to the above list
  out2 = sapply( cdata, n_dnm_inv, nint = nint)
  
  ## part 3. Frequency distribution of pooled mutations
  # split frequency data by class
  cdata = split( x = fdata, f = info$class)
  # apply function to both matrices and save in the pooling list
  out3 = lapply( cdata, maxf_inv, col_inxs = inxs)
  
  out = list(daf=out1,mut=out2,frq=out3)
  return(out)
}

# apply above function to given populations
rs <- lapply( pops, function(x) 
  fun( infile = paste0(x,sfx), tp = xgen) )
names(rs) <- pops

### extraction and transformations

## DAF 
daf_ls = list()
for(i in 1:nclass){
  # non-mutator matrix, average and SD
  daf_mat_n = sapply( rs[npops], function(X) X[[1]][i,])
  daf_avg_n = rowMeans(daf_mat_n)
  daf_sd_n = apply( daf_mat_n, 1, sd)
  # mutator matrix, average and SD
  daf_mat_m = sapply( rs[mpops], function(X) X[[1]][i,])
  daf_avg_m = rowMeans(daf_mat_m)
  daf_sd_m = apply( daf_mat_m, 1, sd)
  avgs <- cbind( daf_avg_n, daf_avg_m)
  sds <- cbind( daf_sd_n, daf_sd_m)
  daf_ls[[i]] <- list( avg = avgs, sd = sds, 
                       matn = daf_mat_n, matm = daf_mat_m)
}

## DMF de novo mutation frequency
dmf_ls = list()
for(i in 1:nclass){
  # non-mutator
  dmf_mat_n = sapply( rs[npops], function(X) X[[2]][,i])/nsites[i]
  # mutator
  dmf_mat_m = sapply( rs[mpops], function(X) X[[2]][,i])/nsites[i]
  dmf_ls[[i]] <- list(matn = dmf_mat_n, matm = dmf_mat_m)
}

## DPF distribution of pooled frequencies 
dpf_ls = list()
for(i in 1:nclass){
  # non-mutator
  tmp = lapply( rs[npops], function(X) X[[3]][[i]])
  dpf_n = sapply( 1:nint, function(j) 
    unlist( sapply( tmp, function(X) X[[j]])) )
  # mutator
  tmp = lapply( rs[mpops], function(X) X[[3]][[i]])
  dpf_m = sapply( 1:nint, function(j) 
    unlist( sapply( tmp, function(X) X[[j]])) )
  dpf_ls[[i]] <- list(lsn = dpf_n,lsm = dpf_m)
}

### tests

## DAF p-value time series 

# pvalue vector for nonmutators and mutators
pvn <- matrix( NA, nclass, ngen )
pvm <- pvn
for(i in 1:ngen){
  # non-mutators
  drbtg <- daf_ls[[1]]$matn[i,]
  drblr <- daf_ls[[2]]$matn[i,]
  drbgr <- daf_ls[[3]]$matn[i,]
  # test b/w TG and LR
  pvn[1,i] <- -log10( suppressWarnings( wilcox.test( x = drbtg, y = drblr, alternative = 'l', paired = T)$p.value))
  # test b/w LR and GR
  pvn[2,i] <- -log10( suppressWarnings( wilcox.test( x = drblr, y = drbgr, alternative = 'l', paired = T)$p.value))
  # test b/w GR and TG
  pvn[3,i] <- -log10( suppressWarnings( wilcox.test( x = drbgr, y = drbtg, alternative = 'g', paired = T)$p.value))
  
  # mutators
  drbtg <- daf_ls[[1]]$matm[i,]
  drblr <- daf_ls[[2]]$matm[i,]
  drbgr <- daf_ls[[3]]$matm[i,]
  # test b/w TG and LR
  pvm[1,i] <- -log10( suppressWarnings( wilcox.test( x = drbtg, y = drblr, alternative = 'l', paired = T)$p.value))
  # test b/w LR and GR
  pvm[2,i] <- -log10( suppressWarnings( wilcox.test( x = drblr, y = drbgr, alternative = 'l', paired = T)$p.value))
  # test b/w GR and TG
  pvm[3,i] <- -log10( suppressWarnings( wilcox.test( x = drbgr, y = drbtg, alternative = 'g', paired = T)$p.value))
}
# loess smoothing
smpvn <- t( apply( pvn, 1, function(x) predict( loess(x ~ gens, span=0.2)) ))
smpvm <- t( apply( pvm, 1, function(x) predict( loess(x ~ gens, span=0.2)) ))
daf_test <- list(n = smpvn, m = smpvm)

## DMF

# p-value matrices
pvm_n <- matrix( data = NA, nrow = nint, ncol = choose(nclass,2), 
               dimnames = list( ghigh, c('TGvLR','LRvGR','GRvTG')) )
pvm_m <- pvm_n
for(i in 1:nint){
  ## non-mutator
  drbtg = dmf_ls[[1]]$matn[i,]
  drblr = dmf_ls[[2]]$matn[i,]
  drbgr = dmf_ls[[3]]$matn[i,]
  # TG vs TF
  pvm_n[i,1] <- suppressWarnings( 
    wilcox.test( x = drbtg, y = drblr, paired = T, alternative = 'l')$p.value)
  # TF vs GR
  pvm_n[i,2] <- suppressWarnings( 
    wilcox.test( x = drblr, y = drbgr, paired = T, alternative = 'l')$p.value)
  # GR vs TG
  pvm_n[i,3] <- suppressWarnings( 
    wilcox.test( x = drbgr, y = drbtg, paired = T, alternative = 'g')$p.value)

  ## mutator
  drbtg = dmf_ls[[1]]$matm[i,]
  drblr = dmf_ls[[2]]$matm[i,]
  drbgr = dmf_ls[[3]]$matm[i,]
  # TG vs TF
  pvm_m[i,1] <- suppressWarnings( 
    wilcox.test( x = drbtg, y = drblr, paired = T, alternative = 'l')$p.value)
  # TF vs GR
  pvm_m[i,2] <- suppressWarnings( 
    wilcox.test( x = drblr, y = drbgr, paired = T, alternative = 'l')$p.value)
  # GR vs TG
  pvm_m[i,3] <- suppressWarnings( 
    wilcox.test( x = drbgr, y = drbtg, paired = T, alternative = 'g')$p.value)
}
# convert both matrices to -log10P
ss_n <- -log10(pvm_n)
ss_m <- -log10(pvm_m)
dmf_test = list(n = ss_n, m = ss_m)

## DPF
# p-value matrices
pvm_n <- matrix( data = NA, nrow = nint, ncol = choose(nclass,2), 
               dimnames = list( ghigh, c('TGvTF','TFvGR','GRvTG')) )
pvm_m <- pvm_n
for(i in 1:nint){
  ## non-mutators
  drbtg <- dpf_ls[[1]]$lsn[[i]]
  drblr <- dpf_ls[[2]]$lsn[[i]]
  drbgr <- dpf_ls[[3]]$lsn[[i]]
  # TG vs TF
  pvm_n[i,1] <- suppressWarnings( 
    wilcox.test( x = drbtg, y = drblr, alternative = 'l')$p.value)
  # TF vs GR
  pvm_n[i,2] <- suppressWarnings( 
    wilcox.test( x = drblr, y = drbgr, alternative = 'l')$p.value)
  # GR vs TG
  pvm_n[i,3] <- suppressWarnings( 
    wilcox.test( x = drbgr, y = drbtg, alternative = 'g')$p.value)

  ## mutators
  drbtg <- dpf_ls[[1]]$lsm[[i]]
  drblr <- dpf_ls[[2]]$lsm[[i]]
  drbgr <- dpf_ls[[3]]$lsm[[i]]
  # TG vs TF
  pvm_m[i,1] <- suppressWarnings( 
    wilcox.test( x = drbtg, y = drblr, alternative = 'l')$p.value)
  # TF vs GR
  pvm_m[i,2] <- suppressWarnings( 
    wilcox.test( x = drblr, y = drbgr, alternative = 'l')$p.value)
  # GR vs TG
  pvm_m[i,3] <- suppressWarnings( 
    wilcox.test( x = drbgr, y = drbtg, alternative = 'g')$p.value)
}
# convert both matrices to -log10P
ss_n <- -log10(pvm_n)
ss_m <- -log10(pvm_m)
dpf_test = list(n = ss_n, m = ss_m)

## plotting
svg(plotf)
par( mfrow = c(3,2))
# add extra space to right margin of plot within frame
par( mar=c(0, 3, 0, 0.5), oma = c(4,4,3,4 ))
# colors
cols = colorblind[c(4,2,8,3)]
fcols = c( rgb(0,0.6,0.5,trans), rgb(0.9,0.6,0,trans), rgb(0.8,0.6,0.7,trans))
# plot titles
pltl = c('Derived allele frequency', 
         'De novo mutation frequency', 
         'Variant frequency\ndistribution')

## plot 1: DAF for non-mutators
# GR
# means
Y = daf_ls[[3]]$avg[,1]
# error
er = daf_ls[[3]]$sd[,1]
# curve
plot( gens, Y, col=cols[3], pch=20, type='l', lwd=2, ylim=c(0,0.06), 
      xlab='', axes=F, ylab = '')
# error bars
#segments(x0 = gens, y0 = Y-er, x1 = gens, y1 = Y+er, 
#           col = cols[3], lty = 3)
# LR
# means
Y = daf_ls[[2]]$avg[,1]
# error
er = daf_ls[[2]]$sd[,1]
# curve
lines( gens, Y, lwd=2, col=cols[2], pch=20)
# error bars
#segments(x0 = gens, y0 = Y-er, x1 = gens, y1 = Y+er, 
#         col = cols[2], lty = 3)
# TG
# means
Y = daf_ls[[1]]$avg[,1]
# error
er = daf_ls[[1]]$sd[,1]
# curve
lines( gens, Y, lwd=2, col=cols[1], pch=20)
# error bars
#segments(x0 = gens, y0 = Y-er, x1 = gens, y1 = Y+er, 
#         col = cols[1], lty = 3)

# y-axis 1
axis( 2, at = seq( 0, 0.06, 0.02), 
      labels = seq( 0, 6, 2), col="black", las=1) 
box()

# add p-value curve
par(new=T)
plot( daf_test[[1]][1, ], x=gens, type='l', lty = 3, col=cols[1], ylim = c(0,2), 
        axes=F, xlab='', ylab='')
lines( daf_test[[1]][2, ], x=gens, col=cols[2], lty = 3) 
lines( daf_test[[1]][3, ], x=gens, col=cols[3], lty = 3) 
# significance threshold line
segments(x0 = 0, y0 = t1, x1 = 60000, y1 = t1, col = cols[4])
# Add vertical grid
abline(v = seq( 0, max(gens), 5000), col = "lightgrey" )
#segments(x0 = 0,y0 = -100,x1 = 0,y1 = 2, xpd=T, outer=T)

## plot 2: DAF for mutators

# GR
# means
Y = daf_ls[[3]]$avg[,2]
# error
er = daf_ls[[3]]$sd[,2]
# curve
plot( gens, Y, col=cols[3], pch=20, type='l', lwd=2, ylim=c(0,0.2), 
      xlab='', ylab='', axes=F)
# error bars
#segments(x0 = gens, y0 = Y-er, x1 = gens, y1 = Y+er, 
#         col = cols[3], lty = 3)
# LR
# means
Y = daf_ls[[2]]$avg[,2]
# error
er = daf_ls[[2]]$sd[,2]
# curve
lines( gens, Y, lwd=2, col=cols[2], pch=20)
# error bars
#segments(x0 = gens, y0 = Y-er, x1 = gens, y1 = Y+er, 
#         col = cols[2], lty = 3)
# y-axis 1
axis( 2, at = seq(0, 0.2, 0.05), 
      labels = seq(0, 20, 5), col="black", las=1) 

# LR
# means
Y = daf_ls[[1]]$avg[,2]
# error
er = daf_ls[[1]]$sd[,2]
# curve
lines( gens, Y, lwd=2, col=cols[1], pch=20)
# error bars
#segments(x0 = gens, y0 = Y-er, x1 = gens, y1 = Y+er, 
#         col = cols[1], lty = 3)
# y-axis 1
axis( 2, at = seq(0, 0.2, 0.05), 
      labels = seq(0, 20, 5), col="black", las=1) 

# add p-value curve
par(new=T)
plot( daf_test[[2]][1, ], x=gens, type='l', lty = 3, col=cols[1], ylim = c(0,2), 
      axes=F, xlab='', ylab='')
lines( daf_test[[2]][2, ], x=gens, col=cols[2], lty = 3) 
lines( daf_test[[2]][3, ], x=gens, col=cols[3], lty = 3) 

# significance threshold line
segments(x0 = 0, y0 = t1, x1 = 60000, y1 = t1, col = cols[4])
# Add vertical grid
abline(v = seq( 0, max(gens), 5000), col = "lightgrey")
# second y-axis (significance)
axis( side = 4, ylim = c(0,2) )
box()

## plot 3: DMF for non-mutators
# transparency
trans = 0.3
# color-blind friendly colors
# fill colors
fcols = c( rgb(0,0.6,0.5,trans), rgb(0.9,0.6,0,trans), rgb(0.8,0.6,0.7,trans))
# border colors
bcols = colorblind[c(4,2,8)]

# TG
boxplot( t(dmf_ls[[1]]$matn), ylim = c( 0, 0.0002), xlim = c(0,24), 
         at = seq(1,23,4),  
         col = fcols[1], border = cols[1], outline = F, 
         main='', axes=F, ylab = '')
# non-GR
boxplot( t(dmf_ls[[2]]$matn), col = fcols[2], border = cols[2], outline = F, 
         add = T, at = seq(2,23,4), 
         axes = F)
# GR
boxplot( t(dmf_ls[[3]]$matn), col = fcols[3], border = cols[3], outline = F, 
         add = T, at = seq(3,23,4), axes = F)
# y-axis
axis( side = 2, at = seq( 0, 0.00016, by = 0.00004),
      labels = seq( 0, 16, 4))
# add a line for zero
segments(x0 = 0, y0 = 0, x1 = 24, y1 = 0, lwd = 0.5)
# add p-value curve
par(new=T)
X = sort( c(seq(1,23,4), seq(2,23,4), seq(3,23,4)))
Y = as.numeric(t(dmf_test$n))
plot( x=X, y=Y, pch = 20, cex = 0.8, col=cols[4], xlim = c(0,24), ylim = c(0,2), 
      axes=F, xlab='', ylab='')
# significance threshold line
segments(x0 = 0, y0 = t1, x1 = 24, y1 = t1, col = cols[4])
box()

## plot 4: DMF for mutators
# TG
boxplot( t(dmf_ls[[1]]$matm), ylim = c( 0, 0.0012), xlim = c(0,24), 
         at = seq(1,23,4),  
         col = fcols[1], border = cols[1], outline = F, 
         main='', axes=F )
# non-GR
boxplot( t(dmf_ls[[2]]$matm), col = fcols[2], border = cols[2], 
         outline = F, add = T, at = seq(2,23,4), 
         axes = F)
# GR
boxplot( t(dmf_ls[[3]]$matm), col = fcols[3], border = cols[3], 
         outline = F, add = T, at = seq(3,23,4), 
         axes = F)
# y-axis
axis( side = 2, at = seq( 0, 0.001, by = 0.00025),
      labels = seq( 0, 100, 25))
# add a line for zero
segments(x0 = 0, y0 = 0, x1 = 24, y1 = 0, lwd = 0.5)
# add p-value curve
par(new=T)
X = sort( c(seq(1,23,4), seq(2,23,4), seq(3,23,4)))
Y = as.numeric(t(dmf_test$m))
plot( x=X, y=Y, pch = 20, cex = 0.8, col=cols[4], xlim = c(0,24), ylim = c(0,2), 
      axes=F, xlab='', ylab='')
# significance threshold line
segments(x0 = 0, y0 = t1, x1 = 24, y1 = t1, col = cols[4])
# second y-axis (significance)
axis( side = 4, ylim = c(0,2) )
box()

## plot 5: DPF for non-mutators
# TG
boxplot( dpf_ls[[1]]$lsn, ylim = c(0,1), xlim = c(0,24), 
         at = seq(1,23,4), 
         col = fcols[1], border = cols[1], outline = F, 
         main='', axes=F,
         ylab = pltl[3], cex.lab = 1.2)
# non-GR 
boxplot( dpf_ls[[2]]$lsn, col = fcols[2], border = cols[2], 
         outline = F, add = T, at = seq(2,23,4), 
         axes = F)
# GR
boxplot( dpf_ls[[3]]$lsn, col = fcols[3], border = cols[3], 
         outline = F, add = T, at = seq(3,23,4), 
         axes = F)
# x-axis
axis( side = 1, at = seq(3,23,4), labels = ghigh)
# y-axis
axis( side = 2, at = seq( 0, 1, by = 0.2))
# add p-value curve
par(new=T)
X = sort( c(seq(1,23,4), seq(2,23,4), seq(3,23,4)))
Y = as.numeric(t(dpf_test$n))
plot( x=X, y=Y, pch = 20, cex = 0.8, col=cols[4], xlim = c(0,24), ylim = c(0,2), 
      axes=F, xlab='', ylab='')
# significance threshold line
segments(x0 = 0, y0 = t1, x1 = 24, y1 = t1, col = cols[4])
box()

## plot 6: DPF for mutators
# TG
boxplot( dpf_ls[[1]]$lsm, ylim = c(0,1), xlim = c(0,24), at = seq(1,23,4), 
         col = fcols[1], border = cols[1], outline = F, 
         main='', axes=F )
# non-GR 
boxplot( dpf_ls[[2]]$lsm, col = fcols[2], border = cols[2], outline = F, add = T, 
         at = seq(2,23,4), axes = F)
# GR
boxplot( dpf_ls[[3]]$lsm, col = fcols[3], border = cols[3], outline = F, add = T, 
         at = seq(3,23,4), axes = F)
# x-axis
axis( side = 1, at = seq(3,23,4), labels = ghigh)
# y-axis
axis( side = 2, at = seq( 0, 1, by = 0.2))
# add p-value curve
par(new=T)
X = sort( c(seq(1,23,4), seq(2,23,4), seq(3,23,4)))
Y = as.numeric(t(dpf_test$m))
plot( x=X, y=Y, pch = 20, cex = 0.8, col=cols[4], xlim = c(0,24), ylim = c(0,2), 
      axes=F, xlab='', ylab='')
# significance threshold line
segments(x0 = 0, y0 = t1, x1 = 24, y1 = t1, col = cols[4])
# second y-axis (significance)
axis( side = 4, ylim = c(0,2) )
box()

# x-axis label
mtext("Generations",side=1,col="black",line=2.5,outer = T)
# y-axis labels
mtext( text = pltl[3], side = 2, line = 3, outer = T, at = .15, 
       padj = 1)

mtext( text = pltl[2], side = 2, line = 2, outer = T, at = .5)
mtext( text = bquote( paste( '(x',10^-5,')')),
       side = 2, line = 0, outer = T, at = .5)

mtext( text = pltl[1], side = 2, line = 2, outer = T, at = .84)
mtext( text = bquote( paste( '(x',10^-2,')')),
       side = 2, line = 0, outer = T, at = .84)

mtext( expression( paste('-log'[10],'P')), side=4, line=2.5, outer=T)
# legend
legend( x = -32, y = 6.8, ncol = 7, bty = 'n',
        legend = c('TG','LR','GR', 'TGvLR', 'LRvGR', 'GRvTG', expression(italic('P'))), 
        lty = rep(c(1,3), c(3,3)), col = c( rep( cols[1:3], 2), cols[4]), 
        xpd = NA, seg.len = 0.8, cex = 1.3, lwd = 2)
dev.off()

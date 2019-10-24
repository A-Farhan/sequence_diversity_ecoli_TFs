# infer genes presence/absence based on depth of coverage for each gene region

# command line arguments
args <- commandArgs(T)
dr <- args[1]

# paths
fcb <- file.path( dr, 'cov_breadth.bed')
fcd <- file.path( dr, 'cov_depth.bed')
fgp <- file.path( dr, 'genes_present.txt')
fmg <- file.path( dr, 'missing_genes.txt')
fex <- file.path( dr, 'excess_depth_genes.txt')

### breadth of coverage ###
dcb <- read.table( file = fcb, header = F, sep = '\t', stringsAsFactors = F)
# genes with zero breadth of coverage
gcb0 <- dcb$V4[ dcb$V10 == 0]

# minimum coverage threshold,
mbt = 0.6 # fraction of gene length overlapping at least one read
gpb <- dcb$V4[ dcb$V10 >= mbt]

### depth of coverage ###
dcd <- read.table( file = fcd, header = F, sep = '\t', stringsAsFactors = F)
# add column for per base coverage
dcd$V8 <- dcd$V7/(dcd$V3 - dcd$V2)
# genes with zero depth of coverage
gcd0 <- dcd$V4[ dcd$V8 == 0]

# find the depth cut-off to subset dcd
# first extract non-zero values
dp <- dcd$V8[ dcd$V8 > 0]
# find their mean
md <- mean(dp)
# further, keep only values smaller than twice of the above mean
dp <- dp[ dp < 2*md]
# recalculate mean
md <- mean(dp) 
# get standard deviation
std <- sd(dp)
# minimum depth threshold, for a gaussian distribution 3SD covers 99.7% data
mdt = md - 3*std
# extract genes satsifying above threshold
gpd <- dcd$V4[ dcd$V8 >= mdt & dcd$V8 > 0]
# genes with excess depth of coverage
E <- dcd$V4[ dcd$V8 >= 2*md ]

# take intersection of genes selected for breadth & depth of coverage
gp <- intersect(gpb,gpd)
# take union of genes with either zero breadth or zero depth
g0 <- union( gcb0, gcd0)

# write above vector of genes to files
write( x = gp, file = fgp)
write( x = g0, file = fmg)
write( x = E, file = fex)

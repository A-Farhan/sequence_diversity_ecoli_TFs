## Function to get the trn excluding tfs
rm_tf_trn <- function(trn){
  # tfs
  fs <- names(trn)
  # no. of tfs
  nf = length(fs)

  trn_wo_f <- list()
  for(i in 1:nf){
    tus = trn[[i]]
    ntu = length(tus)
    # initialize vector to hold modified tus
    vec <- NULL
    for(j in 1:ntu){
      # break tu into individual genes
      genes <- unlist( strsplit(tus[j],','))
      # exclude tfs
      genes <- genes[!genes %in% fs]
      # rejoin the genes
      genes <- paste(genes,collapse=',')
      vec <- c(vec,genes)
    }
    # remove empty elements from vector
    vec = vec[vec != ""]
    # if there is at least one tu, then make an entry
    if(length(vec)!=0) trn_wo_f[[ fs[i] ]] <- vec
  }
  return(trn_wo_f)
}
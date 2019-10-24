# function to get stars for p values
pv_stars <- function(pv){ 
  if(is.na(pv)) return('')
  else{ 
    if(pv > 0.05) s = 'ns'
    else if(pv >= 0.01 & pv <= 0.05) s = '*'
    else s = '**'
    return(s)
  }
}
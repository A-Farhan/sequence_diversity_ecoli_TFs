scinot <- function(num,s=1){
  # scientific notation E form
  sci = formatC( num, format = "e", digits = s)
  # separate components
  parts = unlist( strsplit( x = sci, split = 'e', fixed = T))
  # decimal part
  d = parts[1]
  # exponent part
  p = as.numeric( parts[2])
  # generate x10 form and return components
  sci = paste0( d, ' x 10')
  out = c( sci, p)
  return(out)
}
# Function not implemented in R. This implementation is identical to the formula in Matlab with flag=0, to give consistent results
skewness <- function(x) {  
  m3 <- mean((x-mean(x))^3)
  skew <- m3/ ( sqrt( ( mean((x-mean(x))^2) ) ) ^3)
  return(skew)
}

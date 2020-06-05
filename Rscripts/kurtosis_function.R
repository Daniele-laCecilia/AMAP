# Function not implemented in R. This implementation is identical to the formula in Matlab with flag=0, to give consistent results
kurtosis <- function(x) {  
  m4 <- mean((x-mean(x))^4)
  kurt <- m4/ ( ( mean((x-mean(x))^2) ) ^2)
  return(kurt)
}

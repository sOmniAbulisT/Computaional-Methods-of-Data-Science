rm(list = ls())

m <- 50000; p <- 20
sigma <- sqrt((1:p)/2)
rhalfnorm <- function(n, s){
  s * abs(rnorm(n))
}

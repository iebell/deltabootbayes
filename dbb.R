## Problem 1

# Bernoulli sampler function
bSample <- function (rep, n, p) {
  x <- replicate(rep, rbinom(n, 1, p), simplify = "matrix")
  return(x)
}

#Delta method calculations
ods <- function(p){
  as.numeric((p/(1-p)))
}

delta <- function(n, d, p) {
  sqrt(n)*(ods(d)-ods(p)) 
}

#Problem 1 data and results

sampleSize <- c(10, 30, 50, 100, 500, 1000)
#Each sample of size N is one column of the arrays in p1list
p1List <- lapply(sampleSize, function (x) {bSample(5000, x, 0.5)})
#
p1Means <- lapply(p1List, function (y) {apply(y, 2, mean)})
p1Delta <- p1Means
for (i in 1:6) {
  p1Delta[[i]] <- vapply(p1Means[[i]], 
    function (x) {delta(sampleSize[i], x, 0.5)}, numeric(1))
}
hist(p1Delta[[6]], breaks = 100)

##Problem 2, Delta Method
s2 <- c(10, 30, 100, 1000)
p2 <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
sims <- expand.grid(
  n = s2,
  p = p2
)
p2List <- apply(sims, 1, function(x) {bSample(10, x[1], x[2])})
p2Means <- lapply(p2list, function (y) {apply(y, 2, mean)})
p2Exp <- p2Means
for (r in 1:length(sims[,1])) {
  
}
p2Var <-

##Problem 2, Bootstrap
  
  
  
  
  
  
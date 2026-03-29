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
  p1Delta[[i]] <- sapply(p1Delta[[i]], 
    function (x) {delta(sampleSize[i], x, 0.5)})
}

##Problem 2, Delta Method
s2 <- c(10, 30, 100, 1000)
p2 <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
sims <- expand.grid(
  n = s2,
  p = p2
)
# Samples generated, do not re-reun
# p2List <- apply(sims, 1, function(x) {bSample(100, x[1], x[2])})
saveRDS(p2List, file = "p2list.rds")
p2ListPerm <- readRDS("p2list.rds")

#Filter the lists to get 

p2ListFiltered <-

#Obtain the expectation and variances
expectation <- function(x) {
  mean(x)/(1-mean(x))
}
deltaVar <- function (y) {
  (1/(1-mean(y)))^4 * var(y)
}
p2Exp <- lapply(p2ListPerm, function (y) {apply(y, 2, expectation)})
p2Var <- lapply(p2ListPerm, function (y) {apply(y, 2, deltaVar)})


##Problem 2, Bootstrap


bootSamples <- function (x, b) {
    replicate(b, sample(x, length(x), replace = TRUE), simplify = TRUE)
}


bootstrapMeans <- function(x, b) {
  
  samples <- bootSamples(x, b)
  
#If any of the samples are all 1's, our estimates will be useless.
#We replace any all one's column with a column of n-1 1's and one 0.
  
  replace(bootSamples, which(colSums(samples) == nrow(samples), 
                             c(r, ))
  
  bMean <- function (y) {
    y <- ifelse(is.infinite(expectation(y)), 10000, expectation(y))
  }
  
  # Calculate the mean
  bSampleMeans <- mean((apply(samples, 2, bMean)))
  return(bSampleMeans)
  
}

  
  
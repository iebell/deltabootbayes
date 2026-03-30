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
#p2List <- apply(sims, 1, function(x) {bSample(100, x[1], x[2])})
saveRDS(p2List, file = "p2list.rds")
p2ListPerm <- readRDS("p2list.rds")

remCol <- function(x, replaceZeros = TRUE) {
  
  #Removes and replaces all 1 and all 0 rows
  
  allOne <- matrix(c(replicate(nrow(x) - 1, 1), 0), nrow = nrow(x))
  allZero <- matrix(c(replicate(nrow(x) - 1, 0), 1), nrow = nrow(x))
  
  x[, which(colSums(x) == nrow(x))] <- allOne
  
  if(replaceZeros){
    x[, which(colSums(x) == 0)] <- allZero
  }
  
  return(x)
}

#Filter the lists to remove any samples of all 1's or all 0's, which
#result in inappropriate estimation of the expectation and variance.

p2ListFiltered <- lapply(p2ListPerm, function (x) {remCol(x, replaceZeros = FALSE)})

#Obtain the expectation and variances
expectation <- function(x) {
  mean(x)/(1-mean(x))
}
deltaVar <- function (x) {
  s <- length(x)
  m <- (1/(1-mean(x)))^4
  variance <- m * (mean(x)*(1-mean(x))) * (1/s)
  return(variance)
}

p2ExpDelta <- lapply(p2ListFiltered, function (y) {apply(y, 2, expectation)})
p2VarDelta <- lapply(p2ListFiltered, function (y) {apply(y, 2, deltaVar)})
p2ExpDeltaP <- lapply(p2ListPerm, function (y) {apply(y, 2, expectation)})
p2VarDeltaP <- lapply(p2ListPerm, function (y) {apply(y, 2, deltaVar)})

v_delta <- apply(sims, 1, function(x) { x[2] / (x[1] * (1 - x[2])^3)}, simplify = "array")

averageVarsFiltered <- lapply(p2VarDelta, function(x) {mean(x[is.finite(x)])})
averageVarsPerm <- lapply(p2VarDeltaP, function(x) {mean(x[is.finite(x)])})
averageExpFiltered <- lapply(p2ExpDelta, function(x) {mean(x[is.finite(x)])})
averageExpPerm  <- lapply(p2ExpDeltaP, function(x) {mean(x[is.finite(x)])})

cbind(averageVarsFiltered, averageVarsPerm, averageExpFiltered, averageExpPerm)

##Problem 2, Bootstrap Estimation

bootSamples <- function (x, b) {
    replicate(b, sample(x, length(x), replace = TRUE), simplify = TRUE)
}


bootstrapMeans <- function(samples) {
  
#We replace any all 1's samples with a column of n-1 1's and one 0.
  
  #samples <- remCol(samples, replaceZeros = FALSE)
  
  # Calculate the mean
  bSampleMeans <- apply(samples, 2, expectation)
  bSampleMeans <- bSampleMeans[is.finite(bSampleMeans)]
  if(length(bSampleMeans) == 0) return(NA)
  mean(bSampleMeans)
  
}

bootstrapVariance <- function(samples) {
  
  #We replace any all 1's samples with a column of n-1 1's and one 0.
  
  #samples <- remCol(samples, replaceZeros = FALSE)

  n <- nrow(samples)
  
  bSampleMeans <- apply(samples, 2, expectation)
  bSampleMeans <- apply(samples, 2, expectation)
  bSampleMeans <- bSampleMeans[is.finite(bSampleMeans)]
  if(length(bSampleMeans) < 2) return(NA)
  var(bSampleMeans)
  
}

#Gather bootstrap samples

p2BootSamples <- lapply(p2ListPerm, function (x){
    apply(x, 2, function (z) {
      bootSamples(z, 1000)
      }, simplify = FALSE)
    })

p2ExpBoot <- lapply(p2BootSamples, function (x){
              sapply(x, bootstrapMeans, simplify = TRUE)
            } 
          )



p2VarBoot <- lapply(p2BootSamples, function (x){
      sapply(x, bootstrapVariance, simplify = TRUE)
    } 
  )

bootVarAverage <- lapply(p2VarBoot, function(x) {mean(x[is.finite(x)])})
bootAverage <- lapply(p2ExpBoot, function(x) {mean(x[is.finite(x)])})

cbind(averageVarsFiltered, bootVarAverage)

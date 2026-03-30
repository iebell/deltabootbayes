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

#Filter the lists to remove any samples of all 1's, which
#result in non-existent estimation of the vmean and variance.

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

##Problem 3, Bayesian estimation

#Produce Bayesian Estimates of the Mean and Variance

betaEst <- function(x) {
  a <- sum(x == 1)
  b <- sum(x == 0)
  betaDraws <- rbeta(1000, (0.5 + a), (0.5 + b))
  oddsDraws <- betaDraws / (1-betaDraws)
  if ((b+0.5) < 2){
    c(mean = NA_real_, var = NA_real_)
  } else {
    c(mean(oddsDraws), var(oddsDraws))
  }
}

#Testing to see how often there is (at least) no posterior variance is available

testBayes <- lapply(p2ListPerm, function(x) {apply(x, 2, betaEst)})
totalNA <- lapply(testBayes, function(x) {sum(is.na(x[1,]))})

#Alternative method: report mode and 95% confidence interval with coverage rate

betaBad <- function(x) {
  
  a <- sum(x == 1) + 0.5
  b <- sum(x == 0) + 0.5
  m <- ifelse(a  < 1, 0, (a - 1)/(b + 1))
  bConf <- qbeta(c(0.025, 0.975), a, b)
  oddsConf <- bConf / (1 - bConf)
  return(c(mode = m, confLow = oddsConf[1], confHigh = oddsConf[2]))
  
}

badBayes <- lapply(p2ListPerm, function(x) {apply(x, 2, betaBad)})

#Part 5, Plots and graphs

##Box Plots of Means

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 1:4){
  boxplot(p2ExpDelta[[i]], p2ExpBoot[[i]], testBayes[[i]][1,], 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
  abline(h = (0.01)/(1-0.01), lty = 2)
}
mtext("Estimated Expectations of Odds, p = 0.01", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 5:8){
  boxplot(p2ExpDelta[[i]], p2ExpBoot[[i]], testBayes[[i]][1,], 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-4]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
  abline(h = (0.1)/(1-0.1), lty = 2)
}
mtext("Estimated Expectations of Odds, p = 0.10", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 9:12){
  boxplot(p2ExpDelta[[i]], p2ExpBoot[[i]], testBayes[[i]][1,], 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-8]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
  abline(h = (0.25)/(1-0.25), lty = 2)
}
mtext("Estimated Expectations of Odds, p = 0.25", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 13:16){
  boxplot(p2ExpDelta[[i]], p2ExpBoot[[i]], testBayes[[i]][1,], 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-12]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
  abline(h = (0.5)/(1-0.5), lty = 2)
}
mtext("Estimated Expectations of Odds, p = 0.50", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 17:20){
  boxplot(p2ExpDelta[[i]], p2ExpBoot[[i]], if(i == 17){badBayes[[17]][1,]}
          else{testBayes[[i]][1,]}, 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-16]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
  abline(h = (0.75)/(1-0.75), lty = 2)
}
mtext("Estimated Expectations of Odds, p = 0.75", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 21:24){
  boxplot(p2ExpDelta[[i]], p2ExpBoot[[i]], if(i == 22){badBayes[[22]][1,]}
          else if(i == 21){badBayes[[21]][1,]}
          else{testBayes[[i]][1,]}, 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-20]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
  abline(h = (0.9)/(1-0.9), lty = 2)
}
mtext("Estimated Expectations of Odds, p = 0.90", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 25:28){
  boxplot(p2ExpDelta[[i]], p2ExpBoot[[i]], if(i == 25){badBayes[[25]][1,]}
          else if(i == 26){badBayes[[26]][1,]}
          else if(i == 27){badBayes[[27]][1,]}
          else{testBayes[[i]][1,]}, 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-24]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
  abline(h = (0.99)/(1-0.99), lty = 2)
}
mtext("Estimated Expectations of Odds, p = 0.99", outer = TRUE, cex = 1.5, font = 2)

#Boxplots of Variances

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 1:4){
  boxplot(p2VarDelta[[i]], p2VarBoot[[i]], testBayes[[i]][2,], 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
}
mtext("Variances of Odds, p = 0.01", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 5:8){
  boxplot(p2VarDelta[[i]], p2VarBoot[[i]], testBayes[[i]][2,], 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-4]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
}
mtext("Variances of Odds, p = 0.10", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 9:12){
  boxplot(p2VarDelta[[i]], p2VarBoot[[i]], testBayes[[i]][2,], 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-8]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
}
mtext("Variances of Odds, p = 0.25", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 13:16){
  boxplot(p2VarDelta[[i]], p2VarBoot[[i]], testBayes[[i]][2,], 
          names = c("Delta Method", "Bootstrap", "Bayesian Analysis"), 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-12]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
}
mtext("Variances of Odds, p = 0.50", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 17:20){
  boxplot(p2VarDelta[[i]], p2VarBoot[[i]], if(i == 17){NULL}
          else{testBayes[[i]][2,]}, 
          names = if(i ==17){c("Delta Method", "Bootstrap", "Bayesian Analysis*")}
          else{c("Delta Method", "Bootstrap", "Bayesian Analysis")}, 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-16]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
}
mtext("Variances of Odds, p = 0.75", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 21:24){
  boxplot(p2VarDelta[[i]], p2VarBoot[[i]], if(i == 22 | i == 21){}
          else{testBayes[[i]][2,]}, 
          names = if(i == 22| i == 21){c("Delta Method", "Bootstrap", "Bayesian Analysis*")}
          else{c("Delta Method", "Bootstrap", "Bayesian Analysis")}, 
          main = paste0("n = ", c(10, 30, 100, 1000)[i-20]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
}
mtext("Variances of Odds, p = 0.90", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(2, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
for(i in 25:28){
  boxplot(p2VarDelta[[i]], p2VarBoot[[i]], if(i == 25 | i == 26 | i == 27){}
          else{testBayes[[i]][2,]}, 
          names = if(i == 25|i == 26|i == 27){c("Delta Method", "Bootstrap", "Bayesian Analysis*")}
          else{c("Delta Method", "Bootstrap", "Bayesian Analysis")},
          main = paste0("n = ", c(10, 30, 100, 1000)[i-24]), col = c("#5DA5A4", "#A3C4BC", "#E8D8C3"))
}
mtext("Variances of Odds, p = 0.99", outer = TRUE, cex = 1.5, font = 2)

# Table of Confidence Intervals and MAPs

n <- c(10, 10, 30, 10, 30, 100)
p <- c(0.75, 0.90, 0.90, 0.99, 0.99, 0.99)
odds <- p / (1-p)
averageMap <- c(mean(badBayes[[17]][1,]), mean(badBayes[[21]][1,]), mean(badBayes[[22]][1,]),
                mean(badBayes[[25]][1,]), mean(badBayes[[26]][1,]), mean(badBayes[[27]][1,]))
coverageRate <- function(b, odds){
  cover <- apply(b, 2, function (x){ifelse(x[2] < odds && odds < x[3],TRUE,FALSE)}, simplify = TRUE)
  return(sum(cover == TRUE) / 100)
}
cR <- c(coverageRate(badBayes[[17]], odds[1]), 
        coverageRate(badBayes[[21]], odds[2]), 
        coverageRate(badBayes[[22]], odds[3]),
        coverageRate(badBayes[[25]], odds[4]), 
        coverageRate(badBayes[[26]], odds[5]), 
        coverageRate(badBayes[[27]], odds[6])
        )
dfBb <- data.frame(
  "Sample Size" = n,
  "Probability" = p,
  "Odds" = odds,
  "Mode" = averageMap,
  "Coverage Rate" = cR
)
  
#Delta Method Histograms

par(mfrow = c(3, 2), mar=c(4, 4, 2, 1), oma=c(0, 0, 4,0))
dig <- seq(-10, 25, length.out = 1000)
norm <- dnorm(dig, 0, sqrt(16*(0.5)*0.5))
hist(p1Delta[[1]], main = "n = 10", xlab = "", probability = TRUE, breaks = 25)
lines(dig, norm, col = "red", lwd = 2)
hist(p1Delta[[2]], main = "n = 30", breaks = 10, xlab = "", probability = TRUE)
lines(dig, norm, col = "red", lwd = 2)
hist(p1Delta[[3]], main = "n = 50", breaks = 20, xlab = "", probability = TRUE)
lines(dig, norm, col = "red", lwd = 2)
hist(p1Delta[[4]], main = "n = 100", breaks = 20, xlab = "", probability = TRUE)
lines(dig, norm, col = "red", lwd = 2)
hist(p1Delta[[5]], main = "n = 500", breaks = 20, xlab = "", probability = TRUE)
lines(dig, norm, col = "red", lwd = 2)
hist(p1Delta[[6]], main = "n = 1000", breaks = 50, xlab = "", probability = TRUE)
lines(dig, norm, col = "red", lwd = 2)
mtext("Histograms of √N (g(X̄) − g(p))", outer = TRUE, cex=1, font = 2)


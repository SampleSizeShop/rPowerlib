#
# Routines to simulate power for a general linear
# hypothesis for a general linear
# multivariate model with one or more Gaussian
# predictors.  
#
# Author: Sarah Kreidler
#
#
#
library(car)
set.seed(1066)

getTestResults <- function(x){
	test <- x$test
	if ((!is.null(x$singular)) && x$singular){
		warning("the error SSP matrix is singular; multivariate tests are unavailable")
		return(invisible(x))
	}
	SSPE.qr <- qr(x$SSPE)
	# the following code is adapted from summary.manova
	eigs <- Re(eigen(qr.coef(SSPE.qr, x$SSPH), symmetric = FALSE)$values)
	tests <- matrix(NA, 4, 4)
	rownames(tests) <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
	if ("Pillai" %in% test)
		tests[1, 1:4] <- stats:::Pillai(eigs, x$df, x$df.residual)
	if ("Wilks" %in% test)
		tests[2, 1:4] <- stats:::Wilks(eigs, x$df, x$df.residual)
	if ("Hotelling-Lawley" %in% test)
		tests[3, 1:4] <- stats:::HL(eigs, x$df, x$df.residual)
	if ("Roy" %in% test)
		tests[4, 1:4] <- stats:::Roy(eigs, x$df, x$df.residual)
	tests <- na.omit(tests)
	ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
	ok <- !is.na(ok) & ok
	tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3], tests[ok, 4],
					lower.tail = FALSE))
	colnames(tests) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)")
	tests <- structure(as.data.frame(tests),
			heading = paste("\nMultivariate Test",
					if (nrow(tests) > 1) "s", ": ", x$title, sep=""),
			class = c("anova", "data.frame"))
	return(tests)
}



#
# Generate a sample from a multivariate normal
# of dimension rows x cols with the specified
# correlation structure
#
#
getMultivariateGaussianSample <- function(rows, correlation) {
  cols = ncol(correlation)
  cholMatrix = chol(correlation)
  randomNormalMatrix = matrix(nrow=rows,rnorm(rows*cols));
  return (randomNormalMatrix %*% cholMatrix)
}

#
# Get a single instance of the design matrix
# by generating the random predictors based
# on the specified correlation
#
getDesignMatrix <- function(XFixed,sigmaG) {
  randomX = getMultivariateGaussianSample(nrow(XFixed),sigmaG)
  return (cbind(XFixed,randomX))
}

#
# Build the error matrix
#
#
getSigmaError <- function(sigmaGInverse,sigmaY,sigmaYG) {
  sigmaGY = t(sigmaYG)
  sigmaE = sigmaY - (sigmaYG %*% sigmaGInverse %*% sigmaGY);
  return(sigmaE)
}

#
# Simulate errors and fit the model
#
simulateAndFitModel <- function(typeIError, N, X, Beta, sigmaE, C, U) {
  E = getMultivariateGaussianSample(N,sigmaE)
  Y = X %*% Beta + E
  model.fit = lm(Y ~ X-1)
  glh.test = linearHypothesis(model.fit,hypothesis.matrix=C,test="Hotelling-Lawley",P=U)
  pvalue = getTestResults(glh.test)[1,6]
  return (as.numeric(pvalue <= typeIError))
}

#
# Simulate power
#
simulatePower <- function(typeIError, XFixed,
	Beta, sigmaG, sigmaE, C, U, numIterations=1000) {

  X = getDesignMatrix(XFixed,sigmaG)
  N = nrow(X)

  # for each sample, simulate errors for 
  # numIterations iterations.
  #
  numRejections = 0
  for(i in 1:numIterations) {
    numRejections = numRejections + simulateAndFitModel(typeIError, N, X, Beta,
			  sigmaE, C, U)
  }

  return (numRejections / numIterations)
}

#
# Simulate unconditional power.
#
#
#
#
simulateUnconditionalPower <- function(typeIError=0.05,
      XFixed,betaFixed,
	sigmaG,sigmaY,sigmaYG,
	C,U,numSamples=1000,numIterations=1000) {

  # build the random component of beta
  sigmaGInverse = solve(sigmaG)
  sigmaGY = t(sigmaYG)
  Beta = rbind(betaFixed,(sigmaGInverse %*% sigmaGY))
  # get the overall error covariance for the model
  sigmaE = sigmaY - (sigmaYG %*% sigmaGInverse %*% sigmaGY)

  # generate several possible samples and simulate power
  powerSample = replicate(numSamples,
			  simulatePower(typeIError, XFixed, Beta,
			  sigmaG, sigmaE, C, U, numIterations))
  return(powerSample)

}


#
#
#
unconditionalPower.example1 <- function() {
  essX = diag(2)
  XFixed = essX %x% matrix(rep(1,10))
  sigmaG = diag(c(1,4))
  sigmaY = matrix(ncol=3,byrow=1,c(2,0.4,0.2,0.4,2,0.4,0.2,0.4,2))
  sigmaYG = matrix(ncol=2,c(0.2,0.2,0.2,0.1,0.1,0.1)) %*% diag(c(2,8))

  betaFixed = matrix(ncol=3,byrow=1,c(2,0,0,0,0,0))
  C = matrix(ncol=4,byrow=1,c(1,-1,0,0))
  U = t(cbind(matrix(rep(1,2)),-1*diag(2)))

  powerSample = simulateUnconditionalPower(0.05,
      XFixed,betaFixed,
	sigmaG,sigmaY,sigmaYG,
	C,U,1000,1000)
}


#
# test code
#

getMultivariateGaussianSample.test <- function() {
  sigmaG = matrix(ncol=2,byrow=1,c(0.2,0.1,0.1,0.4))
  sample = getMultivariateGaussianSample(100,sigmaG)
  sample
  cov(sample)
}

getDesignMatrix.test <- function() {
  sigmaG = matrix(ncol=2,byrow=1,c(0.2,0.1,0.1,0.4))
  essX = diag(2)
  XFixed = essX %x% matrix(rep(1,10))
  sample = getDesignMatrix(XFixed,sigmaG)
  print(sample)
  cov(sample[,3:4])
}
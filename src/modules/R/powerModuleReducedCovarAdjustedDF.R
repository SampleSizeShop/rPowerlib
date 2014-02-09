#
# Routines to calculate power for a GLMM(F,G) using
# the proposed Kreidler method  
#
# 1. Set the covariance to
# sigmaE = sigmaY - (sigmaYG %*% sigmaGInverse %*% sigmaGY)
# 2. Calculate power for the HLT using XFixed, sigmaE, etc.
#    but with degrees of freedom based on the full X
#
# Author: Sarah Kreidler
#
#
#

source("powerGLMMFGExampleStudyDesigns.R")


#
# Build the error matrix
#
#
getSigmaError <- function(sigmaG,sigmaY,sigmaYG) {
  sigmaGInverse = solve(sigmaG)
  sigmaGY = t(sigmaYG)
  sigmaE = sigmaY - (sigmaYG %*% sigmaGInverse %*% sigmaGY);
  return(sigmaE)
}



powerPillaiBartlettTrace <- function(alpha, N, nMinusR, a, b, eInverseH) {

  # calculate numerator degrees of freedom
  ndf = a*b
  # calculate denominator degrees of freedon
  s = min(a,b)
  ddf = s*(nMinusR+s-b)
  # get critical F
  Fcrit = qf(1-alpha,ndf,ddf)

  # get Pillai Bartlett trace
  pbt = sum(diag(eInverseH %*% solve(diag(b)+eInverseH)))

  # get the noncentrality
  ncp = N*s*pbt/(s-pbt)

  # calculate power
  power = 1 - pf(Fcrit,ndf,ddf,ncp)

  return(data.frame(test="Pillai Bartlett Trace",power,Fcrit,ndf,ddf,ncp))
}

powerHotellingLawleyTracePillaiSamson <- function(alpha, N, nMinusR, a, b, eInverseH) {


  # calculate numerator degrees of freedom
  ndf = a*b
  # calculate denominator degrees of freedon
  s = min(a,b)
  ddf = s*(nMinusR-b-1)+2
  # get critical F
  Fcrit = qf(1-alpha,ndf,ddf)

  # get Hotelling Lawley trace
  hlt = sum(diag(eInverseH))

  # get the noncentrality
  ncp = N*hlt

  # calculate power
  power = 1 - pf(Fcrit,ndf,ddf,ncp)

  return(data.frame(test="Hotelling Lawley Trace (Pillai & Samson)",power,Fcrit,ndf,ddf,ncp))

}

powerHotellingLawleyTraceMcKeon <- function(alpha, N, nMinusR, a, b, eInverseH) {

  # calculate numerator degrees of freedom
  ndf = a*b
  # calculate denominator degrees of freedom
  t1 = nMinusR^2 - nMinusR *(2 * b + 3) + b * (b + 3);
  t2 = (nMinusR * (a + b + 1) - (a + 2 * b + b * b - 1));
  ddf = 4 + (a * b + 2) * (t1/t2);

  # get critical F
  Fcrit = qf(1-alpha,ndf,ddf)

  # get Hotelling Lawley trace
  hlt = sum(diag(eInverseH))

  # get the noncentrality
  ncp = N*hlt

  # calculate power
  power = 1 - pf(Fcrit,ndf,ddf,ncp)

  return(data.frame(test="Hotelling Lawley Trace (McKeon)",power,Fcrit,ndf,ddf,ncp))
}

#
# Power for the Wilks lambda using the method
# of Shieh
#
powerWilks <- function(alpha, N, nMinusR, a, b, eInverseH) {

  # get numerator degrees of freedom
  ndf = a*b

  # get denominator degrees of freedom
  t = 1
  if (ndf >= 4) {
    t = sqrt((ndf^2-4)/(a^2+b^2-5))
  }
  ddf = t*(nMinusR-(b-a+1)/2)-(ndf-2)/2

  # get critical F
  Fcrit = qf(1-alpha,ndf,ddf)

  lambda = det(solve(diag(b) + eInverseH))
  ncp = N*t*(lambda^(-1/t)-1)
  power = 1 - pf(Fcrit,ndf,ddf,ncp)

  return(data.frame(test="Wilks' Lambda",power,Fcrit,ndf,ddf,ncp))

}






#
# Calculate power for a general linear multivariate
# model with one or more random predictors
#
calculatePowerGLMMFG <- function(alpha, XFixed, betaFixed,
	sigmaG, sigmaY, sigmaYG,
	C, U, thetaNull) {

  # get the reduced covariance of errors
  sigmaE = getSigmaError(sigmaG,sigmaY,sigmaYG)

  # initialize some dimensions
  N = nrow(XFixed)
  r = nrow(betaFixed) + ncol(sigmaG)
  p = ncol(betaFixed)
  a = nrow(C)
  b = ncol(U)
  s = min(a,b)
  ab = a*b

  # adjust the degrees of freedom to match the true
  # rank of the design matrix
  nMinusR = N - nrow(betaFixed)

  # calculate the observed theta
  theta = C %*% betaFixed %*% U

  # calculate the sum of squares/products for the error
  SSPE = t(U) %*% sigmaE %*% U
  # calculate the sum of squares/ products for the hypothesis
  SSPH = t(theta - thetaNull) %*% solve(C %*% solve(t(XFixed)%*%XFixed) %*% t(C)) %*% (theta - thetaNull)

  # HE^-1
  eInverseH = solve(SSPE) %*% SSPH

  ## create an empty data frame to hold the results
  power = data.frame(test=character(0),power=double(0),Fcrit=double(0),
			   ndf=double(0),ddf=double(0),ncp=double(0))
  ## calculate power for each of the tests
  power = rbind(
	powerWilks(alpha, N, nMinusR, a, b, eInverseH),
	powerPillaiBartlettTrace(alpha, N, nMinusR, a, b, eInverseH),
	powerHotellingLawleyTracePillaiSamson(alpha, N, nMinusR, a, b, eInverseH),
	powerHotellingLawleyTraceMcKeon(alpha, N, nMinusR, a, b, eInverseH)
	)
  return(power)
}



#
# Unit test
#
calculatePowerGLMMFG.test <- function() {

  design1 = powerGLMMFG.example1()
  calculatePowerGLMMFG(design@alpha, 
	design@XFixed, 
	design@betaFixed, 
	design@sigmaG, 
	design@sigmaY, 
	design@sigmaYG, 
	design@CMatrix, 
	design@UMatrix, 
	design@thetaNull)

}
calculatePowerGLMMFG.test()


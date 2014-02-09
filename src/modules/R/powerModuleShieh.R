#
# Power for the GLMM(F,G) using the method of Shieh.  Calculations
# based on the algorithm described in
#
# Shieh, G. (2007). A unified approach to power calculation
# and sample size determination for random regression models. 
# Psychometrika, 72(3), 347-360.
#
# Notation based on 
#
# Muller, K. E., & Stewart, P. W. (2006). Linear model theory: 
# univariate, multivariate, and mixed models. Hoboken, 
# New Jersey: John Wiley and Sons.
#
#
#


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


calculatePowerShieh <- function(alpha, N, KStar, Beta, C, U, sigmaE, thetaNull) {

  # initialize some dimensions
  r = nrow(Beta)
  p = ncol(Beta)
  a = nrow(C)
  b = ncol(U)
  s = min(a,b)
  ab = a*b
  nMinusR = N - r

  # calculate the observed theta
  theta = C %*% Beta %*% U

  # calculate the sum of squares/products for the error
  SSPE = t(U) %*% sigmaE %*% U
  # calculate the sum of squares/ products for the hypothesis
  SSPH = t(theta - thetaNull) %*% solve(C %*% solve(KStar) %*% t(C)) %*% (theta - thetaNull)

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
calculatePowerShieh.test <- function() {

 alpha = 0.05
 N = 110

 C = matrix(ncol=4,byrow=1,c(0,1,0,0,0,0,1,0,0,0,0,1));
 U = cbind((1/sqrt(2))*matrix(c(-1,0,1)),(1/sqrt(6))*matrix(c(1,-2,1)))
 
 Beta=matrix(ncol=3,byrow=1,c(114.46,104.66,98.83,2.88,8.77,10.67,-0.71,-0.90,-1.30,-0.21,-0.54,-0.72));
 sigmaE=matrix(ncol=3,byrow=1,c(218.48, 83.66, 72.19, 83.66, 251.92, 158.60, 72.19, 158.60, 244.58));
 KStar=matrix(ncol=4,byrow=1,c(1, 0, 1, 0, 0, 1, 0, 3, 1, 0, 3, 0, 0, 3, 0, 15));
 thetaNull = matrix(ncol=2,c(0,0,0,0,0,0))
 calculatePowerShieh(alpha, N, KStar, Beta, C, U, sigmaE, thetaNull)


}
calculatePowerShieh.test()
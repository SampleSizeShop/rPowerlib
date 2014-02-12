
#
# design.glmmFG
#
# Class describing the general linear model with fixed predictors and
# one or more Gaussian covariates
#
setClass (
  "design.glmmFG",
  representation ( name = "character",
                   description = "character",
                   XEssence = "matrix",
                   perGroupN = "numeric",
                   Beta = "matrix",
                   SigmaY = "matrix",
                   SigmaG = "matrix",
                   SigmaYG = "matrix"
  ),
  prototype ( name ="",
              description ="",
              XEssence = diag(2),
              perGroupN = 10,
              Beta = matrix(c(1,0),nrow=2),
              SigmaY = matrix(c(1)),
              SigmaG = matrix(c(1)),
              SigmaYG = matrix(c(0.5))
  ),
  validity = function(object) {
    #
    # Note, the class definition will already enforce that the
    # matrices are non-null 
    #
    
    # make sure that Sigma Y is square, symmetric and positive definite
    if (!isSymmetric(object@SigmaY)) {
      stop("SigmaY matrix is not symmetric")
    } else {
      sigmaYEigenValues = eigen(object@SigmaY)
      if (sum(sigmaYEigenValues$values > 0) < length(sigmaYEigenValues$values)) {
        stop("SigmaY matrix is not positive definite")
      }
    }
    # make sure that Sigma G is square, symmetric and positive definite
    if (!isSymmetric(object@SigmaG)) {
      stop("SigmaG matrix is not symmetric")
    } else {
      sigmaGEigenValues = eigen(object@SigmaG)
      if (sum(sigmaGEigenValues$values > 0) < length(sigmaGEigenValues$values)) {
        stop("SigmaY matrix is not positive definite")
      }
    }
    
    # check matrix conformance
    if (ncol(object@XEssence) != nrow(object@Beta)) {
      stop("The number of columns in the essence X matrix does not match the 
           number of columns in the Beta matrix")
    } else if (ncol(object@Beta) != nrow(object@SigmaY)) {
      stop("The number of columns in the Beta matrix does not match the 
           number of rows in the SigmaY matrix")    
    } else if (ncol(object@SigmaYG) != nrow(object@SigmaG)) {
      stop("The number of columns in the SigmaYG matrix does not match the
           number of rows in the SigmaG matrix")
    } else if (nrow(object@SigmaYG) != nrow(object@SigmaY)) {
      stop("The number of rows in the SigmaYG matrix does not match the
           number of rows in the SigmaY matrix")
    }
    
    return(TRUE)
    }
      )




#
# glh
#
# Class describing the general linear hypothesis
#
setClass (
  "glh",
  representation ( alpha = "numeric",
                   betweenContrast = "matrix",
                   withinContrast = "matrix",
                   thetaNull = "matrix",
                   test = "character"
  ),
  prototype ( alpha = 0.05,
              betweenContrast = matrix(c(1,-1), nrow=1),
              withinContrast = matrix(c(1)),
              thetaNull = matrix(c(0)),
              test = "Hotelling-Lawley"
  ),
  validity = function(object) {
    # make sure thetaNull conforms with the between and within contrasts
    if (nrow(object@betweenContrast) != nrow(object@thetaNull)) {
      stop("The number of rows in the between contrast must match the number of rows of thetaNull")
    }
    if (ncol(object@withinContrast) != ncol(object@thetaNull)) {
      stop("The number of columns in the within contrast must match the number of columns of thetaNull")
    }
    
    return(TRUE)
  }
)



#
# Adjustment factor for the M matrix based on the
# expected value of the projection matrix 
# of the covariates: G(G'G)^-1G'
#
mAdjust.fixedVsRandomDF <- function(totalN, numFixedPredictors, numRandomPredictors) {
  return((totalN - numFixedPredictors) / (totalN - numFixedPredictors - numRandomPredictors))
}

#
# Adjustment factor for the M matrix based on the
# expected value of the projection matrix 
# of the covariates: G(G'G)^-1G'
#
mAdjust.expectedProjection <- function(totalN, numFixedPredictors, numRandomPredictors) {
  return ((totalN) / (totalN - numRandomPredictors)) 
}

#
# Adjustment factor of 1, that is, no adjustment to M
#
mAdjust.noAdjust <- function(totalN, numFixedPredictors, numRandomPredictors) { return(1) } 

#
# Calculate power using the covariate adjustment approach proposed
# by Kreidler et al.
#
glmmPower.covariateAdjusted = function(design, hypothesis, mAdjust=mAdjust.expectedProjection) {
  # TODO: make sure model and glh conform
  if (class(design) != "design.glmmFG") {
    stop("The specified design is not a linear model with fixed and random covariates")
  }
  if (class(hypothesis) != "glh") {
    stop("the specified hypothesis is not a general linear hypothesis")
  }
  if (hypothesis@test != "Hotelling-Lawley") {
    stop("At present, this function only supports the Hotelling Lawley Trace test")
  }
  
  # form sigmaE
  SigmaE = design@SigmaY - design@SigmaYG %*% solve(design@SigmaG) %*% t(design@SigmaYG)
  
  # get total sample size and nuE
  totalN = design@perGroupN * nrow(design@XEssence)
  rank = qr(design@XEssence)$rank
  nuE = totalN - (rank + ncol(design@SigmaG))
  
  # calculate the fixed portion of X'X inverse
  XtXinverse = (1/design@perGroupN) * solve(t(design@XEssence) %*% design@XEssence)
  
  #XtXinverse = XtXinverse - totalN * design@
  #
  # variable naming conventions are from Muller & Stewart
  #
  #### NOTE: at present, we only support the Hotelling Lawley Trace test ####
  
  # calculate the numerator degrees of freedom
  a = nrow(hypothesis@betweenContrast)
  b = ncol(hypothesis@withinContrast)
  s = min(a,b)
  ndf = a * b
  
  # calculate the denominator degrees of freedom
  t1 = nuE * nuE - nuE * (2 * b + 3) + b * (b + 3);
  t2 = (nuE * (a  + b + 1) - (a + 2 * b + b * b - 1));
  ddf = 4 + (a * b + 2) * (t1/t2);
  
  # calculate the error sums of squares 
  SSPE = (t(hypothesis@withinContrast) %*% SigmaE %*% hypothesis@withinContrast) * nuE
  # calculate the hypothesis sum of squares for the fixed portion of the design
  C = matrix(hypothesis@betweenContrast[1:nrow(hypothesis@betweenContrast),1:nrow(design@Beta)], nrow=nrow(hypothesis@betweenContrast))
  
  thetaObserved = C %*% design@Beta %*% hypothesis@withinContrast
  thetaDiff = thetaObserved - hypothesis@thetaNull
  M = C %*% XtXinverse %*% t(C) * mAdjust(totalN, rank, ncol(design@SigmaG))
  
  SSPH = t(thetaDiff) %*% solve(M) %*% thetaDiff
  # calculate the Hotelling-Lawley trace
  hlt = sum(diag(SSPH %*% solve(SSPE)))
  
  # calculate the critical F
  Fcrit = qf(1 - hypothesis@alpha, ndf, ddf)
  
  # compute the noncentrality
  omega = hlt * nuE
  
  # calculate power
  power = 1 - pf(Fcrit, ndf, ddf, ncp=omega)
  
  return(power)
}




#
# Define the study design for the proposed trial of
# salivary biomarkers in oral cancer
#
design.elashoff = new("design.glmmFG", 
                      XEssence=diag(2),
                      perGroupN=35,
                      Beta=matrix(c(2.14,2.14,0,2.14), nrow=2), 
                      SigmaY=4.84*matrix(c(1,0.4,0.4,1), nrow=2), 
                      SigmaG=diag(c(4.84,161.29)), 
                      SigmaYG=matrix(c(1.94,2.79,0.97,2.79), nrow=2, byrow=TRUE))

#
# Define the hypothesis
#
glh.elashoff = new("glh", betweenContrast=matrix(c(1,-1,0,0),nrow=1), 
                   withinContrast=matrix(c(1,-1), nrow=2), 
                   thetaNull=matrix(c(0)))

# Calculate power, adjusted for covariates
power = glmmPower.covariateAdjusted(design.elashoff, glh.elashoff)
# print the result
power

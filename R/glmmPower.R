# 
#  Package glmmPower calculates power for the general linear 
#  multivariate model
#  Copyright (C) 2013 Sarah Kreidler.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# 

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
glmmPower.covariateAdjusted = function(design, hypothesis, mAdjust=mAdjust.fixedVsRandomDF) {
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
# Calculate power based on the asymptotic approach 
# developed by Shieh
#
glmmPower.shieh = function(design, glh) {
  # validate
  # TODO: make sure model and glh conform
  if (class(design) != "design.glmmFG") {
    stop("The specified design is not a linear model with fixed and random covariates")
  }
  if (class(glh) != "glh") {
    stop("the specified hypothesis is not a general linear hypothesis")
  }
  
  # form sigmaE
  sigmaE = design@SigmaY - design@SigmaYG %*% solve(design@SigmaG) %*% t(design@SigmaYG)
  
  # get total sample size and nuE
  totalN = design@perGroupN * nrow(design@XEssence)
  rank = qr(design@XEssence)$rank
  nuE = totalN - (rank + ncol(design@SigmaG))
  
  # determine the Kstar matrix from the design
  XtX = t(design@XEssence) %*% design@XEssence
  Kstar = adiag(design@perGroupN * XtX, totalN * design@SigmaG)
  
  # form the error and hypothesis sum of squares 
  betaFull = rbind(design@Beta, solve(design@SigmaG) %*% t(design@SigmaYG))  
  E = t(glh@withinContrast) %*% sigmaE %*% glh@withinContrast * nuE
  thetaDiff = (glh@betweenContrast %*% betaFull %*% glh@withinContrast) - glh@thetaNull
  H = (t(thetaDiff) %*% 
         solve(glh@betweenContrast %*% solve(Kstar) %*% t(glh@betweenContrast)) %*% 
         thetaDiff)
  EinvH = solve(E) %*% H
  
  # calculate some components of the target F distribution parameters
  a = nrow(glh@betweenContrast)
  b = ncol(glh@withinContrast)
  ab = a*b 
  s = min(a,b)
  
  # calculate the numerator degrees of freedom
  ndf = ab 
  # calculate the denominator degrees of freedom,
  # and the noncentrality parameter for the given test
  if (glh@test == "Wilks Lambda") {
    T = ifelse(ab < 4, 1, sqrt((ab^2 - 4)/(a^2 + b^2 - 5)))
    ddf = T * (nuE - ((b - a + 1) / 2)) - ((ab - 2)/2);
    lambda = det(solve(diag(b) + EinvH))
    noncentrality = totalN * T * (lambda^(-1/T) - 1)
    
  } else if (glh@test == "Pillai-Bartlett Trace") {
    ddf = s * (nuE + s - b)
    pbt = sum(diag(EinvH %*% solve(diag(b) + EinvH)))
    noncentrality = totalN * (s * pbt / (s - pbt)) 
    
  } else if (glh@test == "Hotelling-Lawley Trace, Pillai-Sampson 1959") {
    ddf = s * (nuE - b - 1) + 2
    hlt = sum(diag(EinvH))
    noncentrality = nuE * hlt
    
  } else if (glh@test == "Hotelling-Lawley") {
    g = ((nuE^2 - nuE*(2*b + 3) + b*(b+3)) / 
           (nuE*(a + b + 1) - (a + 2*b + b^2 - 1))
    )
    ddf = 4 + (ab + 2)*g
    hlt = sum(diag(EinvH))
    noncentrality = nuE * hlt  
  }
  
  # get the critical value under the null
  Fcrit = qf(1 - glh@alpha, ndf, ddf)
  
  # get the power under the alternative
  return(1 - pf(Fcrit, ndf, ddf, noncentrality))
  
}


#
# Calculate power for linear models with a single
# Gaussian covariate
#
glmmPower.unconditionalSingleCovariate = function(design, hypothesis) {
  # validate
  # TODO: make sure model and glh conform
  if (class(design) != "design.glmmFG") {
    stop("The specified design is not a linear model with fixed and random covariates")
  }
  if (class(hypothesis) != "glh") {
    stop("the specified hypothesis is not a general linear hypothesis")
  }
  
  obj=.jnew("com/kreidles/PowerCalculator")
  power = .jcall(obj, "D", "calculateUnconditionalSingleCovariatePower", 
                 designToJSON(design, FALSE), toJSON(hypothesis))
  return(power)
}

#
# Calculate power for linear models with fixed predictors
# only
#
glmmPower.fixed = function(design, hypothesis) {
  # validate
  # TODO: make sure model and glh conform
  if (class(design) != "design.glmmF") {
    stop("The specified design is not a linear model with fixed predictors only")
  }
  if (class(hypothesis) != "glh") {
    stop("the specified hypothesis is not a general linear hypothesis")
  }
  
  obj=.jnew("com/kreidles/PowerCalculator")
  power = .jcall(obj, "D", "calculateConditionalPower", 
                 designToJSON(design, FALSE), toJSON(hypothesis))
  return(power)
}


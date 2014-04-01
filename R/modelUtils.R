#####################################################################
# 
#  Package rPowerlib calculates power for the general linear 
#  multivariate model with and without Gaussian covariates
#  Copyright (C) 2014 University of Colorado Denver.
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
#####################################################################

### Copied from car package since not exported
Pillai <- function (eig, q, df.res) {
  test <- sum(eig/(1 + eig))
  p <- length(eig)
  s <- min(p, q)
  n <- 0.5 * (df.res - p - 1)
  m <- 0.5 * (abs(p - q) - 1)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * n + s + 1
  c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s * tmp2)
}

Wilks <- function (eig, q, df.res) {
  test <- prod(1/(1 + eig))
  p <- length(eig)
  tmp1 <- df.res - 0.5 * (p - q + 1)
  tmp2 <- (p * q - 2)/4
  tmp3 <- p^2 + q^2 - 5
  tmp3 <- if (tmp3 > 0) 
    sqrt(((p * q)^2 - 4)/tmp3)
  else 1
  c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q, 
    p * q, tmp1 * tmp3 - 2 * tmp2)
}

HL <- function (eig, q, df.res) {
  test <- sum(eig)
  p <- length(eig)
  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (df.res - p - 1)
  s <- min(p, q)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * (s * n + 1)
  c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
}

Roy <- function (eig, q, df.res) {
  p <- length(eig)
  test <- max(eig)
  tmp1 <- max(p, q)
  tmp2 <- df.res - tmp1 + q
  c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
}

#
# utility function to extract p-value from linearHypothesis call
#
#getTestResults <- function(x) {
getTestResults <- function(test, SSPH, SSPE, df, df.residual) {
  
  test = test
  SSPE.qr <- qr(SSPE)
  # the following code is adapted from summary.manova
  eigs <- Re(eigen(qr.coef(SSPE.qr, SSPH), symmetric = FALSE)$values)
  
  if ("Pillai" %in% test)
    tests <- Pillai(eigs, df, df.residual)
  else if ("Wilks" %in% test)
    tests <- Wilks(eigs, df, df.residual)
  else if ("Hotelling-Lawley" %in% test)
    tests <- HL(eigs, df, df.residual)
  else if ("Roy" %in% test)
    tests <- Roy(eigs, df, df.residual)
  ok <- tests[2] >= 0 & tests[3] > 0 & tests[4] > 0
  ok <- !is.na(ok) & ok
  testResults <- data.frame(test=test,
                            Df=df, 
                            testStat=tests[1],
                            approxF=tests[2],
                            ndf=tests[3],
                            ddf=tests[4], 
                            ProbF=pf(tests[2], tests[3], tests[4],
                                     lower.tail = FALSE))
  
  return (testResults)
}

#
# Fit the model and test the specified glh
# using a cached (X'X)^-1
#
# This function replaces lm() in the simulation because
# it was just too slow.
#
#
fitModelAndTest <- function(XtXInverse, X, Y, glh, N, rank) {
  # fit the model
  betaHat = XtXInverse %*% t(X) %*% Y
  yHat = X %*% betaHat
  yDiff = Y - yHat
  sigmaE = t(yDiff) %*% yDiff / (N - rank)
  mInv = solve(glh@betweenContrast %*% XtXInverse %*% t(glh@betweenContrast))
  
  # calculate the hypothesis sum of squares
  thetaHat = glh@betweenContrast %*% betaHat %*% glh@withinContrast
  thetaDiff = thetaHat - glh@thetaNull
  SSPH = t(thetaDiff) %*% mInv %*% thetaDiff
  
  # calculate the error sum of squares
  SSPE = (N-rank) * (t(glh@withinContrast) %*% sigmaE %*% glh@withinContrast)
  
  return(getTestResults(test, SSPH, SSPE, nrow(glh@betweenContrast), N-rank))
  
}



#
# Convert a matrix to JSON
#
#
matrixToJSON <- function(name, m) {
  return(
    paste(
      c(
        "\"", name, "\": ",
        "{\"rows\": ", nrow(m), ",",
        "\"columns\": ", ncol(m), ",",
        "\"data\": [",
        paste(apply(m, 1, function(x) { 
          return(paste("[", paste(x, collapse=","), "]", collapse=""))
        }), collapse=","),
        "]}"
      ),
      collapse=""
    )
  )
}
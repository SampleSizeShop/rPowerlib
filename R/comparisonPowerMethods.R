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
library(magic)
library(car)

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
getTestResults <- function(x) {
  test = x$test
  SSPE.qr <- qr(x$SSPE)
  # the following code is adapted from summary.manova
  eigs <- Re(eigen(qr.coef(SSPE.qr, x$SSPH), symmetric = FALSE)$values)

  if ("Pillai" %in% test)
    tests <- Pillai(eigs, x$df, x$df.residual)
  else if ("Wilks" %in% test)
    tests <- Wilks(eigs, x$df, x$df.residual)
  else if ("Hotelling-Lawley" %in% test)
    tests <- HL(eigs, x$df, x$df.residual)
  else if ("Roy" %in% test)
    tests <- Roy(eigs, x$df, x$df.residual)
  ok <- tests[2] >= 0 & tests[3] > 0 & tests[4] > 0
  ok <- !is.na(ok) & ok
  testResults <- data.frame(test=x$test,
                            Df=x$df, 
                            testStat=tests[1],
                            approxF=tests[2],
                            ndf=tests[3],
                            ddf=tests[4], 
                            ProbF=pf(tests[2], tests[3], tests[4],
                                 lower.tail = FALSE))
  
  return (testResults)
}

###### end copied from car

glmmPower.glueck = function(design, hypothesis) {
  # get the most strongly associated covariate
  
  # form the new SigmaG and SigmaYg
  
  # build the json string to send to glimmpse
  
  # issue a power request to GLIMMPSE
}

#
# Calculate power based on the asymptotic approach 
# developed by Shieh
#
glmmPower.shieh = function(design, glh) {
  # validate
  # TODO: make sure model and glh conform
  if (class(design) != design.glmmFG) {
    stop("The specified design is not a linear model with fixed and random covariates")
  }
  if (class(glh) != "glh") {
    stop("the specified hypothesis is not a general linear hypothesis")
  }
  
  # form sigmaE
  sigmaE = design@SigmaY - design@SigmaYG %*% solve(design@SigmaG) %*% t(design@SigmaYG)
  
  # get total sample size and nuE
  totalN = design@perGroupN * nrow(design@XEssence)
  nuE = totalN - (rank + ncol(design@SigmaG))
  
  # determine the Kstar matrix from the design
  XtX = t(design@XEssence) %*% design@XEssence
  rank = qr(XtX)$rank
  Kstar = totalN * adiag(XtX, rank * design@SigmaG)
  
  # form the error and hypothesis sum of squares 
  E = t(glh@withinContrast) %*% sigmaE %*% glh@withinContrast
  thetaObs = (glh@betweenContrast %*% design@Beta %*% glh@withinContrast) - glh@thetaNull
  H = (t(thetaObs) %*% 
         solve(glh@betweenContrast %*% solve(Kstar) %*% t(glh@betweenContrast)) %*% 
         thetaObs)
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
    ddf = T * (nuE − (b − a + 1) / 2) − ((ab − 2)/2);
    lambda = det(solve(diag(b) + EinvH))
    noncentrality = totalN * T * (lambda^(-1/T) - 1)
    
  } else if (glh@test == "Pillai-Bartlett Trace") {
    ddf = s * (nuE + s - b)
    pbt = sum(diag(EinvH %*% solve(diag(b) + EinvH)))
    noncentrality = totalN * (s * pbt / (s - pbt)) 
    
  } else if (glh@test == "Hotelling-Lawley Trace, Pillai-Sampson 1959") {
    ddf = s * (nuE - b - 1) + 2
    hlt = sum(diag(EinvH))
    noncentrality = totalN * hlt

  } else if (glh@test == "Hotelling-Lawley Trace") {
    g = ((nuE^2 - nuE*(2*b + 3) + b*(b+3)) / 
           (nuE*(a + b + 1) - (a + 2*b + b^2 - 1))
    )
    ddf = 4 + (ab + 2)*g
    hlt = sum(diag(EinvH))
    noncentrality = totalN * hlt  
  }

  # get the critical value under the null
  Fcrit = qf(1 - glh@alpha, ndf, ddf)
  
  # get the power under the alternative
  return(1 - df(Fcrit, ndf, ddf, noncentrality))

}

#
# Calculate empirical power
#
glmmPower.empirical = function(design, glh, replicates=1000, blockSize=1000,
                               outputDir=".", filePrefix="simulatedData",
                               xNames=NA, yNames=NA, realizations=1000) {
  
  # generate data
  simulateData(design, replicates=replicates, blockSize=blockSize,
           outputDir=outputDir, filePrefix=filePrefix,
           xNames=xNames, yNames=yNames, realizations=realizations);
  
  # get placeholder for realization rejection counts
  powerList = rep(0, realizations)
  
  # fit models to data
  files = list.files(outputDir, pattern=cat(filePrefix,"*.csv", sep=""))
  ignore = sapply(files, function(file, powerList) {
    # open the data file
    data = read.csv(file)
    # create the model statement
    model = paste(c(
      paste(names(data[grepl("Y.*", names(data))]), sep=" + "),
      paste(c("0", names(data[grepl("XF.*|XG.*", names(data))])), collapse=" + ")
      ), collapse=" ~ "
    )
    
    ## select each realization
    tmp = sapply(unique(data$realizationID), 
                 function(realization, data, model, alpha, powerList) {
      realizationData = data[data$realizationID == realization,]
      
      ## select each data set within the realization
      numRejections = sum(sapply(unique(realizationData$setID), function(setID, data, model, alpha) {
        setData = data[data$setID == setID,]
        # fit the model
        fit.glmm = lm(model, setData)        
        # test contrasts 
        test = linearHypothesis(fit.glmm,
                                hypothesis.matrix=glh@betweenContrast,
                                test=glh@test,P=glh@withinContrast)
        test.results = getTestResults(test)
        return(test.results$Probf < 0.05)
      }, data, model, alpha)) 
      
      # increment the rejections for this realization
      powerList[realization] = powerList[realization] + powerList                    
    }, data, model, alpha, powerList)
   
  }, powerList=powerList)
  
  # calculate power as percent rejection
  powerList = powerList / 100
  
  # return the empirical unconditional power
  return(mean(powerList))  
  
}


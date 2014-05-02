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

# some utility functions used in JSON conversion
# and power calculations
# source("modelUtils.R")

#
# defineClasses.R
#
# Defines the following classes and related methods
#  design.glmmF - class describing study designs with fixed predictors only
#  design.glmmFG - class describing study designs with fixed and random predictors
#  glh - class describing a general linear hypothesis
#
#

########### CLASS DEFINITIONS ##########

#'
#' design.glmmF
#'
#' Class describing a general linear multivariate model design
#' in which all predictors are fixed and known as part of the study design. 
#'
#' @slot name An optional \code{character} string specifying the name of the study design.
#' @slot description An optional \code{character} string specifying the 
#'      brief description of the study design.
#' @slot perGroupN A \code{numeric} value indicating the number
#'      of participants in each study group.
#' @slot XEssence The design "essence" \code{matrix}.
#' @slot Beta The \code{matrix} of regression coefficients.
#' @slot SigmaError The residual covariance \code{matrix}.      
#'
#' @name design.glmmF 
#' @rdname design.glmmF
#' @note For theoretical details, please see
#' 
#' Muller, K. E., Lavange, L. M., Ramey, S. L., & Ramey, C. T. (1992). 
#' Power Calculations for General Linear Multivariate Models Including 
#' Repeated Measures Applications. Journal of the American Statistical 
#' Association, 87(420), 1209-1226.
#' 
setClass (
  "design.glmmF",
  representation ( name = "character",
                   description = "character",
                   XEssence = "matrix",
                   perGroupN = "numeric",
                   Beta = "matrix",
                   SigmaError = "matrix"
  ),
  prototype ( name ="",
              description ="",
              XEssence = diag(2),
              perGroupN = 10,
              Beta = matrix(c(1,0),nrow=2),
              SigmaError = matrix(c(1))
  ),
  validity = function(object) {
    #
    # Note, the class definition will already enforce that the
    # matrices are non-null 
    #
    
    # make sure the perGroupN is greater than 2
    if (object@perGroupN < 2) {
      stop("perGroupN is not an integer greater than or equal to 2")
    }
    
    # make sure that sigma is square, symmetric and positive definite
    if (!isSymmetric(object@SigmaError)) {
      stop("SigmaError matrix is not symmetric")
    } else {
      sigmaEigenValues = eigen(object@SigmaError)
      if (sum(sigmaEigenValues$values > 0) < length(sigmaEigenValues$values)) {
        stop("SigmaError matrix is not positive definite")
      }
    }
    
    # check matrix conformance
    if (ncol(object@XEssence) != nrow(object@Beta)) {
      stop("The number of columns in the essence X matrix does not match the 
           number of columns in the Beta matrix")
    } else if (ncol(object@Beta) != nrow(object@SigmaError)) {
      stop("The number of columns in the Beta matrix does not match the 
           number of rows in the SigmaError matrix")    
    } 
    
    return(TRUE)
  }
)


#'
#' design.glmmFG
#'
#' Class describing a general linear multivariate model design
#' containing both fixed predictors and one or more Gaussian covariates 
#'
#' @slot name An optional \code{character} string specifying the name of the study design.
#' @slot description An optional \code{character} string specifying the 
#'      brief description of the study design.
#' @slot perGroupN A \code{numeric} value indicating the number
#'      of participants in each study group.
#' @slot XEssence The design "essence" \code{matrix} for the fixed predictors.
#' @slot Beta The \code{matrix} of regression coefficients related to fixed predictors.
#' @slot SigmaY The covariance \code{matrix} for the outcomes, prior to adjusting for covariates.
#' @slot SigmaG The covariance \code{matrix} for the Gaussian predictors.
#' @slot SigmaYG The covariance \code{matrix} relating the outcomes and the Gaussian predictors.      
#'
#' @name design.glmmFG 
#' @rdname design.glmmFG
#' @note For theoretical details, please see
#' 
#' Glueck, D. H., & Muller, K. E. (2003). Adjusting power for a baseline covariate 
#' in linear models. Statistics in Medicine, 22(16), 2535-2551. doi:10.1002/sim.1341
#' 
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
#'
#' glh
#'
#' Class describing a general linear hypothesis for a multivariate model
#'
#' @slot alpha The Type I error rate.  Must be a \code{numeric} value between 0 and 1.
#' @slot betweenContrast The \code{matrix} of between participant contrasts.
#' @slot withinContrast The \code{matrix} of within participant contrasts.
#' @slot thetaNull The \code{matrix} of between participant contrasts.
#' @slot test A \code{character} string indicating the statistical test. At present,
#' only the Hotelling-Lawley trace is supported with value "Hotelling-Lawley".
#'
#' @name glh
#' @rdname glh
#' @note For theoretical details, please see
#' 
#' Muller, K. E., Lavange, L. M., Ramey, S. L., & Ramey, C. T. (1992). 
#' Power Calculations for General Linear Multivariate Models Including 
#' Repeated Measures Applications. Journal of the American Statistical 
#' Association, 87(420), 1209-1226.
#' 
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

########### END CLASS DEFINITIONS ##########

########### GENERIC DEFINITIONS ##########

#' simulateXMatrix
#'
#' Generate a complete X matrix for a given design.  If random
#' covariates are present in the design, the values of the covariates
#' will be randomly generated.
#'
#' @param design An object describing a study design, commonly a \code{design.glmmFG}
#'
#' @return A complete design matrix
#' 
#' @seealso \code{\link{design.glmmFG}}
#' 
#' @export
#' @docType methods
#' @rdname simulateXMatrix-methods
#' @name simulateXMatrix
#'
setGeneric("simulateXMatrix", function(design) standardGeneric("simulateXMatrix"))

#' simulateData
#'
#' Simulate data sets based on the specified design.
#'
#' @param design An object describing a study design, commonly a \code{design.glmmFG}
#' @param replicates the total number of data sets to generate
#' @param blockSize the data sets can be written to multiple files.  The \code{blockSize}
#' is the number of data sets which will be contained in each file.
#' @param outputDir the directory to which the data set files will be written
#' @param filePrefix the filename prefix for each data set file
#' @param xNames the column names for the predictors in the data set
#' @param yNames the column names for the outcomes in the data set 
#'
#' @return Writes one or more CSV files containing data sets 
#' to the output directory.
#' 
#' @seealso \code{\link{design.glmmFG}} and \code{\link{design.glmmF}} 
#' 
#' @export
#' @docType methods
#' @rdname simulateData-methods
#' @name simulateData
#'
setGeneric("simulateData", function(design, replicates=10000, blockSize=1000,
                                    outputDir=".", filePrefix="simulatedData",
                                    xNames=NA, yNames=NA, ...) 
  standardGeneric("simulateData"))


#' designToJSON
#'
#' Generate a character string with the JSON representation of a
#' study design object, optionally with a complete X matrix.  
#' Used to pass complex R objects into the
#' supporting Java library for the \code{rPowerlib} package.
#'
#' @param design An object describing a study design
#' @param expandX If true, the JSON representation will include a complete X matrix.
#' 
#' @return The JSON representation of the design object
#' 
#' @seealso \code{\link{design.glmmFG}} and \code{\link{design.glmmF}} 
#' 
#' @export
#' @docType methods
#' @rdname designToJSON-methods
#' @name designToJSON
#'
setGeneric("designToJSON", function(obj, expandX=FALSE) standardGeneric("designToJSON"))

#' toJSON
#'
#' Generate a character string with the JSON representation of an object,
#' most commonly a \code{glh} object. Used to pass complex R objects into the
#' supporting Java library for the \code{rPowerlib} package.
#'
#' @param design An object describing a study design
#' @param expandX If true, the JSON representation will include a complete X matrix.
#' 
#' @return The JSON representation of the object
#' 
#' @seealso \code{\link{glh}} 
#' 
#' @export
#' @docType methods
#' @rdname toJSON-methods
#' @name toJSON
#'
setGeneric("toJSON", function(obj) standardGeneric("toJSON"))

#' empiricalPower
#'
#' Calculate empirical unconditional power for a given design and hypothesis.  This
#' function is implemented purely in R and may run slowly for complex
#' designs.
#'
#' @param design An object describing a study design
#' @param glh An object describing the general linear hypothesis to be tested
#' @replicates total number of error replicates to produce for a given X matrix
#' @realizations total number of realizations of the X matrix to generate
#' 
#' @return empirical unconditional power
#' 
#' @seealso \code{\link{design.glmmFG}} and \code{\link{glh}} 
#' 
#' @export
#' @docType methods
#' @rdname empiricalPower-methods
#' @name empiricalPower
#'
setGeneric("empiricalPower", 
           function(design, glh, replicates=1000, realizations=1000) standardGeneric("empiricalPower"))

#' fastEmpiricalPower
#'
#' Calculate empiricalPower for a given design and hypothesis.  Utilizes
#' Java code to provide faster computation of empirical power.
#'
#' @param design An object describing a study design
#' @param glh An object describing the general linear hypothesis to be tested
#' @replicates total number of error replicates to produce for a given X matrix
#' @realizations total number of realizations of the X matrix to generate
#' 
#' @return empirical unconditional power
#' 
#' @seealso \code{\link{design.glmmFG}} and \code{\link{glh}} 
#' 
#' @export
#' @docType methods
#' @rdname fastEmpiricalPower-methods
#' @name fastEmpiricalPower
#'
setGeneric("fastEmpiricalPower", 
           function(design, glh, realizations=1000, replicates=1000) standardGeneric("fastEmpiricalPower"))



########### END GENERIC DEFINITIONS ##########

########### METHOD DEFINITIONS ##########

#' @rdname simulateData-methods
#' @aliases simulateData,design.glmmF,design.glmmF-method
setMethod("simulateData", "design.glmmF", 
          function(design, replicates=10000, blockSize=1000,
                   outputDir=".", filePrefix="simulatedData",
                   xNames=NA, yNames=NA) {
            if (is.na(design@perGroupN)) {
              stop("Per group sample size not specified in study design")
            }
            if (!is.na(xNames) && length(xNames) != ncol(design@XEssence)) {
              stop("The number of values in xNames does not match the number of columns in XEssence")
            }
            if (!is.na(yNames) && length(yNames) != ncol(design@Beta)) {
              stop("The number of values in yNames does not match the number of columns in Beta")
            }
            
            ### calculate the mean, XB ###
            # first, calculate the full X matrix
            X = matrix(rep(1,design@perGroupN), nrow=design@perGroupN) %x% design@XEssence
            # calculate XB
            XB = X %*% design@Beta
            # total sample size
            totalN = nrow(X) 
            
            ### determine the number of sets and the size of each set ###
            
            # when the blocksize does not divide evenly into the
            # total replicates, the size of the last set is the remainder
            lastSetSize = replicates %% blockSize;
            numSets = floor(replicates / blockSize)
            if (lastSetSize > 0) {
              numSets = numSets + 1
            }
            
            ### generate column names for the data sets ###
            yPredef = sapply(1:ncol(XB), function(x) {paste("Y.",x,sep="",collapse="")})
            xPredef = sapply(1:ncol(X), function(x) {paste("X.",x,sep="",collapse="")})
            if (!is.na(yNames) && !is.na(xNames)) {
              dataSetNames = c("setID", yNames, xNames)
            } else if (!is.na(yNames) && is.na(xNames)) {
              dataSetNames = c("setID", yNames, xPredef)
            } else if (is.na(yNames) && !is.na(xNames)) {
              dataSetNames = c("setID", yPredef, xNames)
            } else {
              dataSetNames = c("setID", yPredef, xPredef)
            }
            
            ### buil the data set blocks ###
            sapply(1:numSets, function(setNumber) {
              # set start and end iteration numbers for this block
              startIter = ((setNumber-1)*blockSize) + 1
              if (setNumber == numSets && lastSetSize > 0 && lastSetSize != blockSize) {
                numDataSets = lastSetSize
              } else {
                numDataSets = blockSize
              }
              endIter = startIter + numDataSets - 1;
              cat("Generating data sets ", startIter, " to ", endIter, "\n")
              
              # generate numDataSets number of data sets
              dataSet = do.call("rbind", 
                                lapply(1:numDataSets, function(blockNumber) {
                                  errorMatrix = 
                                    mvrnorm(n = totalN, 
                                            mu=rep(0, nrow(design@SigmaError)), 
                                            design@SigmaError)
                                  yData = XB + errorMatrix
                                  dataBlock = data.frame(setID=((setNumber-1)*blockSize+blockNumber),Y=yData, X=X)
                                  return(dataBlock)
                                  
                                }))
              # write to the output directory
              names(dataSet) = dataSetNames
              filename = paste(outputDir,"/", filePrefix,startIter,"to",endIter,".csv",sep="")
              write.csv(dataSet, row.names=FALSE, file=filename)
              return(filename)
            }) # end sapply
          }
)

#' @rdname designToJSON-methods
#' @aliases designToJSON,design.glmmF,design.glmmF-method
setMethod("designToJSON", signature("design.glmmF"), 
          function(obj, expandX=FALSE) {
            if (expandX) {
              XFull = obj@XEssence %x% matrix(rep(1,obj@perGroupN))
              xJSON = matrixToJSON("XFixed", XFull)

            } else {
              xJSON = paste(
                c(matrixToJSON("XEssence", obj@XEssence), 
                  paste(c("\"perGroupN\": ", obj@perGroupN), collapse="")    
                ), collapse=",")
            }
            
            return(
              paste(c(
                "{\"name\": \"", obj@name, "\",",
                "\"description\": \"", obj@description, "\",",
                xJSON, ",",
                matrixToJSON("beta", obj@Beta), ",",
                matrixToJSON("sigmaError", obj@SigmaError), "}"
              ),
                    collapse="")
            )
          })





#' @rdname simulateData-methods
#' @aliases simulateData,design.glmmFG,design.glmmFG-method
setMethod("simulateData", "design.glmmFG", 
          function(design, replicates=1000, blockSize=100,
                   outputDir=".", filePrefix="simulatedData",
                   xNames=NA, yNames=NA, realizations=1000) {
            if (is.na(design@perGroupN)) {
              stop("Per group sample size not specified in study design")
            }
            if (!is.na(xNames) && length(xNames) != (ncol(design@XEssence) + ncol(design@SigmaG))) {
              stop("The number of values in xNames does not match the number of columns in XEssence")
            }
            if (!is.na(yNames) && length(yNames) != ncol(design@Beta)) {
              stop("The number of values in yNames does not match the number of columns in Beta")
            }
            
            ### calculate the mean, XB for fixed predictors ###
            # first, calculate the full X matrix
            X = matrix(rep(1,design@perGroupN), nrow=design@perGroupN) %x% design@XEssence
            # calculate XB
            XB = X %*% design@Beta
            # total sample size
            totalN = nrow(X) 
            
            ### determine the number of sets and the size of each set ###
            
            # when the blocksize does not divide evenly into the
            # total replicates, the size of the last set is the remainder
            lastSetSize = replicates %% blockSize;
            numSets = floor(replicates / blockSize)
            if (lastSetSize > 0) {
              numSets = numSets + 1
            }
            
            ### generate column names for the data sets ###
            yPredef = sapply(1:ncol(XB), function(x) {paste("Y.",x,sep="",collapse="")})
            xPredef = c(sapply(1:ncol(X), function(x) {paste("XF.",x,sep="",collapse="")}),
                        sapply(1:ncol(design@SigmaG), function(x) {paste("XG.",x,sep="",collapse="")}))
            if (!is.na(yNames) && !is.na(xNames)) {
              dataSetNames = c("realizationID", "setID", yNames, xNames)
            } else if (!is.na(yNames) && is.na(xNames)) {
              dataSetNames = c("realizationID", "setID", yNames, xPredef)
            } else if (is.na(yNames) && !is.na(xNames)) {
              dataSetNames = c("realizationID", "setID", yPredef, xNames)
            } else {
              dataSetNames = c("realizationID", "setID", yPredef, xPredef)
            }
            
            #
            # Calculate the covariance of errors per Glueck and Muller
            #  
            SigmaError = design@SigmaY - design@SigmaYG %*% solve(design@SigmaG) %*% t(design@SigmaYG)
            
            ### generate the realizations ###
            sapply(1:realizations, function(realizationID) {
              
              XG = mvrnorm(n = totalN, 
                           mu=rep(0, nrow(design@SigmaG)), 
                           design@SigmaG)
              ### buil the data set blocks ###
              sapply(1:numSets, function(setNumber) {
                # set start and end iteration numbers for this block
                startIter = ((setNumber-1)*blockSize) + 1
                if (setNumber == numSets && lastSetSize > 0 && lastSetSize != blockSize) {
                  numDataSets = lastSetSize
                } else {
                  numDataSets = blockSize
                }
                endIter = startIter + numDataSets - 1;
                cat("Generating data sets ", startIter, " to ", endIter, "\n")
                
                # generate numDataSets number of data sets
                dataSet = do.call("rbind", 
                                  lapply(1:numDataSets, function(blockNumber) {
                                    errorMatrix = 
                                      mvrnorm(n = totalN, 
                                              mu=rep(0, nrow(SigmaError)), 
                                              SigmaError)
                                    yData = XB + errorMatrix
                                    dataBlock = data.frame(realizationID=realizationID,
                                                           setID=((setNumber-1)*blockSize+blockNumber),
                                                           Y=yData, X=X, XG=XG)
                                    return(dataBlock)
                                    
                                  }))
                # write to the output directory
                names(dataSet) = dataSetNames
                filename = paste(outputDir,"/", filePrefix,"Realization", realizationID, "Iter", 
                                 startIter,"to",endIter,".csv",sep="")
                write.csv(dataSet, row.names=FALSE, file=filename)
                return(filename)
              }) # end sapply
            })
          }
)



#' @rdname simulateXMatrix-methods
#' @aliases simulateXMatrix,design.glmmFG,design.glmmFG-method
setMethod("simulateXMatrix", "design.glmmFG", 
          function(design) {
            # first, calculate the full X matrix
            XF = matrix(rep(1,design@perGroupN), nrow=design@perGroupN) %x% design@XEssence
            # total sample size
            totalN = nrow(XF)
            # generate the covariates
            XG = mvrnorm(n = totalN, 
                         mu=rep(0, nrow(design@SigmaG)),
                         design@SigmaG)
            
            X = data.frame(XF,XG)
            names(X) = c(sapply(1:ncol(XF), function(x) {paste("XF.",x,sep="",collapse="")}),
                         sapply(1:ncol(design@SigmaG), function(x) {paste("XG.",x,sep="",collapse="")}))
            
            return(X)
          }
)

#' @rdname designToJSON-methods
#' @aliases designToJSON,design.glmmFG,design.glmmFG-method
setMethod("designToJSON", signature("design.glmmFG"), 
          function(obj, expandX) {
            if (expandX) {
              XFull = obj@XEssence %x% matrix(rep(1,obj@perGroupN))
              xJSON = matrixToJSON("XFixed", XFull)
            } else {
              xJSON = paste(
                c(matrixToJSON("XEssence", obj@XEssence), 
                  paste(c("\"perGroupN\": ", obj@perGroupN), collapse="")),
                collapse=","      
              )
            }
            
            return(
              paste(c(
                "{\"name\": \"", obj@name, "\",",
                "\"description\": \"", obj@description, "\",",
                xJSON, ",",
                matrixToJSON("betaFixed", obj@Beta), ",",
                matrixToJSON("sigmaY", obj@SigmaY), ",",
                matrixToJSON("sigmaG", obj@SigmaG), ",",
                matrixToJSON("sigmaYG", obj@SigmaYG), "}"
              ),
                    collapse="")
              
            )
          })


#' @rdname empiricalPower-methods
#' @aliases empiricalPower,design.glmmFG,design.glmmFG-method
setMethod("empiricalPower", "design.glmmFG", 
          function(design, glh, replicates=1000, realizations=1000) {
            
            # Calculate the covariance of errors   
            invSigmaG = solve(design@SigmaG)
            SigmaError = design@SigmaY - design@SigmaYG %*% invSigmaG %*% t(design@SigmaYG)
            # calculate the random portion of the beta matrix
            BetaRandom = invSigmaG %*% t(design@SigmaYG)     
            # form the complete Beta matrix
            BetaFull = rbind(design@Beta, BetaRandom)
            
            # generate realizations of the X matrix
            powerValues = sapply(1:realizations, 
                                 function(realizationID, replicates, Beta, SigmaError) {
                                   if (realizationID %% 100 == 0) {
                                     print(paste(c("Calculating empirical power for realization", realizationID), 
                                                 collapse=" "))
                                   }
                                   
                                   # get X matrix realization
                                   XMatrix = as.matrix(simulateXMatrix(design))
                                   XtXInverse = solve(t(XMatrix) %*% XMatrix)
                                   rank = qr(XMatrix)$rank
                                   
                                   # calculate empirical power for this realization of X
                                   rejectionList = sapply(1:replicates, 
                                                          function(replicate, XtXInverse, XMatrix, rankX,
                                                                   Beta, SigmaError) {
                                                            totalN = nrow(XMatrix)
                                                            
                                                            # generate a data set replicate
                                                            errorMatrix = mvrnorm(n = totalN, mu=rep(0, nrow(SigmaError)), SigmaError)
                                                            yData = XMatrix %*% Beta + errorMatrix
                                                            
                                                            # fit the model and test the hypothesis
                                                            results <- fitModelAndTest(XtXInverse, XMatrix, yData, glh, totalN, rankX)
                                                            pvalue = results$ProbF
                                                            
                                                            # return 1 if we rejected the null, 0 otherwise
                                                            return(as.numeric(pvalue <= glh@alpha))
                                                            
                                                          }, XtXInverse=XtXInverse, XMatrix=XMatrix, rankX=rank, 
                                                          Beta=BetaFull, SigmaError=SigmaError)
                                   
                                   # empirical power for a given realization is the percentage of rejections
                                   return(sum(rejectionList)/replicates)
                                   
                                 }, replicates=replicates, Beta=BetaFull, SigmaError=SigmaError)
            # print(powerValues)
            
            # obtain empirical unconditional power by averaging across the realizations
            # of X
            return(mean(powerValues))
          }
)

#' @rdname fastEmpiricalPower-methods
#' @aliases fastEmpiricalPower,design.glmmFG,design.glmmFG-method
setMethod("fastEmpiricalPower", "design.glmmFG", 
          function(design, glh, realizations=1000, replicates=1000) {
            
            obj=.jnew("com/kreidles/PowerCalculator")
            power = .jcall(obj, "D", "calculateEmpiricalPower", 
                           designToJSON(design, TRUE), toJSON(glh), 
                           as.integer(realizations), as.integer(replicates))
            return(power)
          }
)

#' @rdname toJSON-methods
#' @aliases toJSON,glh,glh-method
setMethod("toJSON", "glh", 
          function(obj) {            
            return(
              paste(c(
                "{\"alpha\": ", obj@alpha, ",",
                "\"test\": \"", obj@test, "\",",
                matrixToJSON("betweenContrast", obj@betweenContrast), ",",
                matrixToJSON("withinContrast", obj@withinContrast), ",",
                matrixToJSON("thetaNull", obj@thetaNull), "}"
              ),
                    collapse="")
              
            )
            
          })

########### END METHOD DEFINITIONS ##########

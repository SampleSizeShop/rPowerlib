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
# Provides the design.glmmFG class which defines matrices
# for the general linear multivariate or univariate model
# with one or more Gaussian covariates.
#
# For notation and theoretical details, see
# 
# 1. Glueck DH, Muller KE. Adjusting power for a baseline covariate in linear models. 
# Stat Med. 2003;22(16):2535–2551. doi:10.1002/sim.1341.
# 
# 2. Muller KE, Stewart PW. Linear model theory: univariate, multivariate, and mixed models. 
# Hoboken, New Jersey: John Wiley and Sons; 2006.
#
#
# 
# library(magic)
# library(car)
# library(MASS)
# library(rJava)
# 
# classpath = paste(c(getwd(), "/inst/java/com.kreidles.covariatepowersimulation-1.0.0.jar"), collapse="")
# .jinit(classpath=classpath)
# source("modelUtils.R")

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

##
# simulateData: simulate data sets for the specified study design
#
# Args:
#  design (required) - the glmmFG study design object
#  replicates (optional) - the total number of data sets
#  blockSize (optional) - the number of data sets to include in each file
#  outputDir (required) - directory in which to write the data sets
#  filePrefix (optional) - prefix added to filenames
#  xNames (optional) - column names for the predictors
#  yNames (optional) - column names for the outcomes
#
# Outputs:
#  writes simulated data to disk in CSV format.  
#  Filnames follow the pattern <filePrefix><start>To<end>.csv 
#  (if filePrefix unspecified, then simulatedData<start>To<end>.csv)
#  Files have the format
# 
#  dataSetID | Y (outcomes) | X (predictors)
#
#
setGeneric("simulateData", function(design, replicates=10000, blockSize=1000,
                                    outputDir=".", filePrefix="simulatedData",
                                    xNames=NA, yNames=NA, ...) standardGeneric("simulateData"))
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



#
# Obtain a single realization of the X matrix
# with randomly generated covariates
#
#
setGeneric("simulateXMatrix", function(design) standardGeneric("simulateXMatrix"))
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

setGeneric("toJSON.XEssence", function(obj) standardGeneric("toJSON.XEssence"))
setMethod("toJSON.XEssence", "design.glmmFG", 
          function(obj) {
            xJSON = paste(
              c(matrixToJSON("XEssence", obj@XEssence), 
                paste(c("\"perGroupN\": ", obj@perGroupN), collapse="")),
                collapse=","      
              )
            
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

setGeneric("toJSON", function(obj) standardGeneric("toJSON"))
setMethod("toJSON", "design.glmmFG", 
          function(obj) {
            XFull = obj@XEssence %x% matrix(rep(1,obj@perGroupN))
            xJSON = matrixToJSON("XFixed", XFull)
            
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
          } )

#
# Calculate empirical power for the specified GLMM(F,G) and
# corresponding hypothesis
#
setGeneric("empiricalPower", 
           function(design, glh, replicates, realizations) standardGeneric("empiricalPower"))
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

setGeneric("fastEmpiricalPower", 
           function(design, glh, realizations=1000, replicates=1000) standardGeneric("fastEmpiricalPower"))
setMethod("fastEmpiricalPower", "design.glmmFG", 
          function(design, glh, realizations=1000, replicates=1000) {
            
            obj=.jnew("com/kreidles/PowerCalculator")
            power = .jcall(obj, "D", "calculateEmpiricalPower", 
                           toJSON(design), toJSON(glh), 
                           as.integer(realizations), as.integer(replicates))
            return(power)
          }
)

#
# Example multivariate
#
design.cov3 = new("design.glmmFG", 
             Beta=matrix(c(1,0,0,0,0,0), nrow=2), 
             SigmaY=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1), nrow=3), 
             SigmaG=matrix(c(1,0.01,0.01,0.01,1,0.01,0.01,0.01,1), nrow=3), 
             SigmaYG=matrix(rep(c(0.2,0.4,0.5), 3), nrow=3, byrow=TRUE))
design.cov1 = new("design.glmmFG", 
             Beta=matrix(c(1,0,0,0,0,0), nrow=2), 
             SigmaY=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1), nrow=3), 
             SigmaG=matrix(c(1.2), nrow=1), 
             SigmaYG=matrix(rep(c(0.2,0.4,0.5), 1), nrow=3, byrow=TRUE))
# startTime <- proc.time()
#  fastEmpiricalPower(design, glh, 1000,1000)
#  ellapsed = proc.time() - startTime
# startTime <- proc.time()
# empiricalPower(design, glh, replicates=1000, realizations=1000)
# ellapsed = proc.time() - startTime
# 
# startTime <- proc.time()
# simulateData(design, outputDir="../data/", blockSize=1000)
# ellapsed = proc.time() - startTime
# 
# Rprof("prof.out")
# prof.out = profr(empiricalPower(design, glh, replicates=1000, realizations=10), quiet=FALSE)
# Rprof(NULL)
# ellapsed = proc.time() - startTime



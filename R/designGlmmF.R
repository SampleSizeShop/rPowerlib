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
# For notation an theoretical details, see
#
# 1. Muller KE, Stewart PW. Linear model theory: univariate, multivariate, and mixed models. 
# Hoboken, New Jersey: John Wiley and Sons; 2006.
#
#
#
library(MASS)

#
# design.glmmF
#
# Class describing the general linear model with
# only fixed predictors
#
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


# test = new("design.glmmF", name="Pre/Post", 
#              SigmaError=diag(2), XEssence=diag(2), Beta=matrix(c(1,1,0,0), nrow=2)
#            )

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

#
# Generate the JSON representation of the design
#
setMethod("toJSON", "design.glmmF", 
          function(obj) {
            xJSON = paste(
              c(matrixToJSON("XEssence", obj@XEssence), 
                paste(c("\"perGroupN\": ", obj@perGroupN), collapse="")    
              ), collapse=",")
            
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


# test design
designF = new("design.glmmF", 
             Beta=matrix(c(1,0,0,0,0,0), nrow=2), 
             SigmaError=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1), nrow=3))


#
# test power in java stats
#



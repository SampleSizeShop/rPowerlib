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

source("modelUtils.R")

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
# Calculate the hypothesis sum of squares for the given glh, beta, and sigma matrices
#
setGeneric("getHypothesisSumOfSquares", function(glh, beta) standardGeneric("getHypothesisSumOfSquares"))
setMethod("getHypothesisSumOfSquares", "glh", 
          function(glh, beta) {
          }
)

#
# Calculate the error sum of squares for a given GLH and error matrix
#
setGeneric("getErrorSumOfSquares", function(glh, sigmaE) standardGeneric("getErrorSumOfSquares"))
setMethod("getErrorSumOfSquares", "glh",
          function(glh, sigmaE) {
            return(t(glh@withinContrast) %*% sigmaE %*% glh@withinContrast)            
          }
)


#setGeneric("toJSON", function(obj) standardGeneric("toJSON"))
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

glh.cov3 = new("glh", betweenContrast=matrix(c(1,-1,0,0,0),nrow=1), 
          withinContrast=matrix(c(1,1,-1,0,0,-1), nrow=3,byrow=TRUE), 
          thetaNull=matrix(c(0,0), nrow=1))

glh.cov1 = new("glh", betweenContrast=matrix(c(1,-1,0),nrow=1), 
          withinContrast=matrix(c(1,1,-1,0,0,-1), nrow=3,byrow=TRUE), 
          thetaNull=matrix(c(0,0), nrow=1))

glhF = new("glh", betweenContrast=matrix(c(1,-1),nrow=1), 
          withinContrast=matrix(c(1,1,-1,0,0,-1), nrow=3,byrow=TRUE), 
          thetaNull=matrix(c(0,0), nrow=1))


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
# Define classes for the following types of study designs
# 1. General linear univariate or multivariate model with fixed predictors
# 2. General linear unviariate or multivariate model with fixed predictor
#    plus one or more Gaussian covariates
#
# For a description of 
#


#
# design.glmmF
#
# Class describing the general linear model with
# only fixed predictors
#
setClass ("design.glmmF",
          representation ( name = "character",
                           description = "character",
                           XEssence = "matrix",
                           Beta = "matrix",
                           SigmaError = "matrix",
                           C = "matrix",
                           U = "matrix",
                           thetaNull = "matrix"
          ),
          prototype ( name ="two sample t-test",
                      description ="A design which tests for the difference of 
                      means between two independent groups",
                      X = diag(2),
                      beta = matrix(c(1,0),nrow=2),
                      Sigma = matrix(c(1)),
                      C = matrix(c(1,-1), nrow=1),
                      thetaNull = matrix(c(0))
          )
          )

## validation routines ##

setGeneric("validObject", function(object) standardGeneric("validObject"))
setMethod("validObject", "mixedStudyDesign", function(object){
  # TODO
  return(TRUE);
})

## Getters / Setters ##
setGeneric("name", function(object) standardGeneric("name"))
setMethod("name", "mixedStudyDesign", function(object){return(object@name)})
setGeneric("name<-", function(object, value) standardGeneric("name<-"))
setReplaceMethod("name", "mixedStudyDesign", function(object, value){
  object@name <- value
  validObject(object)
  object
})

setGeneric("description", function(object) standardGeneric("description"))
setMethod("description", "mixedStudyDesign", function(object){return(object@description)})
setGeneric("description<-", function(object, value) standardGeneric("description<-"))
setReplaceMethod("description", "mixedStudyDesign", function(object, value){
  object@description <- value
  validObject(object)
  object
})

setGeneric("XMatrix", function(object) standardGeneric("XMatrix"))
setMethod("XMatrix", "mixedStudyDesign", function(object){return(object@X)})
setGeneric("XMatrix<-", function(object, value) standardGeneric("XMatrix<-"))
setReplaceMethod("XMatrix", "mixedStudyDesign", function(object, value){
  object@X <- value
  validObject(object)
  object
})

setGeneric("betaMatrix", function(object) standardGeneric("betaMatrix"))
setMethod("betaMatrix", "mixedStudyDesign", function(object){return(object@beta)})
setGeneric("betaMatrix<-", function(object, value) standardGeneric("betaMatrix<-"))
setReplaceMethod("betaMatrix", "mixedStudyDesign", function(object, value){
  object@beta <- value
  validObject(object)
  object
})

setGeneric("SigmaMatrix", function(object) standardGeneric("SigmaMatrix"))
setMethod("SigmaMatrix", "mixedStudyDesign", function(object){return(object@Sigma)})
setGeneric("SigmaMatrix<-", function(object, value) standardGeneric("SigmaMatrix<-"))
setReplaceMethod("SigmaMatrix", "mixedStudyDesign", function(object, value){
  object@Sigma <- value
  validObject(object)
  object
})

setGeneric("contrast", function(object) standardGeneric("contrast"))
setMethod("contrast", "mixedStudyDesign", function(object){return(object@C)})
setGeneric("contrast<-", function(object, value) standardGeneric("contrast<-"))
setReplaceMethod("contrast", "mixedStudyDesign", function(object, value){
  object@C <- value
  validObject(object)
  object
})


setGeneric("thetaNullMatrix", function(object) standardGeneric("thetaNullMatrix"))
setMethod("thetaNullMatrix", "mixedStudyDesign", function(object){return(object@thetaNull)})
setGeneric("thetaNullMatrix<-", function(object, value) standardGeneric("thetaNullMatrix<-"))
setReplaceMethod("thetaNullMatrix", "mixedStudyDesign", function(object, value){
  object@thetaNull <- value
  validObject(object)
  object
})


setGeneric("sascall", function(object) standardGeneric("sascall"))
setMethod("sascall", "mixedStudyDesign", function(object){return(object@sascall)})
setGeneric("sascall<-", function(object, value) standardGeneric("sascall<-"))
setReplaceMethod("sascall", "mixedStudyDesign", function(object, value){
  object@sascall <- value
  validObject(object)
  object
})

test = new("mixedStudyDesign")



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
# Provides the criteria.glh class which defines a general linear
# hypothesis used to calculate power or sample size
#
# For notation and theoretical details, see
# 
# 1. Muller KE, Stewart PW. Linear model theory: univariate, multivariate, and mixed models. 
# Hoboken, New Jersey: John Wiley and Sons; 2006.
#
#

#
# criteria.glh
#
# Class describing a general linear
# hypothesis used to calculate power or sample size
#
setClass (
  "criteria.glh",
  representation ( alpha = "numeric",
                   nominalPower = "numeric",
                   C = "matrix",
                   U = "matrix",
                   ThetaNull = "matrix"
  ),
  prototype ( alpha = 0.05,
              nominalPower = NA,
              C = matrix(c(1,-1), nrow=1),
              U = diag(1),
              ThetaNull = matrix(c(0))
  ),
  validity = function(object) {
    # make sure alpha is valid
    if (object@alpha < 0 | object@alpha > 1) {
      stop("alpha must be between 0 and 1")
    }
    # check matrix conformance
    if (nrow(object@ThetaNull) != nrow(object@C)) {
      stop("The number of rows in the ThetaNull matrix does not match the
           number of rows in the C contrast matrix")
    } else if (ncol(object@ThetaNull) != ncol(object@U)) {
      stop("The number of columns in the ThetaNull matrix does not match the
           number of columns in the U contrast matrix")
    }
    
    return(TRUE)
  }
)
# 
# else if (ncol(object@C) != nrow(object@Beta)) {
#   stop("The number of columns in the C contrast matrix does not match the 
#            number of rows the Beta matrix")
# } else if (ncol(object@Beta) != nrow(object@U)) {
#   stop("The number of columns in the Beta matrix does not match the 
#            number of rows in the U contrast matrix")
# } else if (nrow(object@U) != nrow(object@SigmaError)) {
#   stop("The number of rows in the U contrast matrix does not match the
#            number of rows in the SigmaError matrix")
# } else if (nrow(object@ThetaNull) != nrow(object@C)) {
#   stop("The number of rows in the ThetaNull matrix does not match the
#            number of rows in the C contrast matrix")
# } else if (ncol(object@ThetaNull) != ncol(object@U)) {
#   stop("The number of columns in the ThetaNull matrix does not match the
#            number of columns in the U contrast matrix")
# }





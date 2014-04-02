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

#
# This script runs power for the applied power example in 
# Section 5 of the manuscript 
# "Calculating Power for the General Linear Multivariate
#  Model With One or More Gaussian Covariates"
#  Authors: Sarah M. Kreidler, Keith E. Muller, and Deborah H. Glueck
#

#
# Applied example based on previous research by
#
# 1. Elashoff D, Zhou H, Reiss J, et al. 
# Prevalidation of salivary biomarkers for oral cancer detection. 
# Cancer Epidemiol Biomarkers Prev. 2012;21(4):664–672. 
#
# The proposed design is a randomized, clinical trial investigating
# two treatments for oral cancer.  The outcome is an oral cancer 
# salivary biomarker, IL-8, which would be measured at baseline,
# six months, and 1 year.  Covariates include age and baseline IL-8 level.
# The hypothesis is the test of the time by treatment interaction, using
# the Hotelling-Lawley Trace with a Type I error rate of 0.05.
#

# Mean differences from Table 3 in Elashoff et al.
il8Data = data.frame(control=c(33.8,30.9,19.6,19.8,18.6),
                     oscc=c(32.4,27.7,17.4,17.7,16.8))
il8Data$diff = il8Data$control - il8Data$oscc
mean(il8Data$diff)

# Covariance from Tables 1 and 3 in Elashoff et al.
ageSD = mean(12.7,12.1,13.4,15.2,11.4,8,9.1,13.5,10,11)
ageSD^2

il8SD = mean(2.2,3.6,2.5,3.9,3.3,2.8,2.9,2.6,4.0,2.2)
il8SD^2

# create the design object for the proposed trial
design.elashoff = new("design.glmmFG", 
                  XEssence=diag(2),
                  perGroupN=35,
                  Beta=matrix(c(2.14,2.14,0,2.14), nrow=2), 
                  SigmaY=4.84*matrix(c(1,0.4,0.4,1), nrow=2), 
                  SigmaG=diag(c(4.84,161.29)), 
                  SigmaYG=matrix(c(1.94,2.79,0.97,2.79), nrow=2, byrow=TRUE))
# create the hypothesis object for the proposed tria
glh.elashoff = new("glh", betweenContrast=matrix(c(1,-1,0,0),nrow=1), 
               withinContrast=matrix(c(1,-1), nrow=2), 
               thetaNull=matrix(c(0)))

# for demonstration purposes, we show the final residual error
SigmaError = design.elashoff@SigmaY - 
  design.elashoff@SigmaYG %*% 
  solve(design.elashoff@SigmaG) %*% t(design.elashoff@SigmaYG)
# this calculation is only for demonstration purposes
# and is not required for the power calculation
SigmaError

# show the magnitude of the interaction effect.
# this calculation is only for demonstration purposes
# and is not required for the power calculation
interactionSize = (glh.elashoff@betweenContrast[1:2] %*% 
                     design.elashoff@Beta %*%
                     glh.elashoff@withinContrast)
  
# calculate power for this design
power = glmmPower.covariateAdjusted(design.elashoff, glh.elashoff)
info = paste(c("Power is ", round(power, 3), 
               " for the Hotelling-Lawley Trace test of time by treatment interaction, ",
               "with an interaction magnitude of ", interactionSize, ".\n"), collapse="")

cat(info)




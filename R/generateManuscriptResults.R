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
# generateManuscriptResults.R
#
# Generate power results, summary tables, and plots for manuscript
# !!! Be sure to set working directory to this file location prior to running !!!

# load libraries and initialize Java
generateManuscriptResults <- function(genDesigns=TRUE, runEmpirical=FALSE) {
  
}

# define functions

#
#
#



#
# generate the designs for the experiment or load from
# previously generated .Rdata file
#

#
# Calculate empirical power or load csv from previous run
#
#
# Generate a data set with empirical power for each design
#
#
empiricalPowerData = data.frame(
  designName=sapply(manuscriptDesignList, function(x) { return(x[[1]]@name)}),
  perGroupN=sapply(manuscriptDesignList, function(x) { return(x[[3]]['perGroupN'])}),
  sigmaYGscale=sapply(manuscriptDesignList, function(x) { return(x[[3]]['sigmaYGscale'])}),
  sigmaGscale=sapply(manuscriptDesignList, function(x) { return(x[[3]]['sigmaGscale'])})
)
# calculate empirical power for each design
empiricalPowerAndTimeList = lapply(manuscriptDesignList, function(x) {
  print(paste(c("Calculating power for '", x[[1]]@name ,
                "', N=", x[[3]]['perGroupN'], 
                ", SigmaYGscale=", x[[3]]['sigmaYGscale'],
                ", SigmaGscale=", x[[3]]['sigmaGscale']),
              collapse=""))
  startTime <- proc.time()
  power = fastEmpiricalPower(x[[1]], x[[2]])
  ellapsed = proc.time() - startTime
  print(paste(c("Done (", ellapsed[[1]], "s)"), collapse=""))
  return(list(power, ellapsed))
})

## add the timing results and the empirical power values to the data
empiricalPowerData = data.frame(empiricalPowerData,
                                empiricalPower=sapply(empiricalPowerAndTimeList, function(x) {return(x[[1]])}),
                                time=sapply(empiricalPowerAndTimeList, function(x) {return(x[[2]][1])})

#
# Add columns for calculated power using the following methods:
# 1. Kreidler covariate adjust with M adjusted by (N-qf)/(N-qf-qr)
# 2. Kreidler covariate adjust with M
#



#
# Write the power results to a csv file
#


#
# Generate the box plots showing deviations from
# empirical power
#



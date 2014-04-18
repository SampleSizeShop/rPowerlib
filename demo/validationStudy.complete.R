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
# This script runs the complete validation study for the manuscript
# "Calculating Power for the General Linear Multivariate Model 
# With One or More Gaussian Covariates" by Sarah M. Kreidler,
# Keith E. Muller, and Deborah H. Glueck
#
# The simulation study evaluates the performance of a power approximation
# for designs with Gaussian covariates.
#
# Please note that the "single covariate" power calculations take 15-30 minutes to run.
# Also, the empirical power calculations require several hours to run (about 5-8 hrs).
#

# generate the designs for the validation study
designList = generateDesignsForManuscript()

# build a data frame with all combinations of parameters
# which are varied for each design.
powerData = data.frame(
  designName=sapply(designList, function(x) { return(x[[1]]@name)}),
  perGroupN=sapply(designList, function(x) { return(x[[3]]['perGroupN'])}),
  sigmaYGscale=sapply(designList, function(x) { return(x[[3]]['sigmaYGscale'])}),
  betaScale=sapply(designList, function(x) { return(x[[3]]['betaScale'])})
)


# calculate empirical power for each design
# !! Requires several hours to run !!
empiricalPowerAndTimeList = list()
for(i in 1:length(designList)) {
  #for(i in 1:10) {
  x = designList[[i]]
  print(paste(c(i, ": Calculating power for '", x[[1]]@name ,
                "', N=", x[[3]]['perGroupN'], 
                ", SigmaYGscale=", x[[3]]['sigmaYGscale'],
                ", SigmaGscale=", x[[3]]['sigmaGscale']),
              collapse=""))
  startTime <- proc.time()
  power = fastEmpiricalPower(x[[1]], x[[2]])
  ellapsed = proc.time() - startTime
  print(paste(c("Done (", ellapsed[[1]], "s). Power=", power), collapse=""))
  empiricalPowerAndTimeList[[i]] = list(power, ellapsed)
}

## add the timing results and the empirical power values to the data set
empiricalPowerData = cbind(powerData, data.frame(
  empiricalPower=sapply(empiricalPowerAndTimeList, function(x) { return(x[[1]])}),
  time=sapply(empiricalPowerAndTimeList, function(x) {return(x[[2]][1])})
))

## write the calculated and empirical power data to disk   
write.csv(empiricalPowerData, 
          file=paste(c(output.data.dir, "empiricalPower.csv"), collapse="/"))




# add covariate adjusted power with df adjustment
approxPowerData = powerData
approxPowerData$power.covarAdj.mAdjDF = sapply(designList, function(x) { 
  return(glmmPower.covariateAdjusted(x[[1]], x[[2]],
                                     mAdjust=mAdjust.fixedVsRandomDF))})
# add covariate adjusted power with expected projection adjustment
approxPowerData$power.covarAdj.mAdjExpProj = sapply(designList, function(x) { 
  return(glmmPower.covariateAdjusted(x[[1]], x[[2]], 
                                     mAdjust=mAdjust.expectedProjection))})
# add power using method described by
#
# Shieh, G. (2005). Power and sample size calculations for multivariate 
# linear models with random explanatory variables. Psychometrika, 70(2), 347–358. 
#
approxPowerData$power.shieh=sapply(designList, function(x) { 
  return(glmmPower.shieh(x[[1]], x[[2]])) })

# add fixed power, which does not account for covariates. This method described by 
#
# Muller, K. E., Lavange, L. M., Ramey, S. L., & Ramey, C. T. (1992). 
# Power Calculations for General Linear Multivariate Models Including 
# Repeated Measures Applications. Journal of the American Statistical 
# Association, 87(420), 1209–1226.
#
approxPowerData$power.fixedOnly=sapply(designList, function(x) {
  newDesign = getFixedDesign(x[[1]])
  newHypothesis = resizeBetweenContrast(x[[2]], nrow(newDesign@Beta))
  return(glmmPower.fixed(newDesign, newHypothesis)) 
})

# add power which accounts for only the strongest covariate, using the method described by
#
# Glueck, D. H., & Muller, K. E. (2003). Adjusting power for a baseline covariate 
# in linear models. Statistics in Medicine, 22(16), 2535–2551. 
#
approxPowerData$power.topCovar=sapply(designList, function(x) { 
  print(x[[1]]@name)
  newDesign = getSingleCovariateDesign(x[[1]])
  newHypothesis = resizeBetweenContrast(x[[2]], nrow(newDesign@Beta) + 1)
  return(glmmPower.unconditionalSingleCovariate(newDesign, newHypothesis)) 
})

# set the output directory
if (is.null(output.data.dir)) {
  output.data.dir = "."  
}
if (is.null(output.figures.dir)) {
  output.figures.dir = "."
}

# write the calculated power values to disk
write.csv(approxPowerData, file=paste(c(output.data.dir, "approximatePower.csv"), collapse="/"))



#
# Produce summary figures and calculate max deviations for each method
#

# reload the data file
powerData = read.csv(paste(c(output.data.dir, "calculatedAndEmpiricalPower.csv"), 
                                    collapse="/"), header=TRUE, stringsAsFactors=FALSE)

# get rid of sigmaG scales if still in data set
powerData$id = 1:length(powerData$designName)
powerData$diff.covar = powerData$power.covarAdj.mAdjExpProj - powerData$empiricalPower 
powerData$diff.shieh = powerData$power.shieh - powerData$empiricalPower 
powerData$diff.topCovar = powerData$power.topCovar - powerData$empiricalPower
powerData$diff.fixed = powerData$power.fixed - powerData$empiricalPower

# get max absolute deviations
max(abs(powerData$diff.covar))
max(abs(powerData$diff.shieh))
max(abs(powerData$diff.topCovar))
max(abs(powerData$diff.fixed))

# convert to long with factor identifying power method
powerDataLong = reshape(powerData, 
                        varying=c(
                          "diff.covar",
                          "diff.shieh",
                          "diff.topCovar",
                          "diff.fixed"
                        ),
                        timevar="method",
                        times = c("covar", "shieh", "topCovar", "fixed"),
                        idvar="id", direction="long")
# make 'method' into a factor so we can sort the boxplots
powerDataLong$method = factor(powerDataLong$method, 
                              levels=c("covar", "shieh", "topCovar", "fixed"),
                              labels=c("New Method", 
                                       "Shieh", 
                                       "Single Covariate",
                                       "Fixed Only"))
# silly me, didn't put the number of covariates in a separate column
powerDataLong$numCovar = 
  ifelse(grepl("1 covar", powerDataLong$designName),
         1, 
         ifelse(grepl("3 covar", powerDataLong$designName),
                3, 6))

# Plot deviation from empirical across all designs
pdf(file=paste(c(output.figures.dir, "PowerBoxPlot_Overall.pdf"), collapse="/"), family="Times")
par(lab=c(3,3,7))
boxplot(diff ~ method, data=powerDataLong, las=1, ylim=c(-0.6,0.2),
        ylab="Deviation from Empirical Power")
dev.off()

# plot by number of covariates
pdf(file=paste(c(output.figures.dir, "PowerBoxPlot_NumCovar.pdf"), collapse="/"), family="Times")
par(mfrow=c(3,1), oma=c(5,1,1,1), mar=c(1,4,0,0), lab=c(3,3,7))
boxplot(diff ~ method, data=powerDataLong[powerDataLong$numCovar==1,],
        xaxt='n', ylim=c(-0.6, 0.2), las=1, 
        ylab="1 Covariate")
boxplot(diff ~ method, data=powerDataLong[powerDataLong$numCovar==3,],
        xaxt='n', ylim=c(-0.6, 0.2), las=1,
        ylab="3 Covariates")
boxplot(diff ~ method, data=powerDataLong[powerDataLong$numCovar==6,],
        ylab="6 Covariates", las=1, ylim=c(-0.6, 0.2))
dev.off()

# plot by small and large sample size
pdf(file=paste(c(output.figures.dir, "PowerBoxPlot_PerGroupN.pdf"), collapse="/"), family="Times")
par(mfrow=c(2,1), oma=c(5,1,1,1), mar=c(1,4,0,0), lab=c(3,3,7))
boxplot(diff ~ method, data=powerDataLong[powerDataLong$perGroupN==10,],
        xaxt='n', ylim=c(-0.6, 0.2), las=1,
        ylab="Per Group N = 10")
boxplot(diff ~ method, data=powerDataLong[powerDataLong$perGroupN==100,],
        ylim=c(-0.6, 0.2), las=1,
        ylab="Per Group N = 100")
dev.off()

# plot by covariate influence (i.e. SigmaYG-scale)
pdf(file=paste(c(output.figures.dir, "PowerBoxPlot_SigmaYG_Scale.pdf"), collapse="/"), family="Times")
par(mfrow=c(4,1), oma=c(5,1,1,1), mar=c(1,4,0,0), lab=c(3,3,7))
boxplot(diff ~ method, data=powerDataLong[powerDataLong$sigmaYGscale==0.5,],
        xaxt='n', ylab=expression(bold(Sigma)[YG]-scale == 0.5), las=1,
        ylim=c(-0.6, 0.2))
boxplot(diff ~ method, data=powerDataLong[powerDataLong$sigmaYGscale==1,],
        xaxt='n', ylab=expression(bold(Sigma)[YG]-scale == 1), las=1,
        ylim=c(-0.6, 0.2))
boxplot(diff ~ method, data=powerDataLong[powerDataLong$sigmaYGscale==1.5,],
        xaxt='n', ylab=expression(bold(Sigma)[YG]-scale == 1.5), las=1,
        ylim=c(-0.6, 0.2))
boxplot(diff ~ method, data=powerDataLong[powerDataLong$sigmaYGscale==2,],
        ylab=expression(bold(Sigma)[YG]-scale == 2), las=1,
        ylim=c(-0.6, 0.2))
dev.off()

#
# See the output.data.dir and output.figures.dir for results
#


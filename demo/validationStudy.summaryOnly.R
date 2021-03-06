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
# For notation an theoretical details, see
#
# 1. Muller KE, Stewart PW. Linear model theory: univariate, multivariate, and mixed models. 
# Hoboken, New Jersey: John Wiley and Sons; 2006.
#

#
# Extract the fixed portion of a GLMM(F,G) design 
#
getFixedDesign = function(design) {
  if (class(design) == "design.glmmF") {
    return(design)
  } else if (class(design) == "design.glmmFG") {
    return(new("design.glmmF",
               name=paste(c("As design.glmmF: ", design@name), collapse=""),
               description=design@description,                 
               XEssence=design@XEssence,
               perGroupN=design@perGroupN,
               Beta=design@Beta,
               SigmaError=design@SigmaY
    ))    
  } else {
    stop("Input design must either be a design.glmmF or a design.glmmFG")
  }
}

#
# Create a design with a single covariate from a full
# GLMM(F,G) design
#
getSingleCovariateDesign = function(design) {
  if (class(design) == "design.glmmF") {
    return(design)
  } else if (class(design) == "design.glmmFG") {
    if (ncol(design@SigmaG) == 1) {
      return(design)
    } else {
      # select the covariate with the strongest covariance with any outcome
      # get the column of sigmaYG with the maximum covairance
      covarColumn = which.max(apply(design@SigmaYG,2,max))
      # extract the appropriate sigmaG element
      sigmaG = matrix(design@SigmaG[covarColumn, covarColumn])
      # extract the appropriate column from SigmaYG
      sigmaYG = matrix(design@SigmaYG[, covarColumn])
      # return the new design
      return(new("design.glmmFG",
                 name=paste(c("As Single Covariate: ", design@name), collapse=""),
                 description=design@description,
                 XEssence = design@XEssence,
                 perGroupN = design@perGroupN,
                 Beta = design@Beta,
                 SigmaY = design@SigmaY,
                 SigmaG = sigmaG,
                 SigmaYG = sigmaYG           
      ))    
    }
  } else {
    stop("Input design must either be a design.glmmF or a design.glmmFG")
  }
}

#
# Reduce the number of columns in the between participant contrast
#
resizeBetweenContrast <- function(glh, newColumns) {
  return(new("glh",
             alpha = glh@alpha,
             betweenContrast = matrix(glh@betweenContrast[,1:newColumns], 
                                      nrow=nrow(glh@betweenContrast), byrow=FALSE),
             withinContrast = glh@withinContrast,
             thetaNull = glh@thetaNull,
             test = glh@test
  ))
}



#
#
# generateManuscriptResults.R
#
# Generate power results, summary tables, and plots for manuscript
# !!! Be sure to set working directory to this file location prior to running !!!

# load libraries and initialize Java
generateManuscriptResults <- function(runEmpirical=FALSE) {
  
  # generate the designs for the study
  designList = generateDesignsForManuscript()
#   SigmaEList = lapply(designList, function(x) {
#     design = x[[1]]
#     return(design@SigmaY - 
#              design@SigmaYG %*% solve(design@SigmaG) %*% t(design@SigmaYG))
#   })
#   isPosDefSigmaE(designList)
  
  
  # calculate empirical power (or load from existing csv file)
  empiricalPowerData = data.frame(
    designName=sapply(designList, function(x) { return(x[[1]]@name)}),
    perGroupN=sapply(designList, function(x) { return(x[[3]]['perGroupN'])}),
    sigmaYGscale=sapply(designList, function(x) { return(x[[3]]['sigmaYGscale'])}),
    #sigmaGscale=sapply(designList, function(x) { return(x[[3]]['sigmaGscale'])}),
    betaScale=sapply(designList, function(x) { return(x[[3]]['betaScale'])})
  )
  
  # add covariate adjusted power with df adjustment
  empiricalPowerData$power.covarAdj.mAdjDF = sapply(designList, function(x) { 
    return(glmmPower.covariateAdjusted(x[[1]], x[[2]],
                                       mAdjust=mAdjust.fixedVsRandomDF))})
  # add covariate adjusted power with expected projection adjustment
  empiricalPowerData$power.covarAdj.mAdjExpProj = sapply(designList, function(x) { 
    return(glmmPower.covariateAdjusted(x[[1]], x[[2]], 
                                       mAdjust=mAdjust.expectedProjection))})
  # add power using method described by Shieh
  empiricalPowerData$power.shieh=sapply(designList, function(x) { 
    return(glmmPower.shieh(x[[1]], x[[2]])) })
  
  # add fixed power
  empiricalPowerData$power.fixedOnly=sapply(designList, function(x) {
    newDesign = getFixedDesign(x[[1]])
    newHypothesis = resizeBetweenContrast(x[[2]], nrow(newDesign@Beta))
    return(glmmPower.fixed(newDesign, newHypothesis)) 
  })
  
  # add power using the strongest covariate
  empiricalPowerData$power.topCovar=sapply(designList, function(x) { 
    print(x[[1]]@name)
    newDesign = getSingleCovariateDesign(x[[1]])
    newHypothesis = resizeBetweenContrast(x[[2]], nrow(newDesign@Beta) + 1)
    return(glmmPower.unconditionalSingleCovariate(newDesign, newHypothesis)) 
  })
  
  # write the calculated power values to disk
  write.csv(empiricalPowerData, file="../data/calculatedPower.csv")
  
  # calculate empirical power for each design
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

  
  ## add the timing results and the empirical power values to the data
  empiricalPowerData$empiricalPower=sapply(empiricalPowerAndTimeList, function(x) {
    return(x[[1]])})
  empiricalPowerData$time=sapply(empiricalPowerAndTimeList, function(x) {
    return(x[[2]][1])})
                                  
  ## write the empirical power data to disk   
  write.csv(empiricalPowerData, file="../data/calculatedAndEmpiricalPower.csv")
  
  
}


summarizeResults = function() {
  
  # load the data
  empiricalPowerData = read.csv("../data/calculatedAndEmpiricalPower.csv",
                                header=TRUE, stringsAsFactors=FALSE)
  
  # get rid of sigmaG scales if still in data set
  drops <- c("sigmaGscale")
  powerDataClean = empiricalPowerData[,!(names(empiricalPowerData) %in% drops)]
  powerDataClean$id = 1:length(powerDataClean$designName)
  powerDataClean$diff.covar = powerDataClean$power.covarAdj.mAdjExpProj - powerDataClean$empiricalPower 
  powerDataClean$diff.shieh = powerDataClean$power.shieh - powerDataClean$empiricalPower 
  powerDataClean$diff.topCovar = powerDataClean$power.topCovar - powerDataClean$empiricalPower
  powerDataClean$diff.fixed = powerDataClean$power.fixed - powerDataClean$empiricalPower
  
  # get max absolute deviations
  max(abs(powerDataClean$diff.covar))
  max(abs(powerDataClean$diff.shieh))
  max(abs(powerDataClean$diff.topCovar))
  max(abs(powerDataClean$diff.fixed))
  
  # convert to long with factor identifying power method
  powerDataLong = reshape(powerDataClean, 
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
  pdf(file="../inst/figures/PowerBoxPlot_Overall.pdf", family="Times")
  par(lab=c(3,3,7))
  boxplot(diff ~ method, data=powerDataLong, las=1, ylim=c(-0.6,0.2),
          ylab="Deviation from Empirical Power")
  dev.off()
  
  # plot by number of covariates
  pdf(file="../inst/figures/PowerBoxPlot_NumCovar.pdf", family="Times")
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
  pdf(file="../inst/figures/PowerBoxPlot_PerGroupN.pdf", family="Times")
  par(mfrow=c(2,1), oma=c(5,1,1,1), mar=c(1,4,0,0), lab=c(3,3,7))
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$perGroupN==10,],
          xaxt='n', ylim=c(-0.6, 0.2), las=1,
          ylab="Per Group N = 10")
  boxplot(diff ~ method, data=powerDataLong[powerDataLong$perGroupN==100,],
          ylim=c(-0.6, 0.2), las=1,
          ylab="Per Group N = 100")
  dev.off()
  
  # plot by covariate influence (i.e. SigmaYG-scale)
  pdf(file="../inst/figures/PowerBoxPlot_SigmaYG_Scale.pdf", family="Times")
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
  
}


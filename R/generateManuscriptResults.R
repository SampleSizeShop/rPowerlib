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
# Script containing functions which reproduce the results int the manuscript
# "Calculating Power for the General Linear Multivariate Model 
#  With One or More Gaussian Covariates"
# by Sarah M. Kreidler, Keith E. Muller, and Deborah H. Glueck
#
#

#' Extract the fixed-only portion of a GLMM(F,G) design. 
#'
#' Creates a design.glmmF object from a design.glmmFG or design.glmmF
#' object, after removing information regarding covariates.  
#'
#' @param design a design.glmmF or design.glmmFG object
#' @return design.glmmF object
#' 
getFixedDesign <- function(design) {
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

#' Create a design with a single covariate from a full GLMM(F,G) design
#'
#' Creates a new design.glmmFG object with a single covariate.  The
#' covariate with the strongest covariance with any of the outcomes is
#' retained in the design  
#'
#' @param design a design.glmmF or design.glmmFG object
#' @return design.glmmFG object with a single covariate
#' 
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

#' Reduce the number of columns in the between participant contrast
#'
#' Convenience routine to trim the columns of the between participant
#' contrast in a general linear hypothesis  
#'
#' @param glh the general linear hypothesis object
#' @param newColumns the updated number of columns
#' @return design.glmmFG object with a single covariate
#' 
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

#' calculateEmpiricalPowerForDesignList
#' 
#' Calculate empirical power values for a list of designs.
#'   
#' @param glh the general linear hypothesis object
#' @param newColumns the updated number of columns
#' @return writes the empirical power values as a CSV
#' @note This function requires several hours to run
#' 
calculateEmpiricalPowerForDesignList <- function(designList, output.data.dir=".") {
  # calculate empirical power 
  empiricalPowerData = data.frame(
    designName=sapply(designList, function(x) { return(x[[1]]@name)}),
    perGroupN=sapply(designList, function(x) { return(x[[3]]['perGroupN'])}),
    sigmaYGscale=sapply(designList, function(x) { return(x[[3]]['sigmaYGscale'])}),
    #sigmaGscale=sapply(designList, function(x) { return(x[[3]]['sigmaGscale'])}),
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
  empiricalPowerData$empiricalPower=sapply(empiricalPowerAndTimeList, function(x) { return(x[[1]])})
  empiricalPowerData$time=sapply(empiricalPowerAndTimeList, function(x) {return(x[[2]][1])})
  
  ## write the calculated and empirical power data to disk   
  save(empiricalPowerData,
       file=paste(c(output.data.dir, "empiricalPower.RData"), collapse="/"))
  
  return(empiricalPowerData)
}

#' calculateApproximatePowerForDesignList
#' 
#' Calculate approximate power values for a list of design/glh pairs.
#'   
#' @param designList list of pairs of design.glmmFG and glh objects
#' @return design.glmmFG object with a single covariate
#' @note This function requires several hours to run
#' 
calculateApproximatePowerForDesignList <- function(designList, output.data.dir=".") {
  # calculate empirical power 
  approxPowerData = data.frame(
    designName=sapply(designList, function(x) { return(x[[1]]@name)}),
    perGroupN=sapply(designList, function(x) { return(x[[3]]['perGroupN'])}),
    sigmaYGscale=sapply(designList, function(x) { return(x[[3]]['sigmaYGscale'])}),
    #sigmaGscale=sapply(designList, function(x) { return(x[[3]]['sigmaGscale'])}),
    betaScale=sapply(designList, function(x) { return(x[[3]]['betaScale'])})
  )
  
  # add covariate adjusted power with df adjustment
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
  
  ## write the approximate power data to disk     
  save(approxPowerData,
       file=paste(c(output.data.dir, "approximatePower.RData"), collapse="/"))
  
  return(approxPowerData)

}

#' summarizeResults
#' 
#' Generate box plots which summarize the deviations from empirical
#' power for the approximate methods.
#'   
#' @param designList list of pairs of design.glmmFG and glh objects
#' @return design.glmmFG object with a single covariate
#' @note This function requires several hours to run
#'
summarizeResults = function(output.data.dir=".", output.figures.dir=".") {
  
  # load the data - creates a variable 'approximateAndEmpiricalPowerData'
  load(paste(c(output.data.dir, "approximateAndEmpiricalPower.RData"), collapse="/"))
  
  # calculate the deviations
  approximateAndEmpiricalPower$id = 1:length(approximateAndEmpiricalPower$designName)
  approximateAndEmpiricalPower$diff.covar = approximateAndEmpiricalPower$power.covarAdj.mAdjExpProj - 
    approximateAndEmpiricalPower$empiricalPower 
  approximateAndEmpiricalPower$diff.shieh = approximateAndEmpiricalPower$power.shieh - 
    approximateAndEmpiricalPower$empiricalPower 
  approximateAndEmpiricalPower$diff.topCovar = approximateAndEmpiricalPower$power.topCovar - 
    approximateAndEmpiricalPower$empiricalPower
  approximateAndEmpiricalPower$diff.fixed = approximateAndEmpiricalPower$power.fixed - 
    approximateAndEmpiricalPower$empiricalPower
  
  # get max absolute deviations
  maxDeviationData = data.frame(method=c("covar", "shieh", "topCovar", "fixed"),
                                maxDeviation=c(max(abs(approximateAndEmpiricalPower$diff.covar)),
                                               max(abs(approximateAndEmpiricalPower$diff.shieh)),
                                               max(abs(approximateAndEmpiricalPower$diff.topCovar)),
                                               max(abs(approximateAndEmpiricalPower$diff.fixed))))
  # save the max deviation info to a file
  save(maxDeviationData,
       file=paste(c(output.data.dir, "maxDeviationsByMethod.RData"), collapse="/"))
  
  # convert to long with factor identifying power method
  powerDataLong = reshape(approximateAndEmpiricalPower, 
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
  
}

#' runSimulationStudy
#' 
#' This function reproduces the simulation study results for the manuscript:\cr
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. Calculating Power for the General 
#' Linear Multivariate Model With One or More Gaussian Covariates, In review.
#'
#' @param study.seed the random number seed (defaults to 7634)
#' @param study.data.dir the directory into which data files are written (defaults to
#' current working directory)
#' @param study.figures.dir the directory into which pdf figures are written (defaults
#' to the current working directory)
#' @param study.runEmpirical if true, empirical power values will be recalculated. If false,
#' existing empirical power values will be loaded from the R package.
#' 
#' @note 
#' The empirical power calculations may take several hours to run
#'
runSimulationStudy <- function(study.seed=7634, study.data.dir=".", study.figures.dir=".",
                               study.runEmpirical=TRUE) {
  # set the random seed
  set.seed(study.seed)
  
  # generate the designs for the validation study
  designList = generateDesignsForManuscript()
  
  # calculate empirical power
  if (study.runEmpirical) {
    # calculate empirical power for each design
    # !! Requires several hours to run !!
    empiricalPowerData = calculateEmpiricalPowerForDesignList(designList, study.data.dir)
    
  } else {
    # load the existing empirical data 
    data("empiricalPower", package="rPowerlib")
  }
  
  # calculate approximate power - runs in about 10-30 minutes
  approxPowerData = calculateApproximatePowerForDesignList(designList, study.data.dir)

  # combine the data into a single data set and write to disk
  approxAndEmpiricalPowerData = data.frame(approxPowerData, 
                                           empiricalPower=empiricalPowerData$empiricalPower,
                                           empiricalTime=empiricalPowerData$time)
  # save to disk
  save(approxAndEmpiricalPowerData,
       file=paste(c(study.data.dir, "approximateAndEmpiricalPower.RData"), collapse="/"))
  
  # Produce summary figures and calculate max deviations for each method
  summarizeResults(study.data.dir, study.figures.dir)

  
}

#' runElashoffExample
#' 
#' This function reproduces the applied example results for the manuscript:\cr
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. Calculating Power for the General 
#' Linear Multivariate Model With One or More Gaussian Covariates, In review.
#'
#' @param study.seed the random number seed (defaults to 7634)
#' @param study.data.dir the directory into which data files are written (defaults to
#' current working directory)
#' @param study.figures.dir the directory into which pdf figures are written (defaults
#' to the current working directory)
#' @param study.runEmpirical if true, empirical power values will be recalculated. If false,
#' existing empirical power values will be loaded from the R package.
#' 
#' @note 
#' The example is based on the study by 
#' 
#' Elashoff, D., Zhou, H., Reiss, J., Wang, J., Xiao, H., Henson, B., … Wong, D. T. W. (2012). 
#' Prevalidation of salivary biomarkers for oral cancer detection. Cancer Epidemiology, 
#' Biomarkers & Prevention, 21(4), 664–672.
#'  
runElashoffExample = function() {
  demo(appliedExample, package="rPowerlib")
}

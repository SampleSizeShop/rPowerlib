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

# Define the X essence matrices
XEssence.trt2 = diag(2)
XEssence.trt4 = diag(4)
perGroupNList = c(10, 100)

# define covariances for covariates (designs with 1, 3, or 6 covariates)
rho.G = 0.01
SigmaG.cov1 = matrix(c(1))
SigmaG.cov3 = rho.G * as.matrix(rep(1,3)) %*% t(as.matrix(rep(1,3))) + diag(3)*(1-rho.G)
SigmaG.cov6 = rho.G * as.matrix(rep(1,6)) %*% t(as.matrix(rep(1,6))) + diag(6)*(1-rho.G)
SigmaG.scaleList = c(1,2)

# define sigma Y for univariate and multivariate designs
rho.Y = 0.04
SigmaY.univariate = matrix(c(1.3))
SigmaY.LEAR = matrix(c(1.0,0.4,0.3,0.2,0.16,
                               0.4,1.0,0.4,0.3,0.2,
                               0.3,0.4,1.0,0.4,0.3,
                               0.2,0.3,0.4,1.0,0.4,
                               0.16,0.2,0.3,0.4,1.0), nrow=5)
SigmaY.CS = rho.Y * as.matrix(rep(1,5)) %*% t(as.matrix(rep(1,5))) + diag(5)*(1-rho.Y)
SigmaY.CSxCS = (
  (0.1 * as.matrix(rep(1,3)) %*% t(as.matrix(rep(1,3))) + diag(3)*(1-0.1)) %x%
    (0.2 * as.matrix(rep(1,2)) %*% t(as.matrix(rep(1,2))) + diag(2)*(1-0.2))
  )

# define beta fixed for unviariate and multivariate designs
betaFixed.trt2Uni = matrix(c(1,0), nrow=2) 
betaFixed.trt2Multi5 = matrix(c(1,rep(0,9)), nrow=2) 
betaFixed.trt2Multi6 = matrix(c(1,rep(0,11)), nrow=2) 
betaFixed.trt4Uni = matrix(c(1,0,0,0), nrow=4) 
betaFixed.trt4Multi5 = matrix(c(1,rep(0,19)), nrow=4) 
betaFixed.trt4Multi6 = matrix(c(1,rep(0,23)), nrow=4) 

# define sigmaYG (p x qg), p=#responses, qg=#covariates for univariate and multivariate designs
rho.YG = 0.6
# one covariate
SigmaYG.cov1Uni = matrix(c(rho.YG))
SigmaYG.cov1Lear = matrix(sapply(seq(1,2,length.out=5), function(x) { rho.YG^(x)}), nrow=5)
SigmaYG.cov1CS = matrix(rep(rho.YG,5), nrow=5)
SigmaYG.cov1CSxCS = (matrix(c(0.3, 0.3, 0.3), nrow=3) %x% matrix(c(0.6, 0.5), nrow=2))
# 3 covariates
SigmaYG.cov3Uni = SigmaYG.cov1Uni %x% matrix(seq(1,0.5, length.out=3), nrow=1)
SigmaYG.cov3Lear = matrix(seq(1,0.5, length.out=3), nrow=1) %x% SigmaYG.cov1Lear
SigmaYG.cov3CS = matrix(seq(1,0.5, length.out=3), nrow=1) %x% SigmaYG.cov1CS
SigmaYG.cov3CSxCS = matrix(seq(1,0.5, length.out=3), nrow=1) %x% SigmaYG.cov1CSxCS
# 6 covariates
SigmaYG.cov6Uni = SigmaYG.cov1Uni %x% matrix(seq(1,0.5,length.out=6), nrow=1)
SigmaYG.cov6Lear = matrix(seq(1,0.5,length.out=6), nrow=1) %x% SigmaYG.cov1Lear
SigmaYG.cov6CS = matrix(seq(1,0.5,length.out=6), nrow=1) %x% SigmaYG.cov1CS
SigmaYG.cov6CSxCS = matrix(seq(1,0.5,length.out=6), nrow=1) %x% SigmaYG.cov1CSxCS
# scale factors for SigmaYG
SigmaYG.scaleList = c(0.5, 1, 2)

# build the hypothesis objects for univariate and multivariate designs
glh.trt2uni = new("glh", 
                  betweenContrast=matrix(c(1,-1,0),nrow=1), 
                  withinContrast=matrix(c(1)), 
                  thetaNull=matrix(c(0), nrow=1))
glh.trt2multi5 = new("glh", 
                  betweenContrast=matrix(c(1,-1),nrow=1), 
                  withinContrast=matrix(poly(1:5,degree=4),nrow=5), 
                  thetaNull=matrix(rep(0,4), ncol=4))
glh.trt2doubly6 = new("glh", 
                     betweenContrast=matrix(c(1,-1),nrow=1), 
                     withinContrast=t(matrix(c(1,-1,0,1,0,-1), nrow=2, byrow=TRUE) %x% matrix(c(1,-1),nrow=1)), 
                     thetaNull=matrix(rep(0,2), ncol=2))
# hypotheses for 4 group
glh.trt4uni = new("glh", 
                  betweenContrast=cbind(rep(1,3), -1*diag(3)), 
                  withinContrast=matrix(c(1)), 
                  thetaNull=matrix(c(rep(0,3)), nrow=3))
glh.trt4multi5 = new("glh", 
                     betweenContrast=cbind(rep(1,3), -1*diag(3), matrix(rep(0,9),nrow=3)), 
                     withinContrast=matrix(poly(1:5,degree=4),nrow=5), 
                     thetaNull=matrix(rep(0,12), ncol=4))
glh.trt4doubly6 = new("glh", 
                      betweenContrast=cbind(rep(1,3), -1*diag(3)), 
                      withinContrast=t(matrix(c(1,-1,0,1,0,-1), nrow=2, byrow=TRUE) %x% matrix(c(1,-1),nrow=1)), 
                      thetaNull=matrix(rep(0,6), nrow=3))

#
# create a data frame with combinations of per group N,
# sigmaY scale, sigmaYG scale, and beta scale
#
paramList = list(perGroupN=c(10,100), sigmaYGscale=c(0.5, 1, 2), sigmaGscale=c(1, 2))
paramComboList = data.frame(expand.grid(paramList))


#
# Update the design and/or glh with the specified parameters
#
setDesignParameters <- function(params, design, glh) {
  
  # set the per group N
  design@perGroupN = params['perGroupN']
  # scale the sigma YG matrix
  design@SigmaYG = params['sigmaYGscale'] * design@SigmaYG
  # scale the sigma G matrix
  design@SigmaG = params['sigmaGscale'] * design@SigmaG

  return(list(design, glh, params))
}

#
# Build the final design and hypothesis list
#
manuscriptDesignList = c(
  ### Univariate designs
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 1 covariate, univariate", 
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Uni, 
            SigmaY=SigmaY.univariate, 
            SigmaG=SigmaG.cov1, 
            SigmaYG=SigmaYG.cov1Uni), glh.trt2uni)),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 3 covariates, univariate",
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Uni, 
            SigmaY=SigmaY.univariate, 
            SigmaG=SigmaG.cov3, 
            SigmaYG=SigmaYG.cov3Uni), glh.trt2uni),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 6 covariates, univariate", 
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Uni, 
            SigmaY=SigmaY.univariate, 
            SigmaG=SigmaG.cov6, 
            SigmaYG=SigmaYG.cov6Uni), glh.trt2uni),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 1 covariate, univariate", 
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Uni, 
            SigmaY=SigmaY.univariate, 
            SigmaG=SigmaG.cov1, 
            SigmaYG=SigmaYG.cov1Uni), glh.trt2uni),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 3 covariates, univariate",
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Uni, 
            SigmaY=SigmaY.univariate, 
            SigmaG=SigmaG.cov3, 
            SigmaYG=SigmaYG.cov3Uni), glh.trt2uni),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 6 covariates, univariate", 
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Uni, 
            SigmaY=SigmaY.univariate, 
            SigmaG=SigmaG.cov6, 
            SigmaYG=SigmaYG.cov6Uni), glh.trt2uni),
  
  ### Multivariate designs with compound symmetric SigmaY
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 1 covariate, multivariate CS", 
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Multi5, 
            SigmaY=SigmaY.CS, 
            SigmaG=SigmaG.cov1, 
            SigmaYG=SigmaYG.cov1CS), glh.trt2multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 3 covariates, multivariate CS",
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Multi5, 
            SigmaY=SigmaY.CS, 
            SigmaG=SigmaG.cov3, 
            SigmaYG=SigmaYG.cov3CS), glh.trt2multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 6 covariates, multivariate CS", 
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Multi5, 
            SigmaY=SigmaY.CS, 
            SigmaG=SigmaG.cov6, 
            SigmaYG=SigmaYG.cov6CS), glh.trt2multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 1 covariate, multivariate CS", 
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Multi5, 
            SigmaY=SigmaY.CS, 
            SigmaG=SigmaG.cov1, 
            SigmaYG=SigmaYG.cov1CS), glh.trt4multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 3 covariates, multivariate CS",
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Multi5, 
            SigmaY=SigmaY.CS, 
            SigmaG=SigmaG.cov3, 
            SigmaYG=SigmaYG.cov3CS), glh.trt4multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 6 covariates, multivariate CS", 
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Multi5, 
            SigmaY=SigmaY.CS, 
            SigmaG=SigmaG.cov6, 
            SigmaYG=SigmaYG.cov6CS), glh.trt4multi5),
  
  
  ### Multivariate designs with LEAR SigmaY
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 1 covariate, multivariate LEAR", 
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Multi5, 
            SigmaY=SigmaY.LEAR, 
            SigmaG=SigmaG.cov1, 
            SigmaYG=SigmaYG.cov1Lear), glh.trt2multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 3 covariates, multivariate LEAR",
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Multi5, 
            SigmaY=SigmaY.LEAR, 
            SigmaG=SigmaG.cov3, 
            SigmaYG=SigmaYG.cov3Lear), glh.trt2multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 6 covariates, multivariate LEAR", 
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Multi5, 
            SigmaY=SigmaY.LEAR, 
            SigmaG=SigmaG.cov6, 
            SigmaYG=SigmaYG.cov6Lear), glh.trt2multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 1 covariate, multivariate LEAR", 
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Multi5, 
            SigmaY=SigmaY.LEAR, 
            SigmaG=SigmaG.cov1, 
            SigmaYG=SigmaYG.cov1Lear), glh.trt4multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 3 covariates, multivariate LEAR",
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Multi5, 
            SigmaY=SigmaY.LEAR, 
            SigmaG=SigmaG.cov3, 
            SigmaYG=SigmaYG.cov3Lear), glh.trt4multi5),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 6 covariates, multivariate LEAR", 
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Multi5, 
            SigmaY=SigmaY.LEAR, 
            SigmaG=SigmaG.cov6, 
            SigmaYG=SigmaYG.cov6Lear), glh.trt4multi5),
  
  ### Doubly Multivariate designs with CS@CS SigmaY
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 1 covariate, multivariate CS@CS", 
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Multi6, 
            SigmaY=SigmaY.CSxCS, 
            SigmaG=SigmaG.cov1, 
            SigmaYG=SigmaYG.cov1CSxCS), glh.trt2doubly6),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 3 covariates, multivariate CS@CS",
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Multi6, 
            SigmaY=SigmaY.CSxCS, 
            SigmaG=SigmaG.cov3, 
            SigmaYG=SigmaYG.cov3CSxCS), glh.trt2doubly6),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Two fixed predictors, 6 covariates, multivariate CS@CS", 
            XEssence=XEssence.trt2,
            Beta=betaFixed.trt2Multi6, 
            SigmaY=SigmaY.CSxCS, 
            SigmaG=SigmaG.cov6, 
            SigmaYG=SigmaYG.cov6CSxCS), glh.trt2doubly6),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 1 covariate, multivariate CS@CS", 
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Multi6, 
            SigmaY=SigmaY.CSxCS, 
            SigmaG=SigmaG.cov1, 
            SigmaYG=SigmaYG.cov1CSxCS), glh.trt4doubly6),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 3 covariates, multivariate CS@CS",
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Multi6, 
            SigmaY=SigmaY.CSxCS, 
            SigmaG=SigmaG.cov3, 
            SigmaYG=SigmaYG.cov3CSxCS), glh.trt4doubly6),
  apply(paramComboList, 1, setDesignParameters, 
        new("design.glmmFG", name="Four fixed predictors, 6 covariates, multivariate CS@CS", 
            XEssence=XEssence.trt4,
            Beta=betaFixed.trt4Multi6, 
            SigmaY=SigmaY.CSxCS, 
            SigmaG=SigmaG.cov6, 
            SigmaYG=SigmaYG.cov6CSxCS), glh.trt4doubly6)

  )


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
                           )




  
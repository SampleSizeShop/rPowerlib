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
# Build the timing table in the manuscript
#
# "Calculating Power for the General Linear Multivariate Model 
#  With One or More Gaussian Covariates"
# by Sarah M. Kreidler, Keith E. Muller, and Deborah H. Glueck
#
#
#
buildTimingTable <- function() {
  # load the data - creates a variable 'approximateAndEmpiricalPowerData'
  load("../data/approximateAndEmpiricalPower.RData")
  
  data = approximateAndEmpiricalPowerData
  # silly me, didn't put the number of covariates in a separate column
  data$numCovar = 
    ifelse(grepl("1 covar", approximateAndEmpiricalPowerData$designName),
           1, 
           ifelse(grepl("3 covar", approximateAndEmpiricalPowerData$designName),
                  3, 6))
  
  data$numRM = 
    ifelse(grepl("CS@CS", approximateAndEmpiricalPowerData$designName),
           6, 
           ifelse(grepl("univariate", approximateAndEmpiricalPowerData$designName), 1, 5))
  
  par(mfrow=c(3,1))
  for(i in unique(data$numRM)) {
    plot(data[data$perGroupN==10 & data$numRM ==i,]$numCovar,
         data[data$perGroupN==10 & data$numRM ==i,]$empiricalTime, col="red", pch=1, 
         xlab="Number of Covariates", ylab="Time (seconds)", 
         ylim=c(0,650))
    points(data[data$perGroupN==100 & data$numRM ==i,]$numCovar,
           data[data$perGroupN==100 & data$numRM ==i,]$empiricalTime, col="blue", pch=1)
    
  }
  
  tmpData = data[data$empiricalTime < 300 & data$empiricalTime > 200,]
  
  
  #
  # Create combinations of study design params for table
  #
  paramList = list(numRM=c(1,5,6), perGroupN=c(10,100), numCovar=c(1,3,6))
  paramComboList = data.frame(expand.grid(paramList))
  paramComboList = paramComboList[,c(3,2,1)]
  
  tmp = aggregate(data$empiricalTime, by=list(data$numRM, data$perGroupN, data$numCovar), FUN=mean)
  paramComboList$empiricalTime = round(tmp$x, digits=1)
  
  tmp = aggregate(data$time.covarAdj.mAdjExpProj, by=list(data$numRM, data$perGroupN, data$numCovar), FUN=mean)
  paramComboList$calculationTime = round((tmp$x / 10^-4), digits=0)

  tmp = aggregate(data$time.covarAdj.mAdjExpProj, by=list(data$numRM, data$perGroupN, data$numCovar), FUN=length)
  paramComboList$numStudies = tmp$x
  
  colnames(paramComboList) = c("\\multicolumn{1}{r}{ \n \\begin{tabular}[t]{c} \n Number of \\\\ covariates \n \\end{tabular} \n}",
                               "\\multicolumn{1}{r}{ \n \\begin{tabular}[t]{c} \n Per group \\\\ sample size \\end{tabular} \n}", 
                               "\\multicolumn{1}{r}{ \n \\begin{tabular}[t]{c} \n Number of \\\\ repeated \\\\ measures \\end{tabular} \n}", 
                               "\\multicolumn{1}{r}{ \n \\begin{tabular}[t]{c} \n Mean \\\\ simulation \\\\ time \\\\ (sec.) \\end{tabular} \n}",
                               "\\multicolumn{1}{r}{ \n \\begin{tabular}[t]{c} \n Mean \\\\ calculation \\\\ time \\\\ (sec. $\\times 10^{-4}$) \\end{tabular} \n}",
                               "\\multicolumn{1}{r}{ \n \\begin{tabular}[t]{c} \n Number of \\\\ designs \\\\ considered \\end{tabular} \n}")
  
  timingTable = xtable(paramComboList, display=c("d","d","d","d","d","d", "d"),
                       caption="Comparison of simulation time to calculation time")
  align( timingTable ) <- c(rep('r',7) )
  print(timingTable, include.rownames=FALSE, sanitize.text.function = function(x){x})
  
}

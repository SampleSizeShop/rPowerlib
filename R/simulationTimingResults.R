#
#
#

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

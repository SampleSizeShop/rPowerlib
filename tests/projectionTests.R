design.test = new("design.glmmFG", name="Four fixed predictors, 3 covariates, multivariate LEAR",
    XEssence=XEssence.trt4,
    Beta=betaFixed.trt4Multi5, 
    SigmaY=SigmaY.LEAR, 
    SigmaG=SigmaG.cov3, 
    SigmaYG=SigmaYG.cov3Lear)

glh.test = glh.trt4multi5


empiricalPower = fastEmpiricalPower(design.test, glh.test)
covarAdjPowerSarah = glmmPower.covariateAdjusted(design.test, glh.test)
covarAdjPowerTheory = glmmPower.covariateAdjustedTheory(design.test, glh.test)

cat("Empirical: ", empiricalPower, " Sarah: ", covarAdjPowerSarah, " Theory: ", covarAdjPowerTheory)

empiricalPower2 = fastEmpiricalPower(design.cov1, glh.cov1)
covarAdjPowerSarah2 = glmmPower.covariateAdjusted(design.cov1, glh.cov1)
covarAdjPowerTheory2 = glmmPower.covariateAdjustedTheory(design.cov1, glh.cov1)

cat("Empirical: ", empiricalPower2, " Sarah: ", covarAdjPowerSarah2, " Theory: ", covarAdjPowerTheory2)

design.test@perGroupN = 100
empiricalPower3 = fastEmpiricalPower(design.test, glh.test)
covarAdjPowerSarah3 = glmmPower.covariateAdjusted(design.test, glh.test)
covarAdjPowerTheory3 = glmmPower.covariateAdjustedTheory(design.test, glh.test)

cat("Empirical: ", empiricalPower3, " Sarah: ", covarAdjPowerSarah3, " Theory: ", covarAdjPowerTheory3)

SigmaG.cov3 = matrix(c(51,3,1.5,3,2,3,1.6,3,13), nrow=3, byrow=TRUE)
G = mvrnorm(n=5, matrix(rep(0,3), nrow=3), SigmaG.cov3)

H = G %*% solve(t(G) %*% G) %*% t(G)
eigen(H)$values

nList = c(10,15,20,25,30,100)
nList = c(50)
iter = 10
Gsum3List = lapply(1:length(nList), function(i) {

  GList = lapply(1:iter, function(x, n) {
    G = mvrnorm(n=n, matrix(rep(0,3), nrow=3), SigmaG.cov3)
    return (diag(n) - G %*% solve(t(G) %*% G) %*% t(G))
  }, nList[i])
  
  return((1/iter) * Reduce("+", GList))
})

Gsum3List2 = lapply(1:length(nList), function(i) {
  
  GList = lapply(1:10000, function(x, n) {
    G = mvrnorm(n=n, matrix(rep(0,3), nrow=3), SigmaG.cov3)
    return (G %*% solve(t(G) %*% G) %*% t(G))
  }, nList[i])
  
  return((1/10000) * Reduce("+", GList))
})

Gsum6List = sapply(1:length(nList), function(i) {
  
  GList = lapply(1:10000, function(x, n) {
    G = mvrnorm(n=n, matrix(rep(0,6), nrow=6), 2*SigmaG.cov6)
    return (diag(n) - G %*% solve(t(G) %*% G) %*% t(G))
  }, nList[i])
  
  return(sum(diag((1/10000) * Reduce("+", GList)))/nList[i])
})


plot(num, (num-1)/num, col="red", "l")
points(nList, Gsum3List)
lines(num, (num-3)/num)
points(nList, Gsum6List, col="green")
lines(num, (num-6)/num, col="green")






#
#

#### Example GLMM fixed ####
# test design
designF = new("design.glmmF", 
              Beta=matrix(c(1,0,0,0,0,0), nrow=2), 
              SigmaError=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1), nrow=3))


#
# Example multivariate GLMM_FG
#
design.cov3 = new("design.glmmFG", 
                  Beta=matrix(c(1,0,0,0,0,0), nrow=2), 
                  SigmaY=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1), nrow=3), 
                  SigmaG=matrix(c(1,0.01,0.01,0.01,1,0.01,0.01,0.01,1), nrow=3), 
                  SigmaYG=matrix(rep(c(0.2,0.4,0.5), 3), nrow=3, byrow=TRUE))
design.cov1 = new("design.glmmFG", 
                  Beta=matrix(c(1,0,0,0,0,0), nrow=2), 
                  SigmaY=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1), nrow=3), 
                  SigmaG=matrix(c(1.2), nrow=1), 
                  SigmaYG=matrix(rep(c(0.2,0.4,0.5), 1), nrow=3, byrow=TRUE))
# startTime <- proc.time()
#  fastEmpiricalPower(design, glh, 1000,1000)
#  ellapsed = proc.time() - startTime
# startTime <- proc.time()
# empiricalPower(design, glh, replicates=1000, realizations=1000)
# ellapsed = proc.time() - startTime
# 
# startTime <- proc.time()
# simulateData(design, outputDir="../data/", blockSize=1000)
# ellapsed = proc.time() - startTime
# 
# Rprof("prof.out")
# prof.out = profr(empiricalPower(design, glh, replicates=1000, realizations=10), quiet=FALSE)
# Rprof(NULL)
# ellapsed = proc.time() - startTime





#### Example hypotheses ####
glh.cov3 = new("glh", betweenContrast=matrix(c(1,-1,0,0,0),nrow=1), 
               withinContrast=matrix(c(1,1,-1,0,0,-1), nrow=3,byrow=TRUE), 
               thetaNull=matrix(c(0,0), nrow=1))

glh.cov1 = new("glh", betweenContrast=matrix(c(1,-1,0),nrow=1), 
               withinContrast=matrix(c(1,1,-1,0,0,-1), nrow=3,byrow=TRUE), 
               thetaNull=matrix(c(0,0), nrow=1))

glhF = new("glh", betweenContrast=matrix(c(1,-1),nrow=1), 
           withinContrast=matrix(c(1,1,-1,0,0,-1), nrow=3,byrow=TRUE), 
           thetaNull=matrix(c(0,0), nrow=1))
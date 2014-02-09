#
# Example study designs for the covariate paper
#
#
#

setClass("StudyDesign", 
  representation(
    alpha="numeric",
    XFixed="matrix",
    KStar="matrix",
    sigmaG="matrix",
    sigmaY="matrix",
    sigmaYG="matrix",
    betaFixed="matrix",
    CMatrix="matrix",
    UMatrix="matrix",
    thetaNull="matrix"
  )
)


powerGLMMFG.example1 <- function() {
  design <- new("StudyDesign",
	alpha=0.01,
	XFixed=(diag(2) %x% matrix(rep(1,10))),
	KStar = (matrix(c(0))),
	sigmaG = (diag(c(1,4))),
	sigmaY = (matrix(ncol=3,byrow=1,c(2,0.4,0.2,0.4,2,0.4,0.2,0.4,2))),
	sigmaYG = (matrix(ncol=2,c(0.2,0.2,0.2,0.1,0.1,0.1)) %*% diag(c(2,8))),
	betaFixed = (matrix(ncol=3,byrow=1,c(0.2,0,0,0,0,0))),
	CMatrix = (matrix(ncol=2,byrow=1,c(1,-1))),
	UMatrix = (t(cbind(matrix(rep(1,2)),-1*diag(2)))),
	thetaNull = (matrix(ncol=2,byrow=1,c(0,0))))

  
  return(design)
}
/*
* Simulates unconditional power for GLMM(F,G) designs
* 
* 
*
*/

proc iml;

/* Simulation parameters */
numSamples = {1000};
numSimulations = {1000};


/* Study design */
typeIError = {0.05};

perGroupN = {20};
/* Es(X) for the fixed predictors */
essenceXFixed = {1 0, 0 1};
/* covariance matrix for Gaussian covariates */
sigmaG = {3 0 0,
		  0 4 0
		  0 0 1};
/* covariance of Gaussian covariates with Y */
sigmaYG = {};
/* covariance of Y without controlling for covariates */
sigmaY = {};
/* choices for parameters */
beta = {};
/* between sampling unit contrast matrix */
C = {};
/* within sampling unit */
U = {};

 
do 1 to numSamples;
  sigmaE = calculateSigmaError(sigmaG, sigmaYG, sigmaY);
  X = getDesignMatrix(essenceXFixed, sigmaG);

  do 1 to numSimulations;
    Y = getResponses(sigmaE);
	pvalue = fitGeneralLinearHypothesis(X, Y, C, U);
  end;
end;

quit;

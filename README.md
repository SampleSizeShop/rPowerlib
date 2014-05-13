The rPowerlib Package
=========================

The rPowerlib package for R (>3.0.0) calculates power for the general
linear multivariate model with or without Gaussian covariates.  

The package provides companion code for the manuscript:

Kreidler, S. M., Muller, K. E., & Glueck, D. H. 
Calculating Power for the General Linear Multivariate Model 
With One or More Gaussian Covariates, In review.

### Instructions for replicating the manuscript results 

The results in the above manuscript were produced using R version 3.0.0. To reproduce the results,
perform the following steps:

* Install R version 3.0.x or higher by following the instructions at http://www.r-project.org
* From the R environment, install and load the "devtools" package
```R
> install.packages("devtools")
> library(devtools)
```
* Install the "rPowerlib" package directly from Github.com
```R
> install_github(repo="rPowerlib", user="SampleSizeShop", ref="develop")
```
* Load the library
```R
> library(rPowerlib)
```
* Run the simulation study (may take several hours to run with empirical calculations). You may specify output directories for data files (study.data.dir) and figures (study.figures.dir). If omitted, both default to the current working directories. To save time, you may optionally skip the empirical calculations by setting study.runEmpirical=FALSE.
```R
> runSimulationStudy(study.data.dir="myDataDir", study.figures.dir="myFiguresDir", study.runEmpirical=TRUE)
```
* Run the applied example
```R
> runElashoffExample()
```

### Background ###

This package is based on the SAS/IML product POWERLIB (http://www.jstatsoft.org/v30/i05)
developed by Keith E. Muller and colleagues.
At present, the rPowerlib package DOES NOT include the full functionality of POWERLIB.
For features such as power confidence intervals and more statistical tests,
please try POWERLIB (http://github.com/samplesizeshop/powerlib) or 
GLIMMPSE (https://glimmpse.samplesizeshop.org). Please
visit http://SampleSizeShop.org for more information about calculating power and
sample size.
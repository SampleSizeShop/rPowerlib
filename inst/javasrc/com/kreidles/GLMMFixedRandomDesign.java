/*
 * Simulation study for the manuscript "Selecting Sample Size for
 * the General Linear Multivariate Model with One or More
 * Gaussian Covariates"
 * 
 * Copyright (C) 2014 Sarah Kreidler.  
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
package com.kreidles;

import java.io.IOException;
import java.util.List;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.StatUtils;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.PowerException;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public class GLMMFixedRandomDesign {
    /* matrix inputs to the power calculation.  */
    // the design essence matrix 
    /* For details please see Muller & Fetterman (2002) "Regression and ANOVA" */
    private RealMatrix XFixed = null;
    // essence matrix
    private RealMatrix XEssence = null;
    // per group sample size
    private int perGroupN = -1;
    // caching of X'X inverse and rank since these are order(n^3) operations
    private RealMatrix XtXInverse = null;

    // beta matrix of regression coefficients associated with fixed predictors
    RealMatrix betaFixed = null;
    // beta matrix of regression coefficients associated with random predictors
    RealMatrix betaRandom = null;
    // combined beta matrix
    RealMatrix betaFull = null;

    // residual error matrix - calculated from the remaining elements
    private RealMatrix sigmaE = null;
    // covariance for outcomes (assuming no covariates)
    private RealMatrix sigmaY = null;
    // covariance of Gaussian predictors
    private RealMatrix sigmaG = null;
    // covariance of outcomes and random predictors
    private RealMatrix sigmaYG = null;

    /**
     * Create a design from a JSON representation of the study.
     * Allows the design to be passed down from R.
     * 
     * @param designAsJSON
     */
    public GLMMFixedRandomDesign(String designAsJSON) 
            throws IOException, ParseException, IllegalArgumentException {
        JSONParser jsonParser = new JSONParser();
        JSONObject jsonObject = (JSONObject) jsonParser.parse(designAsJSON);

        // get the fixed portion of the design matrix
        JSONObject XFixedJSONObject = (JSONObject) jsonObject.get("XFixed");
        JSONObject XEssenceJSONObject = (JSONObject) jsonObject.get("XEssence");
        if (XFixedJSONObject != null) {
            XFixed = JSONUtils.realMatrixFromJsonObject(XFixedJSONObject);
        } else if (XEssenceJSONObject != null) {
            XEssence = JSONUtils.realMatrixFromJsonObject(XEssenceJSONObject);
            // per group sample size
            String perGroupNString = jsonObject.get("perGroupN").toString();
            if (perGroupNString != null) {
                perGroupN = Integer.parseInt(perGroupNString);
                if (perGroupN < 2) {
                    throw new IllegalArgumentException("Invalid per group sample size: must be greater than 1");
                }
            }
        } else {
            throw new IllegalArgumentException("No X matrix specified"); 
        }

        // get the fixed portion of the beta matrix
        betaFixed = JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get("betaFixed"));
        if (betaFixed == null) {
            throw new IllegalArgumentException("No beta matrix specified");
        }
        // get the covariance matrices
        sigmaG = JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get("sigmaG"));
        if (sigmaG == null) {
            throw new IllegalArgumentException("No sigmaG matrix specified");
        }
        sigmaY = JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get("sigmaY"));
        if (sigmaY == null) {
            throw new IllegalArgumentException("No sigmaY matrix specified");
        }
        sigmaYG = JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get("sigmaYG"));
        if (sigmaYG == null) {
            throw new IllegalArgumentException("No sigmaYG matrix specified");
        }

        // Calculate the covariance of errors   
        RealMatrix invSigmaG = new LUDecomposition(this.sigmaG).getSolver().getInverse();
        this.sigmaE = this.sigmaY.subtract(this.sigmaYG.multiply(invSigmaG).multiply(this.sigmaYG.transpose()));
       
        // calculate the random portion of the beta matrix
        this.betaRandom = invSigmaG.multiply(this.sigmaYG.transpose());     
        // form the complete Beta matrix
        this.betaFull = getVerticalAppend(this.betaFixed, this.betaRandom);

    }


    /**
     * Vertically append two matrices
     * @param matrix
     * @param column
     * @return the combined matrix
     * @throws IllegalArgumentException
     */
    private RealMatrix getVerticalAppend(RealMatrix m1, RealMatrix m2)
            throws IllegalArgumentException {
        if (m1 == null || m2 == null)
            throw new IllegalArgumentException("Missing required argument");
        if (m1.getColumnDimension() != m2.getColumnDimension())
            throw new IllegalArgumentException("Column dimensions must be equal");

        RealMatrix newMatrix = 
                new Array2DRowRealMatrix(m1.getRowDimension()+m2.getRowDimension(),
                        m1.getColumnDimension());
        newMatrix.setSubMatrix(m1.getData(), 0, 0);
        newMatrix.setSubMatrix(m2.getData(), m1.getRowDimension(), 0);

        return newMatrix;               
    }

    /**
     * Calculate the empirical power for the design
     * @param rejectionCriteria
     * @return
     */
    public double calculateEmpiricalPower(LinearHypothesis glh,
            int replicates, int realizations) {

        // create a multivariate normal object to generate replicates of the 
        // error matrix
        double[] errorMeansList = new double[this.sigmaE.getRowDimension()];
        for(int i = 0; i < this.sigmaE.getRowDimension(); i++) { errorMeansList[i] = 0; }
        MultivariateNormalDistribution errorMVN = 
                new MultivariateNormalDistribution(errorMeansList, this.sigmaE.getData());

        // create a multivariate normal object to generate replicates of the covariates
        double[] covarMeansList = new double[this.sigmaG.getRowDimension()];
        for(int i = 0; i < this.sigmaG.getRowDimension(); i++) { covarMeansList[i] = 0; }
        MultivariateNormalDistribution covariateMVN = 
                new MultivariateNormalDistribution(covarMeansList, this.sigmaG.getData());

        // allocate memory for the X matrix
        RealMatrix XMatrix = new Array2DRowRealMatrix(this.XFixed.getRowDimension(),
                this.XFixed.getColumnDimension() + this.sigmaG.getColumnDimension());
        // set the sub-matrix with the fixed portion of the design
        XMatrix.setSubMatrix(this.XFixed.getData(), 0, 0);
        // allocate space for the random matrix
        RealMatrix XRandom = new Array2DRowRealMatrix(this.XFixed.getRowDimension(),
                this.sigmaG.getColumnDimension());

        // pre-allocate an error matrix to save memory
        RealMatrix errorMatrix = new Array2DRowRealMatrix(this.XFixed.getRowDimension(),
                this.sigmaY.getColumnDimension());

        // allocate a memory block for the power values
        double[] powerValues = new double[realizations];

        // generate realizations of the X matrix
        for(int i = 0; i < realizations; i++) {

            int rejectionCount = 0;

            // get X matrix realization
            for(int row = 0; row < XMatrix.getRowDimension(); row++) {
                XRandom.setRow(row, covariateMVN.sample());
            }
            XMatrix.setSubMatrix(XRandom.getData(), 0, this.XFixed.getColumnDimension());

            // calculate X'X inverse and the matrix rank
            this.XtXInverse = new LUDecomposition(XMatrix.transpose().multiply(XMatrix)).
                    getSolver().getInverse();
            int rank = new SingularValueDecomposition(XMatrix).getRank();

            // generate replicates of the error matrices
            for(int rep = 0; rep < replicates; rep++) {

                // generate a data set replicate
                for(int row = 0; row < errorMatrix.getRowDimension(); row++) {
                    errorMatrix.setRow(row, errorMVN.sample());
                }
                RealMatrix YMatrix = XMatrix.multiply(this.betaFull).add(errorMatrix);

                // fit the model and test the hypothesis
                LinearModelFit lmFit = new LinearModelFit(this.XtXInverse, XMatrix, YMatrix, rank);
                LinearHypothesisTestResult result = glh.testHypothesis(lmFit);
                if (result.getPValue() < glh.getAlpha()) {
                    rejectionCount++;
                }
            }

            // empirical power for a given realization is the percentage of rejections
            powerValues[i] = (double) rejectionCount / (double) replicates;

        }

        // obtain empirical unconditional power by averaging across the realizations
        // of X
        return (StatUtils.mean(powerValues));
    }
    
    /**
     * Calculate unconditional power for designs with a single Gaussian covariate
     * @param glh
     * @return
     * @throws PowerException
     * @throws IllegalArgumentException
     */
    public double calculateUnconditionalSingleCovariatePower(LinearHypothesis glh) 
            throws PowerException, IllegalArgumentException {
        // convert the design and hypothesis into the GLMMPowerParameters
        // object
        GLMMPowerParameters params = new GLMMPowerParameters();
        // set matrices
        params.setDesignEssence(XEssence);
        params.setBeta(new FixedRandomMatrix(betaFixed.getData(), betaRandom.getData(), false));

        // break the between participant contrast into fixed and random portions
        double[][] betweenContrastFixed = new double[glh.getBetweenContrastMatrix().getRowDimension()][betaFixed.getRowDimension()];
        glh.getBetweenContrastMatrix().copySubMatrix(0, glh.getBetweenContrastMatrix().getRowDimension()-1,
        		0, betaFixed.getRowDimension()-1, betweenContrastFixed);
        double[][] betweenContrastRandom = 
        		new double[glh.getBetweenContrastMatrix().getRowDimension()][glh.getBetweenContrastMatrix().getColumnDimension() - betaFixed.getRowDimension()];
        glh.getBetweenContrastMatrix().copySubMatrix(0, glh.getBetweenContrastMatrix().getRowDimension()-1,
        		betaFixed.getRowDimension(), glh.getBetweenContrastMatrix().getColumnDimension()-1, betweenContrastRandom);
        params.setBetweenSubjectContrast(new FixedRandomMatrix(betweenContrastFixed, betweenContrastRandom, true));
        
        params.setWithinSubjectContrast(glh.getWithinContrastMatrix());
        params.setTheta(glh.getNullHypothesisMatrix());
        params.setSigmaGaussianRandom(sigmaG);
        params.setSigmaOutcome(sigmaY);
        params.setSigmaOutcomeGaussianRandom(sigmaYG);
        // set lists
        params.addAlpha(glh.getAlpha());
        params.addBetaScale(1);
        params.addSigmaScale(1);
        params.addSampleSize(perGroupN);
        params.addTest(GLMMTestFactory.Test.HOTELLING_LAWLEY_TRACE);
        params.addPowerMethod(GLMMPowerParameters.PowerMethod.UNCONDITIONAL_POWER);

        // calculate unconditional power
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        List<Power> powerList = calc.getPower(params);
        return powerList.get(0).getActualPower();

    }

}


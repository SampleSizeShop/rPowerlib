/*
 *
 *  Package rPowerlib calculates power for the general linear 
 *  multivariate model with and without Gaussian covariates
 *  Copyright (C) 2014 University of Colorado Denver.
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */
package com.kreidles;

import java.io.IOException;

import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

public class LinearHypothesis {
    
    private static final String alphaJSONField = "alpha";
    private static final String testJSONField = "test";
    private static final String nullHypothesisJSONField = "nullHypothesis";
    private static final String betweenContrastJSONField = "betweenContrast";
    private static final String withinContrastJSONField = "withinContrast";

    private static final String hotellingLawleyTrace = "Hotelling-Lawley";
    private static final String pillaiBartlettTrace = "Pillai";
    private static final String wilksLambda = "Wilks";
    
    // the between participant contrast
    RealMatrix betweenContrastMatrix;
    
    // the within participant contrast
    RealMatrix withinContrastMatrix;
    
    // the null hypothesis matrix
    RealMatrix nullHypothesisMatrix;
    
    // the alpha level 
    double alpha;
    
    // the statistical test
    String test;
    
    /**
     * Horizontally append two matrices
     * @param matrix
     * @param column
     * @return the combined matrix
     * @throws IllegalArgumentException
     */
    public static RealMatrix getHorizontalAppend(RealMatrix m1, RealMatrix m2)
    throws IllegalArgumentException
    {
        if (m1 == null || m2 == null)
            throw new IllegalArgumentException("Missing required argument");
        if (m1.getRowDimension() != m2.getRowDimension())
            throw new IllegalArgumentException("Row dimensions must be equal");
        
        RealMatrix newMatrix = 
            new Array2DRowRealMatrix(m1.getRowDimension(),
                    m1.getColumnDimension()+m2.getColumnDimension());
        newMatrix.setSubMatrix(m1.getData(), 0, 0);
        newMatrix.setSubMatrix(m2.getData(), 0, m1.getColumnDimension());
        
        return newMatrix;       
    }
    
    public LinearHypothesis(String rejectionCriteriaAsJSON) 
    throws IOException, ParseException, IllegalArgumentException {
        
        JSONParser jsonParser = new JSONParser();
        JSONObject jsonObject = (JSONObject) jsonParser.parse(rejectionCriteriaAsJSON);
        
        // get the alpha value
        
        this.alpha = Double.parseDouble(jsonObject.get(alphaJSONField).toString());
        if (this.alpha <= 0 || this.alpha >= 1) {
            throw new IllegalArgumentException("Invalid alpha: must be between 0 and 1");
        }
        
        // get the test name
        this.test = (String) jsonObject.get(testJSONField);
        
        // get theta null
        this.nullHypothesisMatrix = 
                JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get(nullHypothesisJSONField));

        // get within contrast
        this.withinContrastMatrix = 
                JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get(withinContrastJSONField));
        if (this.withinContrastMatrix == null) {
            throw new IllegalArgumentException("No within participant contrast specified");
        }
        
        // get fixed and random portions of the between contrast
        this.betweenContrastMatrix = 
                JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get(betweenContrastJSONField));
        if (this.betweenContrastMatrix == null) {
            throw new IllegalArgumentException("No between participant contrast specified");
        }
        
        // if theta null is not specified, create one that conforms with the between and
        // within participant contrasts
        if (this.nullHypothesisMatrix == null) {
            int a = this.betweenContrastMatrix.getRowDimension();
            int b = this.withinContrastMatrix.getColumnDimension();
            this.nullHypothesisMatrix = new Array2DRowRealMatrix(a,b);
            for(int row = 0; row < a; row++) {
                for(int col = 0; col < b; col++) {
                    this.nullHypothesisMatrix.setEntry(row, col, 0);
                }
            }
        } else {
            // if theta null is specified, make sure it conforms with the contrasts
            if (this.nullHypothesisMatrix.getRowDimension() != 
                    this.betweenContrastMatrix.getRowDimension() ||
                    this.nullHypothesisMatrix.getColumnDimension() != 
                    this.withinContrastMatrix.getColumnDimension()) {
                throw new IllegalArgumentException("Null hypothesis does not conform with contrast matrices");
            }
        }

    }

    public RealMatrix getBetweenContrastMatrix() {
        return betweenContrastMatrix;
    }

    public void setBetweenContrastMatrix(RealMatrix betweenContrastMatrix) {
        this.betweenContrastMatrix = betweenContrastMatrix;
    }

    public RealMatrix getWithinContrastMatrix() {
        return withinContrastMatrix;
    }

    public void setWithinContrastMatrix(RealMatrix withinContrastMatrix) {
        this.withinContrastMatrix = withinContrastMatrix;
    }

    public RealMatrix getNullHypothesisMatrix() {
        return nullHypothesisMatrix;
    }

    public void setNullHypothesisMatrix(RealMatrix nullHypothesisMatrix) {
        this.nullHypothesisMatrix = nullHypothesisMatrix;
    }

    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public String getTest() {
        return test;
    }

    public void setTest(String test) {
        this.test = test;
    }
    
    /**
     * Calculate the sum of squares hypothesis matrix (the H matrix)
     * @param params matrices input by user
     * @return H matrix
     */
    protected RealMatrix getHypothesisSumOfSquares(RealMatrix C, RealMatrix beta, RealMatrix U,
            RealMatrix thetaNull, RealMatrix XtXInverse)
    {        
        // the M matrix: [C(X'X)-1C']-1
        RealMatrix M = new LUDecomposition(C.multiply(XtXInverse.multiply(C.transpose()))).getSolver().getInverse();
        // thetaHat = C * Beta * U
        RealMatrix thetaHat = C.multiply(beta.multiply(U));
        // thetaHat - thetaNull.  Multiple by negative one to do subtraction
        RealMatrix thetaDiff = thetaHat.subtract(thetaNull);
        
        // calculate the hypothesis sum of squares: (thetaHat - thetaNull)'[C(X'X)-1C'](thetaHat - thetaNull)
        RealMatrix hss = thetaDiff.transpose().multiply(M.multiply(thetaDiff));
        
        return hss;
        
    }
    
    /**
     * Calculate the sum of squares error matrix (the E matrix)
     * 
     * @param params matrices input by the user
     * @return error sum of squares
     */
    protected RealMatrix getErrorSumOfSquares(RealMatrix U, RealMatrix sigmaError,
            double nuError)
    {        
        return U.transpose().multiply(sigmaError.multiply(U)).scalarMultiply(nuError);
    }    
    
    /**
     * Test the specified hypothesis using the Hotelling Lawley Trace test
     * @param SSPH
     * @param SSPE
     * @param modelFit
     * @return linear hypothesis test result
     */
    protected LinearHypothesisTestResult testHotellingLawleyTrace(RealMatrix SSPH, 
            RealMatrix SSPE, LinearModelFit modelFit) {
        
        // calculate nuE
        double nuE = modelFit.getTotalN() - modelFit.getRankX();
        
        // calculate the test statistic
        RealMatrix inverseE = new LUDecomposition(SSPE).getSolver().getInverse();
        RealMatrix HinverseE = SSPH.multiply(inverseE);
        double hlt = HinverseE.getTrace();
                
        // calculate the numerator degrees of freedom
        double a = this.betweenContrastMatrix.getRowDimension();
        double b = this.withinContrastMatrix.getColumnDimension();
        double ndf = a*b;
        
        // calculate the denominator degrees of freedom
        double t1 = nuE * nuE - nuE * (2 * b + 3) + b * (b + 3);
        double t2 = (nuE * (a  + b + 1) - (a + 2 * b + b * b - 1));
        double ddf = 4 + (a * b + 2) * (t1/t2);
        
        // calculate the observed F
        double observedF = hlt * ((nuE-b-1)*ddf) / (ndf*(ddf-2));
        
        // get the p-Value
        FDistribution fdist = new FDistribution(ndf, ddf);
        double pValue = 1 - fdist.cumulativeProbability(observedF);
        
        return(new LinearHypothesisTestResult(LinearHypothesis.hotellingLawleyTrace, 
                ndf, ddf, observedF, pValue));
        
    }
    
    /**
     * Test the specified hypothesis using the Wilks lambda test
     * @param SSPH
     * @param SSPE
     * @param modelFit
     * @return linear hypothesis test result
     */
    protected LinearHypothesisTestResult testWilksLambda(RealMatrix SSPH, 
            RealMatrix SSPE, LinearModelFit modelFit) {

        
        // calculate nuE
        double nuE = modelFit.getTotalN() - modelFit.getRankX();
        
        // calculate the test statistic
        RealMatrix T = SSPH.add(SSPE);
        RealMatrix inverseT = new LUDecomposition(T).getSolver().getInverse();
        RealMatrix EinverseT = SSPE.multiply(inverseT);       
        double lambda = new LUDecomposition(EinverseT).getDeterminant();
                
        // calculate the numerator degrees of freedom
        double a = this.betweenContrastMatrix.getRowDimension();
        double b = this.withinContrastMatrix.getColumnDimension();
        double ndf = a*b;
        
        double ddf = Double.NaN;
        double association = Double.NaN;
        // calculate the denominator degrees of freedom
        if (a*a*b*b <= 4) {
            ddf = nuE - b + 1;
            association = 1 - lambda;
        } else {
            double gDenominator = (a*a + b*b - 5);
            if (gDenominator == 0)
                throw new IllegalArgumentException("Within and between subject contrasts " +
                		"yielded divide by zero: row of C=" + a + ", cols of U=" + b);
            double g = Math.sqrt((a*a*b*b - 4) / gDenominator);
            ddf = (g*(nuE - (b - a +1)/2)) - (a*b - 2)/2;
            association = 1 - Math.pow(lambda, 1/g);
        }
        
        // calculate the observed F
        double observedF = ((association) / ndf) / ((1 - association) / ddf);
                
        // get the p-Value
        FDistribution fdist = new FDistribution(ndf, ddf);
        double pValue = 1 - fdist.cumulativeProbability(observedF);
        
        return(new LinearHypothesisTestResult(LinearHypothesis.wilksLambda, 
                ndf, ddf, observedF, pValue));

    }
    
    /**
     * Test the specified hypothesis using the Pillai Bartlett Trace test
     * @param SSPH
     * @param SSPE
     * @param modelFit
     * @return linear hypothesis test result
     */
    protected LinearHypothesisTestResult testPillaiBartlettTrace(RealMatrix SSPH, 
            RealMatrix SSPE, LinearModelFit modelFit) {

        
        // calculate nuE
        double nuE = modelFit.getTotalN() - modelFit.getRankX();
        
        // calculate the test statistic
        RealMatrix inverseE = new LUDecomposition(SSPE).getSolver().getInverse();
        RealMatrix HinverseE = SSPH.multiply(inverseE);
        double hlt = HinverseE.getTrace();
                
        // calculate the numerator degrees of freedom
        double a = this.betweenContrastMatrix.getRowDimension();
        double b = this.withinContrastMatrix.getColumnDimension();
        double ndf = a*b;
        
        // calculate the denominator degrees of freedom
        double t1 = nuE * nuE - nuE * (2 * b + 3) + b * (b + 3);
        double t2 = (nuE * (a  + b + 1) - (a + 2 * b + b * b - 1));
        double ddf = 4 + (a * b + 2) * (t1/t2);
        
        // calculate the observed F
        double observedF = hlt * ((nuE-b-1)*ddf) / (ndf*(ddf-2));
        
        // get the p-Value
        FDistribution fdist = new FDistribution(ndf, ddf);
        double pValue = 1 - fdist.cumulativeProbability(observedF);
        
        return(new LinearHypothesisTestResult(LinearHypothesis.hotellingLawleyTrace, 
                ndf, ddf, observedF, pValue));
    }
     
    public LinearHypothesisTestResult testHypothesis(LinearModelFit modelFit) {
        
        RealMatrix SSPH = getHypothesisSumOfSquares(this.betweenContrastMatrix, 
                modelFit.getBetaHat(), this.withinContrastMatrix,
                this.nullHypothesisMatrix, modelFit.getXtXInverse());
        RealMatrix SSPE = getErrorSumOfSquares(this.withinContrastMatrix, 
                modelFit.getSigmaHat(), modelFit.getTotalN() - modelFit.getRankX());
        
        // note, it's easier to use string here so we don't have to deal with 
        // enums from R calls.  Maybe I will fix this someday.
        if (LinearHypothesis.hotellingLawleyTrace.equals(this.test)) {
            return testHotellingLawleyTrace(SSPH, SSPE, modelFit);
        } else if (LinearHypothesis.wilksLambda.equals(this.test)) {
            return testWilksLambda(SSPH, SSPE, modelFit);
        } else if (LinearHypothesis.pillaiBartlettTrace.equals(this.test)) {
            return testPillaiBartlettTrace(SSPH, SSPE, modelFit);
        } else {
            throw new IllegalArgumentException("Invalid test specified");
        }
    }
}

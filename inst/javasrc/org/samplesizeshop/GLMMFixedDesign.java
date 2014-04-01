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
import java.util.List;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import edu.cudenver.bios.matrix.FixedRandomMatrix;
import edu.cudenver.bios.power.GLMMPowerCalculator;
import edu.cudenver.bios.power.Power;
import edu.cudenver.bios.power.PowerException;
import edu.cudenver.bios.power.glmm.GLMMTestFactory;
import edu.cudenver.bios.power.parameters.GLMMPowerParameters;

public class GLMMFixedDesign {
     
    /* matrix inputs to the power calculation.  */
    // the design essence matrix 
    /* For details please see Muller & Fetterman (2002) "Regression and ANOVA" */
    RealMatrix designEssence = null;
    // caching of X'X inverse and rank since these are order(n^3) operations
    RealMatrix XtXInverse = null;
    int designRank = -1;
    // per group sample size
    int perGroupN = -1;

    // beta matrix of regression coefficients
    RealMatrix beta = null;
    // residual error matrix 
    RealMatrix sigmaError = null;

    /**
     * Create a design from a JSON representation of the study.
     * Allows the design to be passed down from R.
     * 
     * @param designAsJSON
     */
    public GLMMFixedDesign(String designAsJSON) 
            throws IOException, ParseException, IllegalArgumentException {
        JSONParser jsonParser = new JSONParser();
        JSONObject jsonObject = (JSONObject) jsonParser.parse(designAsJSON);

        // get the fixed portion of the design matrix
        designEssence = JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get("XEssence"));
        if (designEssence == null) {
            throw new IllegalArgumentException("No X matrix specified");
        }
        // get the fixed portion of the beta matrix
        beta = JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get("beta"));
        if (beta == null) {
            throw new IllegalArgumentException("No beta matrix specified");
        }
        // get the covariance matrices
        sigmaError = JSONUtils.realMatrixFromJsonObject((JSONObject) jsonObject.get("sigmaError"));
        if (sigmaError == null) {
            throw new IllegalArgumentException("No sigmaError matrix specified");
        }
        
        // per group sample size
        perGroupN = Integer.parseInt(jsonObject.get("perGroupN").toString());
        if (perGroupN < 2) {
            throw new IllegalArgumentException("Invalid per group sample size: must be greater than 1");
        }

        // precalculate X'X inverse
        XtXInverse = new LUDecomposition(designEssence.transpose().multiply(designEssence)).getSolver().getInverse();
        // calculate the rank of X
        designRank = new SingularValueDecomposition(designEssence).getRank();
    }
    
    /**
     * Get the conditional power for this design
     * 
     * @param glh linear hypothesis
     */
    public double calculateConditionalPower(LinearHypothesis glh) 
    throws PowerException {
        // convert the design and hypothesis into the GLMMPowerParameters
        // object
        GLMMPowerParameters params = new GLMMPowerParameters();
        // set matrices
        params.setDesignEssence(designEssence);
        params.setBeta(new FixedRandomMatrix(beta.getData(), null, false));
        params.setBetweenSubjectContrast(new FixedRandomMatrix(glh.getBetweenContrastMatrix().getData(), null, true));
        params.setWithinSubjectContrast(glh.getWithinContrastMatrix());
        params.setTheta(glh.getNullHypothesisMatrix());
        params.setSigmaError(sigmaError);
        // set lists
        params.addAlpha(glh.getAlpha());
        params.addBetaScale(1);
        params.addSigmaScale(1);
        params.addSampleSize(perGroupN);
        params.addPowerMethod(GLMMPowerParameters.PowerMethod.CONDITIONAL_POWER);
        params.addTest(GLMMTestFactory.Test.HOTELLING_LAWLEY_TRACE);
        
        // calculate power
        GLMMPowerCalculator calc = new GLMMPowerCalculator();
        List<Power> powerList = calc.getPower(params);
        if (powerList.size() <= 0) {
            throw new IllegalArgumentException("No power returned");
        }
        return powerList.get(0).getActualPower();
    }

}


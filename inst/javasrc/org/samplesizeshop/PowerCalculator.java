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
package org.samplesizeshop;

import java.io.IOException;

import org.json.simple.parser.ParseException;

import edu.cudenver.bios.power.PowerException;

public class PowerCalculator {
    
    /**
     * Calculate the empirical power for the given design and linear hypothesis
     * @param designAsJSON - json representation of a GLMM with fixed and random predictors
     * @param glhAsJSON - json representation of a general linear hypothesis
     * 
     * @return empirical power
     * 
     * @throws IOException
     * @throws ParseException
     * @throws IllegalArgumentException
     */
    public static double calculateEmpiricalPower(String designAsJSON, String glhAsJSON,
            int realizations, int replicates) 
    throws IOException, ParseException, IllegalArgumentException {

        LinearHypothesis glh = new LinearHypothesis(glhAsJSON);
        GLMMFixedRandomDesign design = new GLMMFixedRandomDesign(designAsJSON);

        return design.calculateEmpiricalPower(glh, replicates, realizations);
    }
    
    /**
     * Calculate the conditional power for the given design and linear hypothesis
     * @param designAsJSON - json representation of a GLMM with fixed predictors only
     * @param glhAsJSON - json representation of a general linear hypothesis
     * 
     * @return conditional power
     * 
     * @throws IOException
     * @throws ParseException
     * @throws IllegalArgumentException
     */
    public static double calculateConditionalPower(String designAsJSON, String glhAsJSON) 
    throws IOException, ParseException, PowerException, IllegalArgumentException {

        LinearHypothesis glh = new LinearHypothesis(glhAsJSON);
        GLMMFixedDesign design = new GLMMFixedDesign(designAsJSON);

        return design.calculateConditionalPower(glh);
    }
    
    /**
     * Calculate the unconditional power for the given design and linear hypothesis.
     * 
     * 
     * @param designAsJSON - json representation of a GLMM with fixed and random predictors
     * @param glhAsJSON - json representation of a general linear hypothesis
     * 
     * @return conditional power
     * 
     * @throws IOException
     * @throws ParseException
     * @throws IllegalArgumentException
     */
    public static double calculateUnconditionalSingleCovariatePower(String designAsJSON, String glhAsJSON) 
    throws IOException, ParseException, PowerException, IllegalArgumentException {

        LinearHypothesis glh = new LinearHypothesis(glhAsJSON);
        GLMMFixedRandomDesign design = new GLMMFixedRandomDesign(designAsJSON);

        return design.calculateUnconditionalSingleCovariatePower(glh);
       
    }
    
}

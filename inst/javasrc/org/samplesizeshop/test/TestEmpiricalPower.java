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
package org.samplesizeshop.test;

import java.io.IOException;

import org.json.simple.parser.ParseException;

import org.samplesizeshop.PowerCalculator;
import org.samplesizeshop.GLMMFixedRandomDesign;
import org.samplesizeshop.LinearHypothesis;

import junit.framework.TestCase;

public class TestEmpiricalPower extends TestCase {

    private static final String glhJSON = "{\"alpha\":\"0.05\",\"test\":\"Hotelling-Lawley\"," +
    		"\"nullHypothesisMatrix\": {\"rows\": 1, \"columns\": 1, \"data\": [[0]]}, " + 
    		"\"betweenContrast\": {\"rows\": 1, \"columns\": 3, \"data\": [[1,-1, 0]]}," +    		
            "\"withinContrast\": {\"rows\": 1, \"columns\": 1, \"data\": [[1]]}}";

    
    
    private String designJSON = "{\"name\": \"\",\"description\": \"\",\"XFixed\": {\"rows\": 20,\"columns\": 2,\"data\": [[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ]]},\"betaFixed\": {\"rows\": 2,\"columns\": 3,\"data\": [[ 1,0,0 ],[ 0,0,0 ]]},\"sigmaY\": {\"rows\": 3,\"columns\": 3,\"data\": [[ 1,0.2,0.2 ],[ 0.2,1,0.2 ],[ 0.2,0.2,1 ]]},\"sigmaG\": {\"rows\": 3,\"columns\": 3,\"data\": [[ 1,0.2,0.2 ],[ 0.2,1,0.2 ],[ 0.2,0.2,1 ]]},\"sigmaYG\": {\"rows\": 3,\"columns\": 3,\"data\": [[ 1,0.2,0.2 ],[ 0.2,1,0.2 ],[ 0.2,0.2,1 ]]}}";
    
    private String designJSON2 = 
    		"{\"name\": \"\",\"description\": \"\",\"XFixed\": {\"rows\": 20,\"columns\": 2,\"data\": [[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 1,0 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ],[ 0,1 ]]},\"betaFixed\": {\"rows\": 2,\"columns\": 3,\"data\": [[ 1,0,0 ],[ 0,0,0 ]]},\"sigmaY\": {\"rows\": 3,\"columns\": 3,\"data\": [[ 1,0.2,0.2 ],[ 0.2,1,0.2 ],[ 0.2,0.2,1 ]]},\"sigmaG\": {\"rows\": 1,\"columns\": 1,\"data\": [[ 1.2 ]]},\"sigmaYG\": {\"rows\": 3,\"columns\": 1,\"data\": [[ 0.2 ],[ 0.4 ],[ 0.5 ]]}}";
    
    private String glhJSON2 = 
    		"{\"alpha\": 0.05,\"test\": \"Hotelling-Lawley\",\"betweenContrast\": {\"rows\": 1,\"columns\": 3,\"data\": [[ 1,-1,0 ]]},\"withinContrast\": {\"rows\": 3,\"columns\": 2,\"data\": [[ 1,1 ],[ -1,0 ],[ 0,-1 ]]},\"thetaNull\": {\"rows\": 1,\"columns\": 2,\"data\": [[ 0,0 ]]}}";
    
    
    public void testParsing() {
        try {
            double power = PowerCalculator.calculateEmpiricalPower(designJSON2, 
                    glhJSON2, 1000, 1000);
//            LinearHypothesis glh = new LinearHypothesis(glhJSON);
//            GLMMFixedRandomDesign design = new GLMMFixedRandomDesign(designJSON);
//            double power = design.calculateEmpiricalPower(glh, 1000, 1000);
            System.out.println("Empirical power: " + power);
        } catch (ParseException pe) {
            System.out.println("Parse exception: " + pe.getMessage());
            fail();
        } catch (IOException ioe) {
            System.out.println("IO exception: " + ioe.getMessage());
            fail();
        } catch (IllegalArgumentException iae) {
            System.out.println("Illegal argument: " + iae.getMessage());
            fail();
        } catch (Exception e) {
            System.out.println("Exception: " + e.getMessage());
            fail();
        }

    }
}

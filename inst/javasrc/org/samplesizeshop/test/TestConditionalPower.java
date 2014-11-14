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

import junit.framework.TestCase;

import org.json.simple.parser.ParseException;

import org.samplesizeshop.PowerCalculator;

public class TestConditionalPower extends TestCase {
    private String designJSON = "{\"name\": \"\",\"description\": \"\"," +
    		"\"XEssence\": {\"rows\": 2,\"columns\": 2,\"data\": [[ 1,0 ],[ 0,1 ]]}," +
    		"\"perGroupN\": 20," +
    		"\"beta\": {\"rows\": 2,\"columns\": 3,\"data\": [[ 1,0,0 ],[ 0,0,0 ]]}," +
    		"\"sigmaError\": {\"rows\": 3,\"columns\": 3,\"data\": [[ 1,0.2,0.2 ],[ 0.2,1,0.2 ],[ 0.2,0.2,1 ]]}}";
    private String glhJSON = "{\"alpha\": 0.05,\"test\": \"Hotelling-Lawley\"," +
    		"\"betweenContrast\": {\"rows\": 1,\"columns\": 2,\"data\": [[ 1,-1 ]]}," +
    		"\"withinContrast\": {\"rows\": 3,\"columns\": 2,\"data\": [[ 1,1 ],[ -1,0 ],[ 0,-1 ]]}," +
    		"\"thetaNull\": {\"rows\": 1,\"columns\": 2,\"data\": [[ 0,0 ]]}}";
    
    
    public void testConditionalPower() {
        try {
            double power = PowerCalculator.calculateConditionalPower(designJSON, 
                    glhJSON);
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

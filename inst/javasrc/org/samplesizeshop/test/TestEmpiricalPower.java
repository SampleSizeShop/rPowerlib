package com.kreidles.test;

import java.io.IOException;

import org.json.simple.parser.ParseException;

import com.kreidles.PowerCalculator;
import com.kreidles.GLMMFixedRandomDesign;
import com.kreidles.LinearHypothesis;

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

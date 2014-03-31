package com.kreidles.test;

import java.io.IOException;

import junit.framework.TestCase;

import org.json.simple.parser.ParseException;

import com.kreidles.PowerCalculator;

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

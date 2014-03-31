package com.kreidles.test;

import java.io.IOException;

import junit.framework.TestCase;

import org.json.simple.parser.ParseException;

import com.kreidles.PowerCalculator;

public class TestUnconditionalSingleCovariatePower extends TestCase {
    private String designJSON = "{\"name\": \"\",\"description\": \"\",\"XEssence\": {\"rows\": 2,\"columns\": 2,\"data\": [[ 1,0 ],[ 0,1 ]]},\"perGroupN\": 10,\"betaFixed\": {\"rows\": 2,\"columns\": 3,\"data\": [[ 1,0,0 ],[ 0,0,0 ]]},\"sigmaY\": {\"rows\": 3,\"columns\": 3,\"data\": [[ 1,0.2,0.2 ],[ 0.2,1,0.2 ],[ 0.2,0.2,1 ]]},\"sigmaG\": {\"rows\": 1,\"columns\": 1,\"data\": [[ 1.2 ]]},\"sigmaYG\": {\"rows\": 3,\"columns\": 1,\"data\": [[ 0.2 ],[ 0.4 ],[ 0.5 ]]}}";
    private String glhJSON = "{\"alpha\": 0.05,\"test\": \"Hotelling-Lawley\",\"betweenContrast\": {\"rows\": 1,\"columns\": 3,\"data\": [[ 1,-1,0 ]]},\"withinContrast\": {\"rows\": 3,\"columns\": 2,\"data\": [[ 1,1 ],[ -1,0 ],[ 0,-1 ]]},\"thetaNull\": {\"rows\": 1,\"columns\": 2,\"data\": [[ 0,0 ]]}}";

    public void testConditionalPower() {
        try {
            double power = PowerCalculator.calculateUnconditionalSingleCovariatePower(designJSON, 
                    glhJSON);
            System.out.println("Unconditional power: " + power);
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

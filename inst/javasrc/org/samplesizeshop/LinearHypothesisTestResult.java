package com.kreidles;

public class LinearHypothesisTestResult {

    String test;
    
    double ndf;
    
    double ddf;
    
    double FValue;
    
    double PValue;

    public LinearHypothesisTestResult(String test, double ndf, double ddf, 
            double FValue, double PValue) {
        this.test = test;
        this.ndf = ndf;
        this.ddf = ddf;
        this.FValue = FValue;
        this.PValue = PValue;
    }
    
    public String getTest() {
        return test;
    }

    public void setTest(String test) {
        this.test = test;
    }

    public double getNdf() {
        return ndf;
    }

    public void setNdf(double ndf) {
        this.ndf = ndf;
    }

    public double getDdf() {
        return ddf;
    }

    public void setDdf(double ddf) {
        this.ddf = ddf;
    }

    public double getFValue() {
        return FValue;
    }

    public void setFValue(double fValue) {
        FValue = fValue;
    }

    public double getPValue() {
        return PValue;
    }

    public void setPValue(double pValue) {
        PValue = pValue;
    }

    
    
}

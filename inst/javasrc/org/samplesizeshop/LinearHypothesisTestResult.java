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

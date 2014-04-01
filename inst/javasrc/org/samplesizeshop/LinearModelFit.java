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

import org.apache.commons.math3.linear.RealMatrix;

public class LinearModelFit {

    // the beta coefficients from the fit
    RealMatrix betaHat;
    // the estimated covariance
    RealMatrix sigmaHat;
    
    // store the XtXInverse matrix for use in hypothesis tests
    RealMatrix XtXInverse;
    // store the rank of X
    int rankX;
    // store the total number of independent sampling units
    int totalN;
    
    /**
     * Fits a linear model for the specified X and Y matrices.
     * 
     * @param XtXInverse
     * @param XMatrix
     * @param YMatrix
     */
    public LinearModelFit(RealMatrix XtXInverse, RealMatrix XMatrix, RealMatrix YMatrix, int rankX) {
        
        // store some info about the design
        this.XtXInverse = XtXInverse;
        this.rankX = rankX;
        this.totalN = XMatrix.getRowDimension();
        
        // calculate the regression coefficients
        this.betaHat = XtXInverse.multiply(XMatrix.transpose()).multiply(YMatrix);
        
        // estimate the covariance
        RealMatrix yHat = XMatrix.multiply(this.betaHat);
        RealMatrix yDiff = YMatrix.subtract(yHat);
        this.sigmaHat = yDiff.transpose().multiply(yDiff).scalarMultiply((double) 1/ (double) (this.totalN - rankX));
    }
    
    /**
     * Get the estimated beta coefficients
     * @return
     */
    public RealMatrix getBetaHat() {
        return betaHat;
    }

    /**
     * Get the estimated covariance matrix
     * @return
     */
    public RealMatrix getSigmaHat() {
        return sigmaHat;
    }

    /**
     * Get X'X inverse matrix
     * @return
     */
    public RealMatrix getXtXInverse() {
        return XtXInverse;
    }

    /**
     * Get the rank of the X matrix for this model 
     * @return
     */
    public int getRankX() {
        return rankX;
    }
    
    /**
     * Get the total number of sampling units for this model 
     * @return
     */
    public int getTotalN() {
        return totalN;
    }
    
    
}

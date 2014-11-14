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

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Matrix of random normals with the specified covariance matrix
 */
public class GaussianErrorMatrix
{
	private double symmetryThreshold = 
	    CholeskyDecomposition.DEFAULT_RELATIVE_SYMMETRY_THRESHOLD;
	private double positivityThreshold = 
		CholeskyDecomposition.DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD;
	
	protected long seed = 1234;
	protected RealMatrix matrix;
	protected RealMatrix sigma;
	protected CholeskyDecomposition cholesky = null;
	protected RealMatrix sqrtMatrix = null;
	protected NormalDistribution normalDist;

	public GaussianErrorMatrix(int rows, int cols, RealMatrix sigma)
	{
		this.matrix = new Array2DRowRealMatrix(rows, cols);
		this.sigma = sigma;
		this.normalDist = new NormalDistribution();
		this.normalDist.reseedRandomGenerator(this.seed);
	}

	public void setSeed(long seed)
	{
		this.seed = seed;
		normalDist.reseedRandomGenerator(seed);
	}
	
    /**
     * Set the positivity threshold for Cholesky decomposition.  This
     * allows Cholesky decomposition for matrices with very small negative values.
     */
	public void setPositivityThreshold(double positivityThreshold)
	{
		this.positivityThreshold = positivityThreshold;
		this.cholesky = null;
	}
	
    /**
     * Set the symmetric threshold for Cholesky decomposition.  This
     * allows Cholesky decomposition for matrices with very small differences
     * between symmetric cells.
     */
	public void setSymmetryThreshold(double symmetryThreshold)
	{
		this.symmetryThreshold = symmetryThreshold;
		this.cholesky = null;
	}
	
    /**
     * Simulate the error matrix in the Y = X * beta + e
     * 
     * @param normalDist normal distribution object for generating random samples
     * @param error matrix to hold random normal samples
     * @param rows number of rows in the random normal samples matrix
     * @param columns number of columns in the random normal samples matrix
     * @param sigma the covariance matrix for the error term
     * @return a random instance of the 'e' matrix in the model
     */
	public RealMatrix random()
	throws IllegalArgumentException
	{
        // build a matrix of random values from a standard normal
        // the number of rows = #subjects (rows) in the full design matrix
        // the number of columns = #outcome variables (i.e. columns in beta)
        for(int rowIndex = 0; rowIndex < matrix.getRowDimension(); rowIndex++)
        {
            for(int columnIndex = 0; columnIndex < matrix.getColumnDimension(); columnIndex++)
            {
            	matrix.setEntry(rowIndex, columnIndex, normalDist.sample()); 
            }
        }
        
        // take the square root of the sigma matrix via cholesky decomposition
        try
        {            
        	if (this.cholesky == null)
        	{
        		this.cholesky = new CholeskyDecomposition(sigma, 
                        symmetryThreshold,
                        positivityThreshold);
                sqrtMatrix = cholesky.getLT();
        	}
            return matrix.multiply(sqrtMatrix); 
        }
        catch (Exception e)
        {
            throw new IllegalArgumentException(e);
        }
	}

}

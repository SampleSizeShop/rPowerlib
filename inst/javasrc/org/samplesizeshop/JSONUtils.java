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

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

public class JSONUtils {
    
    /**
     * Read a matrix from a JSON object
     * @param jsonObj
     * @return
     */
    public static RealMatrix realMatrixFromJsonObject(JSONObject jsonObj) {
        Array2DRowRealMatrix matrix = null;
        if (jsonObj != null) {
            int rows = ((Long) jsonObj.get("rows")).intValue();
            int columns = ((Long) jsonObj.get("columns")).intValue();
            
            if (rows > 0 && columns > 0) {
                // allocate the matrix
                matrix = new Array2DRowRealMatrix(rows, columns);
                JSONArray data = (JSONArray) jsonObj.get("data");
                for(int row = 0; row < rows; row++) {
                    JSONArray dataRow = (JSONArray) data.get(row);
                    for(int col = 0; col < columns; col++) {
                        matrix.setEntry(row, col, Double.parseDouble(dataRow.get(col).toString()));                    
                    }
                }
            }
        }
        return matrix;
    }
}

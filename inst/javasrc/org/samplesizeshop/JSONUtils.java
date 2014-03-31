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

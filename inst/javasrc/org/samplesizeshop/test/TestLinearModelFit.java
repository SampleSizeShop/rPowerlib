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
package com.kreidles.test;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import com.kreidles.JSONUtils;
import com.kreidles.LinearHypothesis;
import com.kreidles.LinearHypothesisTestResult;
import com.kreidles.LinearModelFit;

import junit.framework.TestCase;

public class TestLinearModelFit  extends TestCase {
    private static final String XMatrixJSON = "{\"rows\": 20,\"columns\": 5,\"data\": " +
    		"[[ 1,0,0.944138201678652,-0.533421994093656,0.363046523883714 ]," +
    		"[ 0,1,-0.522222561936891,1.17063662173391,1.44574276139972 ]," +
    		"[ 1,0,1.0482893640533,-0.15553609216374,-1.29299836004313 ]," +
    		"[ 0,1,-0.37955304385142,-1.35725375842286,1.63005368086564 ]," +
    		"[ 1,0,1.06148997239156,0.169649022647193,0.60555104312437 ]," +
    		"[ 0,1,0.769218242581803,0.294405872388656,-0.110532833610114 ]," +
    		"[ 1,0,0.0427469079132732,-2.90986877642217,-0.388902871645494 ]," +
    		"[ 0,1,-0.550237315455219,-0.556399961141128,1.59583428462821 ]," +
    		"[ 1,0,-0.774400422186136,-0.691140676933056,-0.0799403970168287 ]," +
    		"[ 0,1,-0.156264163430493,-0.35749990761085,0.382518040411337 ]," +
    		"[ 1,0,0.155855760566286,0.0868367176510585,-1.52058282452452 ]," +
    		"[ 0,1,-1.92564105290871,-0.908094664162244,-1.21837121513043 ]," +
    		"[ 1,0,0.791379143846431,-0.397536522380725,-0.617245828101027 ]," +
    		"[ 0,1,-0.0547750063937297,2.21245213540497,1.1589896922944 ]," +
    		"[ 1,0,-0.941225024031432,0.149055739987833,-1.54182230489007 ]," +
    		"[ 0,1,1.32994553683777,1.0618794707451,0.686720348575536 ]," +
    		"[ 1,0,1.4705384371699,0.697514113743949,1.79730394401647 ]," +
    		"[ 0,1,-2.0762357740981,-0.584129409779393,0.319691342520041 ]," +
    		"[ 1,0,-0.220353886695309,0.575474541798893,2.09836275838748 ]," +
    		"[ 0,1,1.57832547493123,2.24523550211472,1.13465428921809 ]]}";
    
    private static final String YMatrixJSON = "{\"rows\": 20,\"columns\": 3,\"data\": " +
    		"[[ 10.5382263890912,0.325515067883075,-0.765446243355672 ]," +
    		"[ 1.03914925046587,2.43011976587977,0.0579213538549594 ]," +
    		"[ 9.47959854028668,-0.271046432267618,-0.357746486750807 ]," +
    		"[ 0.279757783308449,-0.658731645798706,0.754895586373453 ]," +
    		"[ 10.1232512479531,2.1163390842408,-0.383797999331926 ]," +
    		"[ 0.528169665787211,0.837905711623569,-0.401161698484598 ]," +
    		"[ 8.53320313834908,-1.88352411762637,-0.515149460250139 ]," +
    		"[ -0.0473814953760626,1.21687055128431,-0.203765026222487 ]," +
    		"[ 9.27076333416819,0.764077946703127,-0.1821437220054 ]," +
    		"[ 0.282764488949744,-0.152118295425824,0.936794807121796 ]," +
    		"[ 7.56848905088999,-1.46370038686487,1.2042014371056 ]," +
    		"[ -2.15180404003024,-0.141746136445409,-1.87108387169488 ]," +
    		"[ 8.89663479042918,0.198399614818792,-0.279589715934205 ]," +
    		"[ 2.44257854554352,1.78737481730465,-0.427536745965828 ]," +
    		"[ 9.15638769814498,-0.550417695472878,-0.667944209517552 ]," +
    		"[ 0.999957726595971,1.10303134741868,1.17335963394059 ]," +
    		"[ 11.9794330761105,0.918895975290666,1.37767612203888 ]," +
    		"[ 0.0209727613601934,-0.761333751132666,-0.627174324765553 ]," +
    		"[ 10.2385048338779,1.7908909726996,1.46619577839758 ]," +
    		"[ 1.31418783223407,1.24191273149443,2.97916165094121 ]]}";
    
    private static final String glhJSON = "{\"alpha\": 0.05,\"test\": \"Hotelling-Lawley\"," +
    		"\"betweenContrast\": {\"rows\": 1,\"columns\": 5,\"data\": " +
    		"[[ 1,-1,0,0,0 ]]},\"withinContrast\": " +
    		"{\"rows\": 3,\"columns\": 2,\"data\": [[ 1,1 ],[ -1,0 ],[ 0,-1 ]]}," +
    		"\"thetaNull\": {\"rows\": 1,\"columns\": 2,\"data\": [[ 0,0 ]]}}";
    
    public void testModelFit() {
        try {
            JSONParser jsonParser = new JSONParser();
            JSONObject XJsonObject = (JSONObject) jsonParser.parse(XMatrixJSON);
            JSONObject YJsonObject = (JSONObject) jsonParser.parse(YMatrixJSON);
            
            RealMatrix XMatrix = JSONUtils.realMatrixFromJsonObject(XJsonObject);
            RealMatrix YMatrix = JSONUtils.realMatrixFromJsonObject(YJsonObject); 
            
            RealMatrix XtXInverse = 
                    new LUDecomposition(XMatrix.transpose().
                            multiply(XMatrix)).getSolver().getInverse();
        
            int rank = new SingularValueDecomposition(XMatrix).getRank();
            
            LinearModelFit fit = new LinearModelFit(XtXInverse, XMatrix, YMatrix, rank);
            
            LinearHypothesis hypothesis = new LinearHypothesis(glhJSON);
            LinearHypothesisTestResult result = hypothesis.testHypothesis(fit);
            
            System.out.println("Done");
        } catch (Exception e) {
            System.out.println("Parsing failed: " + e.getMessage());
            fail();
        }

        
    }

}

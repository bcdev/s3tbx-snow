package org.esa.s3tbx.snow;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.CurveFitter;
import org.apache.commons.math3.optim.SimpleVectorValueChecker;
import org.apache.commons.math3.optim.nonlinear.vector.MultivariateVectorOptimizer;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.GaussNewtonOptimizer;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;

/**
 * todo: add comment
 * To change this template use File | Settings | File Templates.
 * Date: 25.09.2017
 * Time: 13:58
 *
 * @author olafd
 */
public class SpectralAlbedoFuncFitter extends CurveFitter {
    public SpectralAlbedoFuncFitter() {
//        super(new LevenbergMarquardtOptimizer());
        super(new GaussNewtonOptimizer(new SimpleVectorValueChecker(1.0, 1.0)));
    }

    public static class SigmoidFunc implements ParametricUnivariateFunction {
        public double value(double t, double... parms) {
            // AK's fitting function:
            return parms[0] + parms[1] / (1.0 + Math.exp((t-parms[2])/parms[3]));
        }

        // Jacobian matrix of the above. In this case, this is just an array of
        // partial derivatives of the above function, with one element for each parameter.
        public double[] gradient(double t, double... parms) {
            final double aDeriv = 1.0;
            final double EXP = Math.exp((t - parms[2]) / parms[3]);
            final double bDeriv = 1.0 / (1.0 + EXP);
            final double cDeriv = EXP/(parms[3]*(1+EXP)*(1+EXP));
            final double dDeriv = (parms[3] - t)*EXP/(parms[3]*(1+EXP)*(1+EXP));

            return new double[]{aDeriv, bDeriv, cDeriv, dDeriv};
        }
    }
}



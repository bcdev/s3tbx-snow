package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;

/**
 * Implementation of Sigmoidal function
 *
 */
class SigmoidalFunction  {

    private static double evaluate(double[] coeffs, double x) throws NullArgumentException, NoDataException {
        MathUtils.checkNotNull(coeffs);
        int n = coeffs.length;
        // approach with 4 parameters: f(x) := A + B/(1+exp(C*x + D))
//        if(n != 4) {
//            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
//        } else {
//            // AK's fitting function:
//            // original:
//            // return coeffs[0] + coeffs[1] / (1.0 + Math.exp((x-coeffs[2])/coeffs[3]));
//
//            // set the coeffs this way:
//            return coeffs[0] + coeffs[1] / (1.0 + Math.exp(coeffs[2]*x + coeffs[3]));
//        }

        // approach with just 2 parameters: f(x) := 1/(1+exp(A*x + B))
        if(n != 2) {
            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
        } else {
            // AK's fitting function:
            // original:
            // return coeffs[0] + coeffs[1] / (1.0 + Math.exp((x-coeffs[2])/coeffs[3]));

            // set the coeffs this way:
            return 1.0 / (1.0 + Math.exp(coeffs[0]*x + coeffs[1]));
        }
    }

    static class Parametric implements ParametricUnivariateFunction {
        Parametric() {
        }

        public double[] gradient(double x, double... parms) {
            // 4 parameter approach:
//            final double aDeriv = 1.0;
//            final double EXP = Math.exp(parms[2]*x + parms[3]);
//            final double bDeriv = 1.0 / (1.0 + EXP);
//            final double cDeriv = -parms[1]*x*EXP/((1+EXP)*(1+EXP));
//            final double dDeriv = -parms[1]*EXP/((1+EXP)*(1+EXP));
//
//            return new double[]{aDeriv, bDeriv, cDeriv, dDeriv};

            // 2 parameter approach:
            final double EXP = Math.exp(parms[0]*x + parms[1]);
            final double cDeriv = -x*EXP/((1+EXP)*(1+EXP));
            final double dDeriv = -EXP/((1+EXP)*(1+EXP));

            return new double[]{cDeriv, dDeriv};
        }

        public double value(double x, double... parameters) throws NoDataException {
            return SigmoidalFunction.evaluate(parameters, x);
        }
    }
}

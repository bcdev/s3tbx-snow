package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;

/**
 * Implementation of exponential function with 4 parameters as:
 * f(x; a, b, kappa1, l) := a *(exp(-(2*PI/x) * (kappa1*L + kappa2*L)))^b
 *
 * @author olafd
 */
public class Exp4ParamFunction implements ParametricUnivariateFunction {

    public Exp4ParamFunction() {
    }

    public double[] gradient(double x, double... parms) {
        // 4 parameter approach:
        final double bracket = Math.exp(-(parms[2] + parms[3] / x));
        final double aDeriv = Math.pow(bracket, parms[3]);
        final double bDeriv = -parms[0] * parms[3] * Math.exp(-parms[2] * parms[3] / x) * Math.exp(-parms[1] * parms[3]);
        final double cDeriv = bDeriv / x;
        final double dDeriv = Math.log(parms[0]*bracket) * Math.exp(parms[3]*Math.log(parms[0]*bracket));

        return new double[]{aDeriv, bDeriv, cDeriv, dDeriv};
    }

    public double value(double x, double... parameters) throws NoDataException {
        return Exp4ParamFunction.evaluate(parameters, x);
    }

    private static double evaluate(double[] parms, double x) throws NullArgumentException, NoDataException {
        MathUtils.checkNotNull(parms);
        final int numParms = parms.length;

        if (numParms == 4) {
            // approach with 4 parameters: f(x; a, b, kappa1, L) := a *(exp(- ( (2*PI/x) * (kappa1*L + kappa2*L) ))))^b
            final double a = parms[0];
            final double kappa1 = parms[1];
            final double L = parms[2];
            final double b = parms[3];

            // todo provide kappa2 from spline
            return 0;
//            return a * Math.pow(Math.exp(- ( (2.*Math.PI/x) * (kappa1*L + kappa2*L) )), b);
        } else {
            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
        }
    }

}

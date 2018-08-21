package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;
import org.esa.s3tbx.snow.OlciSnowPropertiesConstants;

/**
 * Implementation of exponential function with 4 parameters as:
 * f(v; a, kappa1, L, b) := a *(exp(-(2*PI/v) * (kappa1*L + kappa2(v)*L)))^b
 *
 * @author olafd
 */
public class Exp4ParamFunction implements ParametricUnivariateFunction {

    public double[] gradient(double v, double... parms) {
        // 4 parameter approach:

        double[] kappaCoeffs = getKappaCoeffs(v);
        final PolynomialFunction kappa2Function = new PolynomialFunction(kappaCoeffs);
        final double kappa2 = kappa2Function.value(v);

        final double bracket = Math.exp(-(2 * Math.PI / v) * (parms[1] * parms[2] + kappa2 * parms[2]));
        final double aDeriv = Math.pow(bracket, parms[3]);
        final double kappa1Deriv = (-2 * Math.PI * parms[0] * parms[2] * parms[3] / v) *
                Math.exp(-2 * Math.PI * parms[2] * parms[3] *(parms[1] + kappa2) / v);
        final double lDeriv = (-2 * Math.PI * parms[0] * parms[2] * parms[3] *(parms[1] + kappa2) / v) *
                Math.exp(-2 * Math.PI * parms[2] * parms[3] *(parms[1] + kappa2) / v);
        final double bDeriv = parms[0] * Math.log( bracket) * Math.exp(parms[3] * Math.log(bracket));

        return new double[]{aDeriv, kappa1Deriv, lDeriv, bDeriv};
    }

    public double value(double x, double... parameters) throws NoDataException {
        return evaluate(parameters, x);
    }

    private static double evaluate(double[] parms, double v) throws NullArgumentException, NoDataException {
        MathUtils.checkNotNull(parms);
        final int numParms = parms.length;

        if (numParms == 4) {
            // approach with 4 parameters: f(v; a, kappa1, L, b) := a *(exp(- ( (2*PI/v) * (kappa1*L + kappa2(v)*L) ))))^b
            final double a = parms[0];
            final double kappa1 = parms[1];
            final double L = parms[2];
            final double b = parms[3];

            double[] kappaCoeffs = getKappaCoeffs(v);
            final PolynomialFunction kappa2Function = new PolynomialFunction(kappaCoeffs);
            final double kappa2 = kappa2Function.value(v);

            return a * Math.pow(Math.exp(-((2. * Math.PI / v) * (kappa1 * L + kappa2 * L))), b);
        } else {
            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
        }
    }

    private static double[] getKappaCoeffs(double x) {
        if (x >= 0.4 && x < 0.525) {
            return OlciSnowPropertiesConstants.KAPPA_2_SPLINE_COEFFS[0];
        } else if (x >= 0.525 && x < 0.7) {
            return OlciSnowPropertiesConstants.KAPPA_2_SPLINE_COEFFS[1];
        } else if (x >= 0.7 && x <= 1.02) {
            return OlciSnowPropertiesConstants.KAPPA_2_SPLINE_COEFFS[2];
        } else {
            throw new IllegalArgumentException("x must be in range [0.4, 1.02");
        }
    }

}

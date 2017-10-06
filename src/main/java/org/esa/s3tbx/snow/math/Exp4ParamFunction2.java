package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;
import org.esa.s3tbx.snow.OlciSnowAlbedoConstants;

/**
 * Implementation of exponential function with 3 parameters as:
 * f(v; a, b, c) := a * exp(-sqrt((b + c*kappa2(v))/v))
 *
 * @author olafd
 */
public class Exp4ParamFunction2 implements ParametricUnivariateFunction {

    public double[] gradient(double v, double... parms) {
        // approach with 3 parameters: f(v; a, b, c) := a * exp(-sqrt((b + c*kappa2(v))/v))

        double[] kappaCoeffs = getKappaCoeffs(v);
        final PolynomialFunction kappa2Function = new PolynomialFunction(kappaCoeffs);
        final double kappa2 = kappa2Function.value(v);


        final double a = parms[0];
        final double b = parms[1];
        final double c = parms[2];
        final double root = Math.sqrt((b + c*kappa2)/v);

        final double aDeriv = Math.exp(-root);
        final double bDeriv = -a * root * Math.exp(-root) / (2.*v);
        final double cDeriv = kappa2 * bDeriv;

        return new double[]{aDeriv, bDeriv, cDeriv};
    }

    public double value(double x, double... parameters) throws NoDataException {
        return evaluate(parameters, x);
    }

    protected static double evaluate(double[] parms, double v) throws NullArgumentException, NoDataException {
        MathUtils.checkNotNull(parms);
        final int numParms = parms.length;

        if (numParms == 3) {
            // approach with 3 parameters: f(v; a, b, c) := a * exp(-sqrt((b + c*kappa2(v))/v))
            final double a = parms[0];
            final double b = parms[1];
            final double c = parms[2];

            double[] kappaCoeffs = getKappaCoeffs(v);
            final PolynomialFunction kappa2Function = new PolynomialFunction(kappaCoeffs);
            final double kappa2 = kappa2Function.value(v);

            return a * Math.exp(-Math.sqrt((b + c*kappa2)/v));
        } else {
            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
        }
    }

    private static double[] getKappaCoeffs(double x) {
        if (x >= 0.4 && x < 0.525) {
            return OlciSnowAlbedoConstants.KAPPA_2_SPLINE_COEFFS[0];
        } else if (x >= 0.525 && x < 0.7) {
            return OlciSnowAlbedoConstants.KAPPA_2_SPLINE_COEFFS[1];
        } else if (x >= 0.7 && x <= 1.02) {
            return OlciSnowAlbedoConstants.KAPPA_2_SPLINE_COEFFS[2];
        } else {
            throw new IllegalArgumentException("x must be in range [0.4, 1.02");
        }
    }

}

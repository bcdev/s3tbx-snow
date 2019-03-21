package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;
import org.esa.s3tbx.snow.SnowUtils;

/**
 * Implementation of 'fun1' function as in breadboard sice.f90:
 *
 * @author olafd
 */
public class SiceFun1InterpolInsideFunction implements ParametricUnivariateFunction {

    private double[] xa;
    private double[] ya;

    public SiceFun1InterpolInsideFunction(double[] xa, double[] ya) {
        this.xa = xa;
        this.ya = ya;
    }

    public double[] gradient(double v, double... parms) {
        // not needed
        return null;
    }

    public double value(double x, double... parameters) throws NoDataException {
        MathUtils.checkNotNull(parameters);
        final int numParms = parameters.length;
        if (numParms == 8) {
            return evaluate(parameters, x);
        } else {
            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
        }
    }

    private double evaluate(double[] parms, double x) throws NullArgumentException, NoDataException {
        final double astra = computeAstraFromInterpolation(x);
        return compute(parms, x, astra);
    }

    private double computeAstraFromInterpolation(double x) {
        return SnowUtils.linearInterpolate(x, xa, ya);
    }

    private static double compute(double[] parms, double x, double astra) {
        // params:
//        double brr400, double effAbsLength, double r0a1Thresh, double sza,
//        double as, double bs, double cs, double planar, double astra
        final double brr400 = parms[0];
        final double effAbsLength = parms[1];
        final double r0a1Thresh = parms[2];
        final double cosSza = parms[3];
        final double as = parms[4];
        final double bs = parms[5];
        final double cs = parms[6];
        final double planar = parms[7];

        final double dega = effAbsLength * 4.0 * Math.PI * astra / x;
        final double sqrtDega = Math.sqrt(dega);
        final double rsd = sqrtDega > 1.E-6 ? Math.exp(-sqrtDega) : 1.0;
        final double um1 = SnowUtils.computeU(cosSza);

        double f1;
        if (brr400 <= r0a1Thresh && x <= 1.02) {
            f1 = as * x * x + bs * x + cs;
        } else {
            f1 = planar == 0.0 ? rsd : Math.pow(rsd, um1);
        }

        final double p0 = 32.38;
        final double p1 = -160140.33;
        final double p2 = 7959.53;
        final double t1 = 85.34 * 1.e-3;
        final double t2 = 401.79 * 1.e-3;

        double funcs;
        if (x <= 0.4) {
            funcs = p0 + p1 * Math.exp(-0.4 / t1) + p2 * Math.exp(-0.4 / t2);
        } else {
            funcs = p0 + p1 * Math.exp(-x / t1) + p2 * Math.exp(-x / t2);
        }

        return f1 * funcs;
    }
}

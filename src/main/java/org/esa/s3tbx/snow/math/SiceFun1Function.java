package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;
import org.esa.s3tbx.snow.OlciSnowPropertiesConstants;
import org.esa.s3tbx.snow.SnowUtils;

/**
 * Implementation of 'fun1' function as in breadboard sice.f90:
 *
 * @author olafd
 */
public class SiceFun1Function implements ParametricUnivariateFunction {

    public double[] gradient(double v, double... parms) {
        // not needed
        return null;
    }

    public double value(double x, double... parameters) throws NoDataException {
        MathUtils.checkNotNull(parameters);
        final int numParms = parameters.length;
        if (numParms == 8) {
            return evaluate_old(parameters, x);
        } else if (numParms == 9) {
            return evaluate_new(parameters, x);
        } else {
            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
        }
    }

    private static double evaluate_old(double[] parms, double x) throws NullArgumentException, NoDataException {
        final double astra = computeAstraFromPolynom(x);
        return compute(parms, x, astra);
    }

    private static double evaluate_new(double[] parms, double x) throws NullArgumentException, NoDataException {
        final double astra = parms[parms.length-1];
        return compute(parms, x, astra);
    }

    private static double computeAstraFromPolynom(double x) {
        double[] a = new double[6];
        final double[][] bbb = OlciSnowPropertiesConstants.BBB_COEFFS_SICE;
        final double[] bbbBounds = OlciSnowPropertiesConstants.BBB_COEFFS_SICE_BOUNDS;
        for (int i = 0; i < 9; i++) {
            for (int j = 0; j < 6; j++) {
                if (x >= bbbBounds[i] && x < bbbBounds[i + 1]) {
                    a[j] = bbb[i][j];
                }
            }
        }

        double astra;
        if (x < 0.4) {
            astra = 2.0E-11;
        } else {
            // obviously we lose some precision in the breadboard with the 'real' numbers
            final double term1 = a[0];
            final double term2 = a[1] * x;
            final double term3 = a[2] * x * x;
            final double term4 = a[3] * x * x * x;
            final double term5 = a[4] * x * x * x * x;
            final double term6 = a[5] * x * x * x * x * x;
            astra = term1 + term2 + term3 + term4 + term5 + term6;
        }
        return astra;
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
            f1 = planar == 1.0 ? rsd : Math.pow(rsd, um1);
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

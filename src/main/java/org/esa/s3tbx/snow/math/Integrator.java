package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;

/**
 * Class providing methods for numerical integration.
 *
 * @author olafd
 */
public class Integrator {

    /**
     * Provides 1D integration with Trapezoidal rule with equally spaced discretization.
     * See e.g. https://en.wikipedia.org/wiki/Trapezoidal_rule
     *
     * @param lower - lower grid boundary
     * @param upper - upper grid boundary
     * @param y     - function values
     * @param x     - grid values
     * @return the integral
     */
    static double integrateTrapezoid(double lower, double upper, double[] y, double[] x) {
        final double dx = x[1] - x[0];
        int lowerIndex = (int) ((lower - x[0]) / dx);
        lowerIndex = Math.min(Math.max(lowerIndex, 0), x.length - 1);
        int upperIndex = (int) ((upper - x[0]) / dx);
        upperIndex = Math.max(0, Math.min(upperIndex, x.length - 1));
        final int n = upperIndex - lowerIndex;

        double sum = 0.0;
        for (int i = lowerIndex; i < upperIndex; i++) {
            sum += 2.0 * y[i];
        }
        return (upper - lower) * (y[lowerIndex] + sum + y[upperIndex]) / (2 * n);
    }

    /**
     * Provides 1D integration with Composite Simpson's rule with equally spaced discretization.
     * See e.g. https://en.wikipedia.org/wiki/Simpson%27s_rule
     *
     * @param lower - lower grid boundary
     * @param upper - upper grid boundary
     * @param y     - function values
     * @param x     - grid values
     * @return the integral
     */
    public static double integrateSimpson(double lower, double upper, double[] y, double[] x) {
        final double dx = x[1] - x[0];
        int lowerIndex = (int) ((lower - x[0]) / dx);
        lowerIndex = Math.min(Math.max(lowerIndex, 0), x.length - 1);
        int upperIndex = (int) ((upper - x[0]) / dx);
        upperIndex = Math.max(0, Math.min(upperIndex, x.length - 1));
        final int n = upperIndex - lowerIndex;

        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int i = lowerIndex; i < upperIndex; i += 2) {
            sum1 += 4.0 * y[i];
        }
        for (int i = lowerIndex + 1; i < upperIndex - 1; i += 2) {
            sum2 += 2.0 * y[i];
        }
        return (upper - lower) * (y[lowerIndex] + sum1 + sum2 + y[upperIndex]) / (3 * n);
    }

    public static double integrateSimpsonSiceAlex(double lower, double upper,
                                                  ParametricUnivariateFunction fun,
                                                  double[] params) {

        // todo: jmax is the most critical number for processor overall performance!!
        // The lower the better, but 7 is the minimum to keep BB differences within ~1% (0.01)
//        final int jmax = 7;
        final int jmax = 8;
        double trapezResult = 0.0;
        double s = 0.0;

        double ost = -1.E30;
        double os = -1.E30;

        for (int i = 1; i < jmax; i++) {
            trapezResult = integrateTrapezSiceAlex(lower, upper, fun, params, trapezResult, i);
            s = (4.0 * trapezResult - ost) / 3.0;
            if (i > 5) {
                if (Math.abs(s - os) < 1.E-3 * Math.abs(os) || (s == 0.0 && os == 0.0)) {
                    break;
                }
            }
            os = s;
            ost = trapezResult;
        }

        return s;
    }

    private static double integrateTrapezSiceAlex(double lower, double upper,
                                                  ParametricUnivariateFunction fun,
                                                  double[] params,
                                                  double tmpSum,
                                                  int n) {

        double s;
        if (n == 1) {
            double funLower = fun.value(lower, params);
            double funUpper = fun.value(upper, params);
            s = 0.5 * (upper - lower) * (funLower + funUpper);
        } else {
            int it = (int) Math.pow(2.0, n - 2);
            double del = (upper - lower) / it;
            double sum = 0.0;
            double xx = lower + 0.5 * del;
            for (int i = 0; i < it; i++) {
                sum += fun.value(xx, params);
                xx += del;
            }
            s = 0.5 * (tmpSum + (upper - lower) * sum / it);
        }
        return s;
    }

}

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
    public static double integrateTrapezoid(double lower, double upper, double[] y, double[] x) {
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

    public static double integrateSimpsonSice(double lower, double upper,
                                              ParametricUnivariateFunction fun,
                                              double[] params,
                                              double[] x) {
        final double dx = x[1] - x[0];
        int lowerIndex = (int) ((lower - x[0]) / dx);
        lowerIndex = Math.min(Math.max(lowerIndex, 0), x.length - 1);
        int upperIndex = (int) ((upper - x[0]) / dx);
        upperIndex = Math.max(0, Math.min(upperIndex, x.length - 1));
        final int n = upperIndex - lowerIndex;

        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int i = lowerIndex; i < upperIndex; i += 2) {
            sum1 += 4.0 * fun.value(x[i], params);
        }
        for (int i = lowerIndex + 1; i < upperIndex - 1; i += 2) {
            sum2 += 2.0 * fun.value(x[i], params);
        }
        return (upper - lower) *
                (fun.value(x[lowerIndex], params) + sum1 + sum2 + fun.value(x[upperIndex], params)) / (3 * n);
    }

    public static double integrateSimpsonSiceOlaf(double lower, double upper,
                                                  ParametricUnivariateFunction fun,
                                                  double[] params,
                                                  double[] x,
                                                  double[] astra) {
        final double dx = x[1] - x[0];
        int lowerIndex = (int) ((lower - x[0]) / dx);
        lowerIndex = Math.min(Math.max(lowerIndex, 0), x.length - 1);
        int upperIndex = (int) ((upper - x[0]) / dx);
        upperIndex = Math.max(0, Math.min(upperIndex, x.length - 1));
        final int n = upperIndex - lowerIndex;

        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int i = lowerIndex; i < upperIndex; i += 2) {
            params[params.length - 1] = astra[i];
            sum1 += 4.0 * fun.value(x[i], params);
        }
        for (int i = lowerIndex + 1; i < upperIndex - 1; i += 2) {
            params[params.length - 1] = astra[i];
            sum2 += 2.0 * fun.value(x[i], params);
        }

        params[params.length - 1] = astra[lowerIndex];
        final double funLower = fun.value(x[lowerIndex], params);
        params[params.length - 1] = astra[upperIndex];
        final double funUpper = fun.value(x[upperIndex], params);
        return (upper - lower) * (funLower + sum1 + sum2 + funUpper) / (3 * n);
    }


    public static double integrateSimpsonSiceAlex(double lower, double upper,
                                                  ParametricUnivariateFunction fun,
                                                  double[] params,
                                                  double[] x) {

        // todo: jmax is the most critical number for processor overall performance!!
        // The lower the better, but 7 is the minimum to keep BB differences within ~1% (0.01)
        final int jmax = 7;
        double trapezResult = 0.0;
        double s = 0.0;

        double ost = -1.E30;
        double os = -1.E30;

        for (int i = 1; i < jmax; i++) {
            trapezResult = integrateTrapezSiceAlex(lower, upper, fun, params, x, trapezResult, i);
//            System.out.println("trapezResult = " + i + ", " + trapezResult);
            s = (4.0 * trapezResult - ost) / 3.0;
            if (i > 5) {
                if (Math.abs(s - os) < 1.E-3 * Math.abs(os) || (s == 0.0 && os == 0.0)) {
                    break;
                }
            }
            os = s;
            ost = trapezResult;
        }

        return trapezResult;
    }

    private static double integrateTrapezSiceAlex(double lower, double upper,
                                                  ParametricUnivariateFunction fun,
                                                  double[] params,
                                                  double[] x,
                                                  double tmpSum,
                                                  int n) {

        final double dx = x[1] - x[0];
        int lowerIndex = (int) ((lower - x[0]) / dx);
        lowerIndex = Math.min(Math.max(lowerIndex, 0), x.length - 1);
        int upperIndex = (int) ((upper - x[0]) / dx);
        upperIndex = Math.max(0, Math.min(upperIndex, x.length - 1));

        double s;
        if (n == 1) {
            double funLower = fun.value(x[lowerIndex], params);
            double funUpper = fun.value(x[upperIndex], params);
            s = 0.5 * (upper - lower) * (funLower + funUpper);
        } else {
            int it = (int) Math.pow(2.0, n - 2);
            double del = (upper - lower) / it;
            double sum = 0.0;
            double xx = x[lowerIndex] + 0.5 * del;
            for (int i = 0; i < it; i++) {
                sum += fun.value(xx, params);
                xx += del;
            }
            s = 0.5 * (tmpSum + (upper - lower) * sum / it);
        }
        return s;
    }

}

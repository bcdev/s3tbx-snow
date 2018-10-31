package org.esa.s3tbx.snow.math;

/**
 * Class providing methods for numerical integration.
 *
 * @author olafd
 *
 */
public class Integrator {

    /**
     * Provides 1D integration with Trapezoidal rule with equally spaced discretization.
     * See e.g. https://en.wikipedia.org/wiki/Trapezoidal_rule
     *
     * @param lower - lower grid boundary
     * @param upper - upper grid boundary
     * @param y - function values
     * @param x - grid values
     *
     * @return the integral
     */
    public static double integrateTrapezoid(double lower, double upper, double[] y, double[] x) {
        final double dx = x[1] - x[0];
        int lowerIndex = (int) ((lower - x[0]) / dx);
        lowerIndex = Math.min(Math.max(lowerIndex, 0), x.length-1);
        int upperIndex = (int) ((upper - x[0]) / dx);
        upperIndex = Math.max(0, Math.min(upperIndex, x.length-1));
        final int n = upperIndex - lowerIndex;

        double sum = 0.0;
        for (int i = lowerIndex; i < upperIndex; i++) {
            sum += 2.0 * y[i];
        }
        return (upper - lower) * (y[lowerIndex] + sum + y[upperIndex]) / (2*n);
    }

    /**
     * Provides 1D integration with Composite Simpson's rule with equally spaced discretization.
     * See e.g. https://en.wikipedia.org/wiki/Simpson%27s_rule
     *
     * @param lower - lower grid boundary
     * @param upper - upper grid boundary
     * @param y - function values
     * @param x - grid values
     *
     * @return the integral
     */
    public static double integrateSimpson(double lower, double upper, double[] y, double[] x) {
        final double dx = x[1] - x[0];
        int lowerIndex = (int) ((lower - x[0]) / dx);
        lowerIndex = Math.min(Math.max(lowerIndex, 0), x.length-1);
        int upperIndex = (int) ((upper - x[0]) / dx);
        upperIndex = Math.max(0, Math.min(upperIndex, x.length-1));
        final int n = upperIndex - lowerIndex;

        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int i = lowerIndex; i < upperIndex; i+=2) {
            sum1 += 4.0 * y[i];
        }
        for (int i = lowerIndex+1; i < upperIndex-1; i+=2) {
            sum2 += 2.0 * y[i];
        }
        return (upper - lower) * (y[lowerIndex] + sum1 + sum2 + y[upperIndex]) / (3*n);
    }

}

package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;

/**
 * Implementation of 'fun2' function as in breadboard sice.f90:
 *
 * @author olafd
 */
public class SiceFun2Function implements ParametricUnivariateFunction {

    public double[] gradient(double v, double... parms) {
        // not needed
        return null;
    }

    public double value(double x, double... parameters) throws NoDataException {
        return evaluate(x);
    }

    private static double evaluate(double x) throws NullArgumentException, NoDataException {
        final double p0 = 32.38;
        final double p1 = -160140.33;
        final double p2 = 7959.53;
        final double t1 = 85.34 * 1.e-3;
        final double t2 = 401.79 * 1.e-3;
        
        return p0 + p1 * Math.exp(-Math.max(x, 0.4) / t1) + p2 * Math.exp(-Math.max(x, 0.4) / t2);
    }
}

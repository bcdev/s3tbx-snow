package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;

/**
 * Implementation of Sigmoidal function. Provides 2 or 4 parameters as:
 * f(x) := 1/(1+exp(A*x + B))
 * f(x) := A + B/(1+exp(C*x + D))
 *
 * @author olafd
 */
public class SigmoidalFunction implements ParametricUnivariateFunction {

    private int numParms;

    public SigmoidalFunction(int numParms) {
        this.numParms = numParms;
    }

    public double[] gradient(double x, double... parms) {
        if (numParms == 4) {
            // 4 parameter approach:
            final double aDeriv = 1.0;
            final double expTerm = Math.exp(parms[2] * x + parms[3]);
            final double bDeriv = 1.0 / (1.0 + expTerm);
            final double cDeriv = -parms[1] * x * expTerm / ((1 + expTerm) * (1 + expTerm));
            final double dDeriv = -parms[1] * expTerm / ((1 + expTerm) * (1 + expTerm));

            return new double[]{aDeriv, bDeriv, cDeriv, dDeriv};
        } else if (numParms == 2) {
            // 2 parameter approach:
            final double expTerm = Math.exp(parms[0] * x + parms[1]);
            final double cDeriv = -x * expTerm / ((1 + expTerm) * (1 + expTerm));
            final double dDeriv = -expTerm / ((1 + expTerm) * (1 + expTerm));

            return new double[]{cDeriv, dDeriv};
        } else {
            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
        }
    }

    public double value(double x, double... parameters) throws NoDataException {
        return SigmoidalFunction.evaluate(parameters, x);
    }

    private static double evaluate(double[] coeffs, double x) throws NullArgumentException, NoDataException {
        MathUtils.checkNotNull(coeffs);
        final int numParms = coeffs.length;

        if (numParms == 4) {
            // approach with 4 parameters: f(x) := A + B/(1+exp(C*x + D))
            return coeffs[0] + coeffs[1] / (1.0 + Math.exp(coeffs[2]*x + coeffs[3]));
        } else if (numParms == 2) {
            // approach with just 2 parameters: f(x) := 1/(1+exp(A*x + B))
            return 1.0 / (1.0 + Math.exp(coeffs[0] * x + coeffs[1]));
        } else {
            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
        }
    }

}

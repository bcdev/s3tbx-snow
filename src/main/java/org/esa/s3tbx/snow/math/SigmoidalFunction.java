package org.esa.s3tbx.snow.math;

import org.apache.axiom.util.UIDGenerator;
import org.apache.commons.math3.analysis.DifferentiableUnivariateFunction;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;

import java.io.Serializable;
import java.util.UUID;

/**
 * Created by Olaf on 26.09.2017.
 */
public class SigmoidalFunction implements UnivariateDifferentiableFunction, DifferentiableUnivariateFunction, Serializable {

    private final double[] coefficients;

    public SigmoidalFunction(double[] c) throws NullArgumentException, NoDataException {
        MathUtils.checkNotNull(c);
        int n = c.length;
        if(n == 0) {
            throw new NoDataException(LocalizedFormats.ARRAY_ZERO_LENGTH_OR_NULL_NOT_ALLOWED);
        } else {
            while(n > 1 && c[n - 1] == 0.0D) {
                --n;
            }

            this.coefficients = new double[n];
            System.arraycopy(c, 0, this.coefficients, 0, n);
        }
    }

    protected static double[] differentiate(double[] coefficients) throws NullArgumentException, NoDataException {

        // todo: taken from PolynominalFunction - adapt
        MathUtils.checkNotNull(coefficients);
        int n = coefficients.length;
        if(n == 0) {
            throw new NoDataException(LocalizedFormats.ARRAY_ZERO_LENGTH_OR_NULL_NOT_ALLOWED);
        } else if(n == 1) {
            return new double[]{0.0D};
        } else {
            double[] result = new double[n - 1];

            for(int i = n - 1; i > 0; --i) {
                result[i - 1] = (double)i * coefficients[i];
            }

            return result;
        }
    }


    public SigmoidalFunction sigmoidalDerivative() {
        return new SigmoidalFunction(differentiate(this.coefficients));
    }

    public UnivariateFunction derivative() {
        return this.sigmoidalDerivative();
    }

    public DerivativeStructure value(DerivativeStructure derivativeStructure) throws DimensionMismatchException {
        return null;
    }

    public double value(double x) {
        return evaluate(this.coefficients, x);
    }

    protected static double evaluate(double[] coefficients, double argument) throws NullArgumentException, NoDataException {
        MathUtils.checkNotNull(coefficients);
        int n = coefficients.length;
        if(n != 4) {
            throw new NoDataException(LocalizedFormats.DIMENSIONS_MISMATCH);
        } else {
            // AK's fitting function:
            return coefficients[0] + coefficients[1] / (1.0 + Math.exp((argument-coefficients[2])/coefficients[3]));
        }
    }


    public static class Parametric implements ParametricUnivariateFunction {
        public Parametric() {
        }

        public double[] gradient(double x, double... parms) {
            // todo: check if this is correct
            final double aDeriv = 1.0;
            final double EXP = Math.exp((x - parms[2]) / parms[3]);
            final double bDeriv = 1.0 / (1.0 + EXP);
            final double cDeriv = EXP/(parms[3]*(1+EXP)*(1+EXP));
            final double dDeriv = (parms[3] - x)*EXP/(parms[3]*(1+EXP)*(1+EXP));

            return new double[]{aDeriv, bDeriv, cDeriv, dDeriv};
        }

        public double value(double x, double... parameters) throws NoDataException {
            return SigmoidalFunction.evaluate(parameters, x);
        }
    }
}

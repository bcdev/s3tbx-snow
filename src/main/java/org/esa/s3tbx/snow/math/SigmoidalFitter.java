package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.fitting.CurveFitter;
import org.apache.commons.math3.optim.nonlinear.vector.MultivariateVectorOptimizer;

/**
 * Implements Apache commons curve fitting for Sigmoidal type function
 *
 * @author olafd
 *
 */
public class SigmoidalFitter extends CurveFitter<SigmoidalFunction> {
    public SigmoidalFitter(MultivariateVectorOptimizer optimizer) {
        super(optimizer);
    }

    /**
     * Provides fitting parameters for sigmoidal fit
     *
     * @param guess - initial guess of parameter array
     * @param numParms - number of parameters, i.e. 2 or 4
     *
     * @return parameter array
     */
    public double[] fit(double[] guess, int numParms) {
        return this.fit(new SigmoidalFunction(numParms), guess);
    }
}
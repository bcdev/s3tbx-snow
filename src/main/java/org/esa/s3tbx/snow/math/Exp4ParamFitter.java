package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.fitting.CurveFitter;
import org.apache.commons.math3.optim.nonlinear.vector.MultivariateVectorOptimizer;

/**
 * Implements Apache commons curve fitting for an exponential type function with 4 parameters
 *
 * @author olafd
 *
 */
//public class Exp4ParamFitter extends CurveFitter<Exp4ParamFunction> {
public class Exp4ParamFitter extends CurveFitter<Exp4ParamFunction2> {
    public Exp4ParamFitter(MultivariateVectorOptimizer optimizer) {
        super(optimizer);
    }

    /**
     * Provides fitting parameters for sigmoidal fit
     *
     * @param guess - initial guess of parameter array
     *
     * @return parameter array
     */
    public double[] fit(double[] guess) {
//        return this.fit(new Exp4ParamFunction(), guess);
        return this.fit(new Exp4ParamFunction2(), guess);
    }
}

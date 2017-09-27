package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.fitting.CurveFitter;
import org.apache.commons.math3.optim.nonlinear.vector.MultivariateVectorOptimizer;

/**
 * Created by Olaf on 26.09.2017.
 */
public class SigmoidalFitter extends CurveFitter<SigmoidalFunction.Parametric> {
    public SigmoidalFitter(MultivariateVectorOptimizer optimizer) {
        super(optimizer);
    }

    public double[] fit(double[] guess) {
        return this.fit(new SigmoidalFunction.Parametric(), guess);
    }
}

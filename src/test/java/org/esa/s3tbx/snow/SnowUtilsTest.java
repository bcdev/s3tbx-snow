package org.esa.s3tbx.snow;

import org.apache.commons.math3.fitting.PolynomialFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.esa.s3tbx.snow.math.Exp4ParamFitter;
import org.esa.s3tbx.snow.math.Exp4ParamFunction;
import org.esa.s3tbx.snow.math.SigmoidalFitter;
import org.junit.Ignore;
import org.junit.Test;

import static junit.framework.Assert.assertNotNull;
import static junit.framework.TestCase.assertEquals;


public class SnowUtilsTest {

    @Test
    @Ignore
    public void testComputeKappa2() throws Exception {

        final double[] kappa2 = SnowUtils.computeKappa2();
        assertNotNull(kappa2);
        assertEquals(21, kappa2.length);
        assertEquals(2.365E-11, kappa2[0], 1.E-12);

    }

}

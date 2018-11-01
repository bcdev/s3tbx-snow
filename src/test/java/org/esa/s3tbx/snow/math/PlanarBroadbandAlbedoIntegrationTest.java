package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class PlanarBroadbandAlbedoIntegrationTest {

    @Test
    public void testIntegrateTrapezoidal() {
        double[] x = new double[1001];
        double[] y = new double[1001];
        for (int i = 0; i < 1001; i++) {
            x[i] = 0.001*i;
            y[i] = x[i]*x[i];
        }
        double result = Integrator.integrateTrapezoid(0.0, 1.0, y, x);
        assertEquals(0.333333, result, 1.E-3);

        for (int i = 0; i < 1001; i++) {
            x[i] = 0.001*i*Math.PI;
            y[i] = Math.cos(x[i]);
        }
        result = Integrator.integrateTrapezoid(0.0, Math.PI, y, x);
        assertEquals(0.0, result, 1.E-2);
    }

    @Test
    public void testIntegrateSimpson() {
        double[] x = new double[1001];
        double[] y = new double[1001];
        for (int i = 0; i < 1001; i++) {
            x[i] = 0.001*i;
            y[i] = x[i]*x[i];
        }
        double result = Integrator.integrateSimpson(0.0, 1.0, y, x);
        assertEquals(0.333333, result, 1.E-3);

        for (int i = 0; i < 1001; i++) {
            x[i] = 0.001*i*Math.PI;
            y[i] = Math.cos(x[i]);
        }
        result = Integrator.integrateSimpson(0.0, Math.PI, y, x);
        assertEquals(0.0, result, 1.E-2);
    }

    @Test
    public void testXSquareIntegrationWithApacheSimpson() {
        // just a simple integration test...
        XSquareFunction xSquareFunction = new XSquareFunction();
        SimpsonIntegrator integrator = new SimpsonIntegrator();
        assertEquals(0.333333, integrator.integrate(SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT,
                xSquareFunction, 0.0, 1.0), 1.E-3);
    }

    private class XSquareFunction implements UnivariateFunction {

        public double value(double x) {
            return x * x;
        }
    }
}

package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class PlanarBroadbandAlbedoIntegrationTest {

    @Test
    public void testPlanarBroadbandAlbedoIntegration() throws Exception {
        final double mu_0 = 0.251;
        final double p = 1013.25;
        SimpsonIntegrator integrator = new SimpsonIntegrator();

        PlanarBroadbandAlbedoIntegrand planarBroadbandAlbedoNumerator =
                new PlanarBroadbandAlbedoIntegrand(mu_0, p, PlanarBroadbandAlbedoIntegrand.NUMERATOR);
        final double numeratorResult = integrator.integrate(256, planarBroadbandAlbedoNumerator, 0.4, 1.02);

        PlanarBroadbandAlbedoIntegrand planarBroadbandAlbedoDenominator =
                new PlanarBroadbandAlbedoIntegrand(mu_0, p, PlanarBroadbandAlbedoIntegrand.DENOMINATOR);
        final double denominatorResult = integrator.integrate(1024, planarBroadbandAlbedoDenominator, 0.4, 4.0);

        System.out.println("denominatorResult = " + denominatorResult);
        System.out.println("numeratorResult = " + numeratorResult);
        System.out.println("total Result = " + numeratorResult/denominatorResult);
    }

    @Test
    public void testXSquareIntegration() throws Exception {
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

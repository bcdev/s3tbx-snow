package org.esa.s3tbx.snow;

import org.apache.commons.math3.fitting.*;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.esa.s3tbx.snow.math.SigmoidalFitter;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertTrue;


public class OlciSnowAlbedoAlgorithmTest {

    @Test
    public void testComputeBroadbandAlbedos() throws Exception {
        double[] spectralAlbedos = new double[]{
//                0.998,
                0.998, 0.998, 0.996, 0.993, 0.99, 0.984, 0.975, 0.964, 0.961,
                0.95, 0.92, 0.89, 0.86, 0.712
        };
        double r21 = 0.71233;     // theoretical brr_21 for snow of 200mu diameter
        double sza = 50.0;
        final OlciSnowAlbedoAlgorithm.SphericalBroadbandAlbedo sbbaTerms =
                OlciSnowAlbedoAlgorithm.computeSphericalBroadbandAlbedoTerms(spectralAlbedos, r21);
        final double sphericalBroadbandAlbedo = sbbaTerms.getR_b1() + sbbaTerms.getR_b2();
        assertEquals(0.8385, sphericalBroadbandAlbedo, 1.E-2);
        final double broadbandPlanarAlbedo =
                OlciSnowAlbedoAlgorithm.computePlanarFromSphericalAlbedo(sphericalBroadbandAlbedo, sza);
        assertEquals(0.8416, broadbandPlanarAlbedo, 1.E-2);
    }

    @Test
    public void testComputeGrainDiameter() throws Exception {
        final double grainDiameter = OlciSnowAlbedoAlgorithm.computeGrainDiameter(0.71233);
        assertEquals(316.765, grainDiameter, 1.E-2);
    }

    @Test
    public void testComputeSpectralAlbedo() throws Exception {
        double brr = 0.71233;
        double sza = 55.5;
        double vza = 31.55;
        double saa = 154.28;
        double vaa = 103.75;
        final double spectralAlbedo = OlciSnowAlbedoAlgorithm.computeSpectralAlbedo(brr, sza, vza, saa, vaa);
        assertEquals(0.7244, spectralAlbedo, 1.E-2);
    }

    @Test
    public void testIntegrateR_b1() throws Exception {
        double[] spectralAlbedos = new double[]{
//                0.998,
                0.998, 0.998, 0.996, 0.993, 0.99, 0.984, 0.975, 0.964, 0.961,
                0.95, 0.92, 0.89, 0.86, 0.712
        };
        double r_b1 = OlciSnowAlbedoAlgorithm.integrateR_b1(spectralAlbedos);
        assertEquals(0.7552, r_b1, 1.E-2);
    }

    @Test
    public void testIntegrateR_b2() throws Exception {

        double r_b2 = OlciSnowAlbedoAlgorithm.integrateR_b2(200.0);
        assertEquals(0.0947, r_b2, 1.E-2);

        r_b2 = OlciSnowAlbedoAlgorithm.integrateR_b2(50.0);
        assertEquals(0.1289, r_b2, 1.E-2);

        r_b2 = OlciSnowAlbedoAlgorithm.integrateR_b2(100.0);
        assertEquals(0.1118, r_b2, 1.E-2);

        r_b2 = OlciSnowAlbedoAlgorithm.integrateR_b2(400.0);
        assertEquals(0.077, r_b2, 1.E-2);

        r_b2 = OlciSnowAlbedoAlgorithm.integrateR_b2(800.0);
        assertEquals(0.0604, r_b2, 1.E-2);
    }

    @Test
    public void testInterpolateSpectralAlbedos() {
        //preparation
        double[] x = {0, 50, 100};
        double[] y = {0, 50, 200};
        double[] xi = {10, 25, 75, 80, 100};

        double[] yi = OlciSnowAlbedoAlgorithm.interpolateSpectralAlbedos(x, y, xi);

        //assertion
        assertEquals(5, yi.length, 1.E-2);
        assertEquals(10.0, yi[0]);
        assertEquals(25.0, yi[1]);
        assertEquals(125.0, yi[2]);
        assertEquals(140.0, yi[3]);
        assertEquals(200.0, yi[4]);
    }

    @Test
    @Ignore
    public void testPolynominalCurveFitting() throws Exception {
        // this works only for apache commons-math3-3.6.1 !!

        // fit to third-order polynom
//        PolynomialCurveFitter curveFitter = PolynomialCurveFitter.create(3);
//        final WeightedObservedPoints obs = new WeightedObservedPoints();
//
//        obs.add(0.4, 1.0);
//        obs.add(0.753, 0.963);
//        obs.add(0.865, 0.922);
//        obs.add(1.02, 0.737);
//        final double[] curveFit = curveFitter.fit(obs.toList());
//
//        for (int i = 0; i < curveFit.length; i++) {
//            System.out.printf("%d,%s%n", i, curveFit[i]);
//        }
        assertTrue(true);
    }

    @Test
    public void testPolynominal3rdOrderCurveFitting() throws Exception {
//        final double[] initialGuess = {0., 0., 0., 0.};
        final double[] initialGuess = {0., 0., 0., 0.,0., 0., 0., 0.};

        PolynomialFitter curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());

        curveFitter.addObservedPoint(0.4, 1.0);
        curveFitter.addObservedPoint(0.753, 0.963);
        curveFitter.addObservedPoint(0.865, 0.922);
        curveFitter.addObservedPoint(1.02, 0.737);
        final double[] fit = curveFitter.fit(initialGuess);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("%d,%s%n", i, fit[i]);
        }

    }

    @Test
    public void testPolynominal7thOrderCurveFitting() throws Exception {
        final double[] initialGuess = {0., 0., 0., 0.};
//        final double[] initialGuess = {0., 0., 0., 0.,0., 0., 0., 0.};

        PolynomialFitter curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());

        curveFitter.addObservedPoint(400., 1.0);
        curveFitter.addObservedPoint(753., 0.963);
        curveFitter.addObservedPoint(865., 0.922);
        curveFitter.addObservedPoint(1020., 0.737);
        final double[] fit = curveFitter.fit(initialGuess);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("3rd order: %d,%s%n", i, fit[i]);
        }

    }

    @Test
    public void testSigmoidalCurveFitting() throws Exception {
//        final double[] initialGuess = {1., 1., 1., 1.};
        final double[] initialGuess = {1., 1.};

        SigmoidalFitter curveFitter = new SigmoidalFitter(new LevenbergMarquardtOptimizer());

        curveFitter.addObservedPoint(0.4, 1.0);
        curveFitter.addObservedPoint(0.753, 0.963);
        curveFitter.addObservedPoint(0.865, 0.922);
        curveFitter.addObservedPoint(1.02, 0.737);
        final double[] fit = curveFitter.fit(initialGuess);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("Sigmoidal: %d,%s%n", i, fit[i]);
        }

    }



}

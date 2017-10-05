package org.esa.s3tbx.snow;

import org.apache.commons.math3.fitting.PolynomialFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.esa.s3tbx.snow.math.Exp4ParamFitter;
import org.esa.s3tbx.snow.math.Exp4ParamFunction;
import org.esa.s3tbx.snow.math.SigmoidalFitter;
import org.junit.Ignore;
import org.junit.Test;

import static junit.framework.TestCase.assertEquals;


public class OlciSnowAlbedoAlgorithmTest {

    @Test
    @Ignore
    public void testComputeBroadbandAlbedos() throws Exception {
        double[] spectralAlbedos = new double[]{
//                0.998,
                0.998, 0.998, 0.996, 0.993, 0.99, 0.984, 0.975, 0.964, 0.961,
                0.95, 0.92, 0.89, 0.86, 0.712
        };
        double sza = 50.0;
        final OlciSnowAlbedoAlgorithm.SphericalBroadbandAlbedo sbbaTerms =
                OlciSnowAlbedoAlgorithm.computeSphericalBroadbandAlbedoTerms(spectralAlbedos);
        final double sphericalBroadbandAlbedo = sbbaTerms.getR_b1() + sbbaTerms.getR_b2();
        assertEquals(0.8385, sphericalBroadbandAlbedo, 1.E-2);
        final double broadbandPlanarAlbedo =
                OlciSnowAlbedoAlgorithm.computePlanarFromSphericalAlbedo(sphericalBroadbandAlbedo, sza);
        assertEquals(0.8416, broadbandPlanarAlbedo, 1.E-2);
    }

    @Test
    public void testComputeGrainDiameter() throws Exception {
        final double grainDiameter = OlciSnowAlbedoAlgorithm.computeGrainDiameter(0.7372415128980274);
        assertEquals(260.18598, grainDiameter, 1.E-2);
    }

//    @Test
//    public void testComputeSpectralAlbedo() throws Exception {
//        double brr = 0.71233;
//        double sza = 55.5;
//        double vza = 31.55;
//        double saa = 154.28;
//        double vaa = 103.75;
//        final double spectralAlbedo = OlciSnowAlbedoAlgorithm.computeSpectralAlbedo_old(brr, sza, vza, saa, vaa);
//        assertEquals(0.7244, spectralAlbedo, 1.E-2);
//    }

//    @Test
//    public void testIntegrateR_b1() throws Exception {
//        double[] spectralAlbedos = new double[]{
////                0.998,
//                0.998, 0.998, 0.996, 0.993, 0.99, 0.984, 0.975, 0.964, 0.961,
//                0.95, 0.92, 0.89, 0.86, 0.712
//        };
//        double r_b1 = OlciSnowAlbedoAlgorithm.integrateR_b1(spectralAlbedos);
//        assertEquals(0.7552, r_b1, 1.E-2);
//    }

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

    // this works only for apache commons-math3-3.6.1 !!
//    @Test
//    public void testPolynominalCurveFitting() throws Exception {
//
//
//        // fit to third-order polynom
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
//        assertTrue(true);
//    }

    @Test
    public void testPolynominalCurveFitting() throws Exception {
        double[] initialGuess = {0., 0., 0., 0.};

        PolynomialFitter curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());

        curveFitter.addObservedPoint(0.4, 1.0);
        curveFitter.addObservedPoint(0.753, 0.963);
        curveFitter.addObservedPoint(0.865, 0.922);
        curveFitter.addObservedPoint(1.02, 0.737);
        double[] fit = curveFitter.fit(initialGuess);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("3rd order fit: %d,%s%n", i, fit[i]);
        }

        initialGuess = new double[]{0., 0., 0., 0., 0., 0., 0., 0.};
        fit = curveFitter.fit(initialGuess);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("7th order fit: %d,%s%n", i, fit[i]);
        }
    }

    @Test
    public void testSigmoidalCurveFitting() throws Exception {

        SigmoidalFitter curveFitter = new SigmoidalFitter(new LevenbergMarquardtOptimizer());

        curveFitter.addObservedPoint(0.4, 1.0);
        curveFitter.addObservedPoint(0.753, 0.963);
        curveFitter.addObservedPoint(0.865, 0.922);
        curveFitter.addObservedPoint(1.02, 0.737);
        double[] initialGuess = {1., 1.};
        double[] fit = curveFitter.fit(initialGuess, 2);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("Sigmoidal 2 parameter fit: %d,%s%n", i, fit[i]);
        }

        initialGuess = new double[]{1., 1., 1., 1.};
        fit = curveFitter.fit(initialGuess, 4);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("Sigmoidal 4 parameter fit: %d,%s%n", i, fit[i]);
        }
    }

    @Test
    public void testAlgoSep22() throws Exception {
        final double[] rhoToa = new double[]{
                0.8732, 0.8710, 0.8738, 0.8583, 0.8255, 0.7338, 0.7208,
                0.8050, 0.8184, 0.8246, 0.8273, 0.8451, 0.1507, 0.3121,
                0.7134, 0.8263, 0.8132, 0.7909, 0.6403, 0.3507, 0.6673
        };
        final double sza = 75.4347;
        final double vza = 26.35964;
        final double saa = -52.43322;
        final double vaa = -118.132866;

        final double[] spectralSphericalAlbedos =
                OlciSnowAlbedoAlgorithm.computeSpectralSphericalAlbedos(rhoToa, sza, vza, SpectralAlbedoMode.SIGMOIDAL);

        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            System.out.printf("Spectral albedos: %f,%s%n", wvl, spectralSphericalAlbedos[i]);
        }
    }

    @Test
    public void testExp4ParamCurveFitting() throws Exception {

        Exp4ParamFitter curveFitter = new Exp4ParamFitter(new LevenbergMarquardtOptimizer());

        curveFitter.addObservedPoint(0.4, 1.0);
        curveFitter.addObservedPoint(0.49, 0.97365);
        curveFitter.addObservedPoint(0.865, 0.843676);
        curveFitter.addObservedPoint(1.02, 0.689266);

        double[] initialGuess = new double[]{1., 1., 1., 1.};
        double[] fit = curveFitter.fit(initialGuess);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("Exp4Param 4 parameter fit: %d,%s%n", i, fit[i]);
        }

        final double[] spectralSphericalAlbedos = new double[OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length];
        final double[] L = new double[OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length];

        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            L[i] = 2.0*OlciSnowAlbedoConstants.KAPPA_2[i]*Math.PI/fit[2];
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            spectralSphericalAlbedos[i] = Math.exp(-(fit[1] + fit[2]/wvl));
            System.out.printf("Exp4Param albedos: %f,%s%n", wvl, spectralSphericalAlbedos[i]);
        }
        System.out.println();
    }

}

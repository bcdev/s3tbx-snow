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
        double[] initialGuess = {0., 0., 0., 0., 0.};

        PolynomialFitter curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());

        // 350-525nm Im(K), https://atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat:
        curveFitter.addObservedPoint(4.000E-001,  2.365E-011);
        curveFitter.addObservedPoint(4.100E-001,  2.669E-011);
        curveFitter.addObservedPoint(4.200E-001,  3.135E-011);
        curveFitter.addObservedPoint(4.300E-001,  4.140E-011);
        curveFitter.addObservedPoint(4.400E-001,  6.268E-011);
        curveFitter.addObservedPoint(4.500E-001,  9.239E-011);
        curveFitter.addObservedPoint(4.600E-001,  1.325E-010);
        curveFitter.addObservedPoint(4.700E-001,  1.956E-010);
        curveFitter.addObservedPoint(4.800E-001,  2.861E-010);
        curveFitter.addObservedPoint(4.900E-001,  4.172E-010);
        curveFitter.addObservedPoint(5.000E-001,  5.889E-010);
        curveFitter.addObservedPoint(5.100E-001,  8.036E-010);
        curveFitter.addObservedPoint(5.200E-001,  1.076E-009);

        double[] fit = curveFitter.fit(initialGuess);
        for (int i = 0; i < fit.length; i++) {
            System.out.printf("4th order fit 400-525nm: %d,%s%n", i, fit[i]);
        }

        // 525-700nm Im(K), https://atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat:
        curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());
        curveFitter.addObservedPoint(5.200E-001,  1.076E-009);
        curveFitter.addObservedPoint(5.300E-001,  1.409E-009);
        curveFitter.addObservedPoint(5.400E-001,  1.813E-009);
        curveFitter.addObservedPoint(5.500E-001,  2.289E-009);
        curveFitter.addObservedPoint(5.600E-001,  2.839E-009);
        curveFitter.addObservedPoint(5.700E-001,  3.461E-009);
        curveFitter.addObservedPoint(5.800E-001,  4.159E-009);
        curveFitter.addObservedPoint(5.900E-001,  4.930E-009);
        curveFitter.addObservedPoint(6.000E-001,  5.730E-009);
        curveFitter.addObservedPoint(6.100E-001,  6.890E-009);
        curveFitter.addObservedPoint(6.200E-001,  8.580E-009);
        curveFitter.addObservedPoint(6.300E-001,  1.040E-008);
        curveFitter.addObservedPoint(6.400E-001,  1.220E-008);
        curveFitter.addObservedPoint(6.500E-001,  1.430E-008);
        curveFitter.addObservedPoint(6.600E-001,  1.660E-008);
        curveFitter.addObservedPoint(6.700E-001,  1.890E-008);
        curveFitter.addObservedPoint(6.800E-001,  2.090E-008);
        curveFitter.addObservedPoint(6.900E-001,  2.400E-008);
        curveFitter.addObservedPoint(7.000E-001,  2.900E-008);

        fit = curveFitter.fit(initialGuess);
        for (int i = 0; i < fit.length; i++) {
            System.out.printf("4th order fit 525-700nm: %d,%s%n", i, fit[i]);
        }

        // 700-1980nm Im(K), https://atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat:
        curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());
        curveFitter.addObservedPoint(7.000E-001, 2.900E-008);
        curveFitter.addObservedPoint(7.100E-001, 3.440E-008);
        curveFitter.addObservedPoint(7.200E-001, 4.030E-008);
        curveFitter.addObservedPoint(7.300E-001, 4.300E-008);
        curveFitter.addObservedPoint(7.400E-001, 4.920E-008);
        curveFitter.addObservedPoint(7.500E-001, 5.870E-008);
        curveFitter.addObservedPoint(7.600E-001, 7.080E-008);
        curveFitter.addObservedPoint(7.700E-001, 8.580E-008);
        curveFitter.addObservedPoint(7.800E-001, 1.020E-007);
        curveFitter.addObservedPoint(7.900E-001, 1.180E-007);
        curveFitter.addObservedPoint(8.000E-001, 1.340E-007);
        curveFitter.addObservedPoint(8.100E-001, 1.400E-007);
        curveFitter.addObservedPoint(8.200E-001, 1.430E-007);
        curveFitter.addObservedPoint(8.300E-001, 1.450E-007);
        curveFitter.addObservedPoint(8.400E-001, 1.510E-007);
        curveFitter.addObservedPoint(8.500E-001, 1.830E-007);
        curveFitter.addObservedPoint(8.600E-001, 2.150E-007);
        curveFitter.addObservedPoint(8.700E-001, 2.650E-007);
        curveFitter.addObservedPoint(8.800E-001, 3.350E-007);
        curveFitter.addObservedPoint(8.900E-001, 3.920E-007);
        curveFitter.addObservedPoint(9.000E-001, 4.200E-007);
        curveFitter.addObservedPoint(9.100E-001, 4.440E-007);
        curveFitter.addObservedPoint(9.200E-001, 4.740E-007);
        curveFitter.addObservedPoint(9.300E-001, 5.110E-007);
        curveFitter.addObservedPoint(9.400E-001, 5.530E-007);
        curveFitter.addObservedPoint(9.500E-001, 6.020E-007);
        curveFitter.addObservedPoint(9.600E-001, 7.550E-007);
        curveFitter.addObservedPoint(9.700E-001, 9.260E-007);
        curveFitter.addObservedPoint(9.800E-001, 1.120E-006);
        curveFitter.addObservedPoint(9.900E-001, 1.330E-006);
        curveFitter.addObservedPoint(1.000E+000, 1.620E-006);
        curveFitter.addObservedPoint(1.010E+000, 2.000E-006);
        curveFitter.addObservedPoint(1.020E+000, 2.250E-006);
        curveFitter.addObservedPoint(1.030E+000, 2.330E-006);

        initialGuess = new double[]{1., 1., 1., 1., 1.};
        fit = curveFitter.fit(initialGuess);
        for (int i = 0; i < fit.length; i++) {
            System.out.printf("4th order fit 700-1020nm: %d,%s%n", i, fit[i]);
        }
    }

    @Test
    public void testSigmoidalCurveFitting() throws Exception {

//        SigmoidalFitter curveFitter = new SigmoidalFitter(new LevenbergMarquardtOptimizer());
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
    public void testSigmoidalCurveFittingPerformance() throws Exception {

//        SigmoidalFitter curveFitter = new SigmoidalFitter(new LevenbergMarquardtOptimizer());
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
    @Ignore
    public void testAlgoSep22() throws Exception {
        final double[] rhoToa = new double[]{
                0.8732, 0.8710, 0.8738, 0.8583, 0.8255, 0.7338, 0.7208,
                0.8050, 0.8184, 0.8246, 0.8273, 0.8451, 0.1507, 0.3121,
                0.7134, 0.8263, 0.8132, 0.7909, 0.6403, 0.3507, 0.6673
        };
        final double sza = 75.4347;
        final double vza = 26.35964;

        final double[] spectralSphericalAlbedos =
                OlciSnowAlbedoAlgorithm.computeSpectralSphericalAlbedos(rhoToa, sza, vza,
                                                                        SpectralAlbedoMode.SIGMOIDAL_FIT);

        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            System.out.printf("Spectral albedos: %f,%s%n", wvl, spectralSphericalAlbedos[i]);
        }
    }

    @Test
    public void testExp4ParamCurveFitting() throws Exception {

        Exp4ParamFitter curveFitter = new Exp4ParamFitter(new LevenbergMarquardtOptimizer());

        curveFitter.addObservedPoint(0.4, 1.0);
        curveFitter.addObservedPoint(0.753, 0.768120321063813);
        curveFitter.addObservedPoint(0.865, 0.7226259625191246);
        curveFitter.addObservedPoint(1.02, 0.48345802606075666);

        // fit does not always work, result strongly depends on initial guess!
        // see AK TN, 20171005
        // f(v; a, kappa1, L, b) := a *(exp(-(2*PI/v) * (kappa1*L + kappa2(v)*L)))^b

        double[] initialGuess = new double[]{1.068312113070488, 0., 85750.9860257803, 0.7266642582580547};
        double[] fit = curveFitter.fit(initialGuess);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("Exp4Param 4 parameter fit: %d,%s%n", i, fit[i]);
        }

        final double[] refl = new double[OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length];
        final double[] spectralSphericalAlbedos = new double[OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length];
        final Exp4ParamFunction exp4ParamFunction = new Exp4ParamFunction();
        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            spectralSphericalAlbedos[i] = exp4ParamFunction.value(wvl, fit);
            System.out.println("Exp4Param albedos: " + wvl + ", " + refl[i] + ", "  + spectralSphericalAlbedos[i]);
        }

        for (int i = 0; i < 40; i++) {
            final double wvl = 0.4 + 0.02*i;
            System.out.println("Exp4Param fit: " + wvl +  ", "  + exp4ParamFunction.value(wvl, fit));
        }

        System.out.println();
    }

    @Test
    public void testOlciGains() throws Exception {
//        final double sza = 75.33;
//        final double vza = 25.495424;
//        final double[] brrGains =
//        final double[] sphericalAlbedosGains =
//                OlciSnowAlbedoAlgorithm.computeSpectralSphericalAlbedos(brrGains, sza, vza, SpectralAlbedoMode.SIGMOIDAL_FIT);
    }
}

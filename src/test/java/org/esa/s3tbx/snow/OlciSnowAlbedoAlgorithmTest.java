package org.esa.s3tbx.snow;

import org.apache.commons.math3.fitting.PolynomialFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.esa.s3tbx.snow.math.Exp4ParamFitter;
import org.esa.s3tbx.snow.math.Exp4ParamFunction;
import org.esa.s3tbx.snow.math.SigmoidalFitter;
import org.esa.snap.core.util.math.MathUtils;
import org.junit.Test;

import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertNotNull;


public class OlciSnowAlbedoAlgorithmTest {

    @Test
    public void testComputeGrainDiameter() throws Exception {
        final double grainDiameter = OlciSnowAlbedoAlgorithm.computeGrainDiameter(0.7372415128980274, 1020.0);
        assertEquals(255.819, grainDiameter, 1.E-2);
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
    public void testPolynominalCurveFitting() throws Exception {
        double[] initialGuess = {0., 0., 0., 0., 0.};

        PolynomialFitter curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());

        // 350-525nm Im(K), https://atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat:
        curveFitter.addObservedPoint(3.000E-001, 2.0E-011);
        curveFitter.addObservedPoint(3.500E-001, 2.0E-011);
        curveFitter.addObservedPoint(3.900E-001, 2.0E-011);
        curveFitter.addObservedPoint(4.000E-001, 2.365E-011);
        curveFitter.addObservedPoint(4.100E-001, 2.669E-011);
        curveFitter.addObservedPoint(4.200E-001, 3.135E-011);
        curveFitter.addObservedPoint(4.300E-001, 4.140E-011);
        curveFitter.addObservedPoint(4.400E-001, 6.268E-011);
        curveFitter.addObservedPoint(4.500E-001, 9.239E-011);
        curveFitter.addObservedPoint(4.600E-001, 1.325E-010);
        curveFitter.addObservedPoint(4.700E-001, 1.956E-010);
        curveFitter.addObservedPoint(4.800E-001, 2.861E-010);
        curveFitter.addObservedPoint(4.900E-001, 4.172E-010);
        curveFitter.addObservedPoint(5.000E-001, 5.889E-010);
        curveFitter.addObservedPoint(5.100E-001, 8.036E-010);
        curveFitter.addObservedPoint(5.200E-001, 1.076E-009);

        double[] fit = curveFitter.fit(initialGuess);
        for (int i = 0; i < fit.length; i++) {
            System.out.printf("4th order fit 400-525nm: %d,%s%n", i, fit[i]);
        }

        // 525-700nm Im(K), https://atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat:
        curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());
        curveFitter.addObservedPoint(5.200E-001, 1.076E-009);
        curveFitter.addObservedPoint(5.300E-001, 1.409E-009);
        curveFitter.addObservedPoint(5.400E-001, 1.813E-009);
        curveFitter.addObservedPoint(5.500E-001, 2.289E-009);
        curveFitter.addObservedPoint(5.600E-001, 2.839E-009);
        curveFitter.addObservedPoint(5.700E-001, 3.461E-009);
        curveFitter.addObservedPoint(5.800E-001, 4.159E-009);
        curveFitter.addObservedPoint(5.900E-001, 4.930E-009);
        curveFitter.addObservedPoint(6.000E-001, 5.730E-009);
        curveFitter.addObservedPoint(6.100E-001, 6.890E-009);
        curveFitter.addObservedPoint(6.200E-001, 8.580E-009);
        curveFitter.addObservedPoint(6.300E-001, 1.040E-008);
        curveFitter.addObservedPoint(6.400E-001, 1.220E-008);
        curveFitter.addObservedPoint(6.500E-001, 1.430E-008);
        curveFitter.addObservedPoint(6.600E-001, 1.660E-008);
        curveFitter.addObservedPoint(6.700E-001, 1.890E-008);
        curveFitter.addObservedPoint(6.800E-001, 2.090E-008);
        curveFitter.addObservedPoint(6.900E-001, 2.400E-008);
        curveFitter.addObservedPoint(7.000E-001, 2.900E-008);

        fit = curveFitter.fit(initialGuess);
        for (int i = 0; i < fit.length; i++) {
            System.out.printf("4th order fit 525-700nm: %d,%s%n", i, fit[i]);
        }

        // 700-1020nm Im(K), https://atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat:
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

        // 1020-3000nm Im(K), https://atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat:
        curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());
        curveFitter.addObservedPoint(1.030E+000, 2.330E-006);
        curveFitter.addObservedPoint(1.040E+000, 2.330E-006);
        curveFitter.addObservedPoint(1.050E+000, 2.170E-006);
        curveFitter.addObservedPoint(1.060E+000, 1.960E-006);
        curveFitter.addObservedPoint(1.070E+000, 1.810E-006);
        curveFitter.addObservedPoint(1.080E+000, 1.740E-006);
        curveFitter.addObservedPoint(1.090E+000, 1.730E-006);
        curveFitter.addObservedPoint(1.100E+000, 1.700E-006);
        curveFitter.addObservedPoint(1.110E+000, 1.760E-006);
        curveFitter.addObservedPoint(1.120E+000, 1.820E-006);
        curveFitter.addObservedPoint(1.130E+000, 2.040E-006);
        curveFitter.addObservedPoint(1.140E+000, 2.250E-006);
        curveFitter.addObservedPoint(1.150E+000, 2.290E-006);
        curveFitter.addObservedPoint(1.160E+000, 3.040E-006);
        curveFitter.addObservedPoint(1.170E+000, 3.840E-006);
        curveFitter.addObservedPoint(1.180E+000, 4.770E-006);
        curveFitter.addObservedPoint(1.190E+000, 5.760E-006);
        curveFitter.addObservedPoint(1.200E+000, 6.710E-006);
        curveFitter.addObservedPoint(1.210E+000, 8.660E-006);
        curveFitter.addObservedPoint(1.220E+000, 1.020E-005);
        curveFitter.addObservedPoint(1.230E+000, 1.130E-005);
        curveFitter.addObservedPoint(1.240E+000, 1.220E-005);
        curveFitter.addObservedPoint(1.250E+000, 1.290E-005);
        curveFitter.addObservedPoint(1.260E+000, 1.320E-005);
        curveFitter.addObservedPoint(1.270E+000, 1.350E-005);
        curveFitter.addObservedPoint(1.280E+000, 1.330E-005);
        curveFitter.addObservedPoint(1.290E+000, 1.320E-005);
        curveFitter.addObservedPoint(1.300E+000, 1.320E-005);
        curveFitter.addObservedPoint(1.310E+000, 1.310E-005);
        curveFitter.addObservedPoint(1.320E+000, 1.320E-005);
        curveFitter.addObservedPoint(1.330E+000, 1.320E-005);
        curveFitter.addObservedPoint(1.340E+000, 1.340E-005);
        curveFitter.addObservedPoint(1.350E+000, 1.390E-005);
        curveFitter.addObservedPoint(1.360E+000, 1.420E-005);
        curveFitter.addObservedPoint(1.370E+000, 1.480E-005);
        curveFitter.addObservedPoint(1.380E+000, 1.580E-005);
        curveFitter.addObservedPoint(1.390E+000, 1.740E-005);
        curveFitter.addObservedPoint(1.400E+000, 1.980E-005);
        curveFitter.addObservedPoint(1.410E+000, 3.442E-005);
        curveFitter.addObservedPoint(1.420E+000, 5.959E-005);
        curveFitter.addObservedPoint(1.430E+000, 1.028E-004);
        curveFitter.addObservedPoint(1.440E+000, 1.516E-004);
        curveFitter.addObservedPoint(1.449E+000, 2.030E-004);
        curveFitter.addObservedPoint(1.460E+000, 2.942E-004);
        curveFitter.addObservedPoint(1.471E+000, 3.987E-004);
        curveFitter.addObservedPoint(1.481E+000, 4.941E-004);
        curveFitter.addObservedPoint(1.493E+000, 5.532E-004);
        curveFitter.addObservedPoint(1.504E+000, 5.373E-004);
        curveFitter.addObservedPoint(1.515E+000, 5.143E-004);
        curveFitter.addObservedPoint(1.527E+000, 4.908E-004);
        curveFitter.addObservedPoint(1.538E+000, 4.594E-004);
        curveFitter.addObservedPoint(1.563E+000, 3.858E-004);
        curveFitter.addObservedPoint(1.587E+000, 3.105E-004);
        curveFitter.addObservedPoint(1.613E+000, 2.659E-004);
        curveFitter.addObservedPoint(1.650E+000, 2.361E-004);
        curveFitter.addObservedPoint(1.680E+000, 2.046E-004);
        curveFitter.addObservedPoint(1.700E+000, 1.875E-004);
        curveFitter.addObservedPoint(1.730E+000, 1.650E-004);
        curveFitter.addObservedPoint(1.760E+000, 1.522E-004);
        curveFitter.addObservedPoint(1.800E+000, 1.411E-004);
        curveFitter.addObservedPoint(1.830E+000, 1.302E-004);
        curveFitter.addObservedPoint(1.840E+000, 1.310E-004);
        curveFitter.addObservedPoint(1.850E+000, 1.339E-004);
        curveFitter.addObservedPoint(1.855E+000, 1.377E-004);
        curveFitter.addObservedPoint(1.860E+000, 1.432E-004);
        curveFitter.addObservedPoint(1.870E+000, 1.632E-004);
        curveFitter.addObservedPoint(1.890E+000, 2.566E-004);
        curveFitter.addObservedPoint(1.905E+000, 4.081E-004);
        curveFitter.addObservedPoint(1.923E+000, 7.060E-004);
        curveFitter.addObservedPoint(1.942E+000, 1.108E-003);
        curveFitter.addObservedPoint(1.961E+000, 1.442E-003);
        curveFitter.addObservedPoint(1.980E+000, 1.614E-003);
        curveFitter.addObservedPoint(2.000E+000, 1.640E-003);
        curveFitter.addObservedPoint(2.020E+000, 1.566E-003);
        curveFitter.addObservedPoint(2.041E+000, 1.458E-003);
        curveFitter.addObservedPoint(2.062E+000, 1.267E-003);
        curveFitter.addObservedPoint(2.083E+000, 1.023E-003);
        curveFitter.addObservedPoint(2.105E+000, 7.586E-004);
        curveFitter.addObservedPoint(2.130E+000, 5.255E-004);
        curveFitter.addObservedPoint(2.150E+000, 4.025E-004);
        curveFitter.addObservedPoint(2.170E+000, 3.235E-004);
        curveFitter.addObservedPoint(2.190E+000, 2.707E-004);
        curveFitter.addObservedPoint(2.220E+000, 2.228E-004);
        curveFitter.addObservedPoint(2.240E+000, 2.037E-004);
        curveFitter.addObservedPoint(2.245E+000, 2.026E-004);
        curveFitter.addObservedPoint(2.250E+000, 2.035E-004);
        curveFitter.addObservedPoint(2.260E+000, 2.078E-004);
        curveFitter.addObservedPoint(2.270E+000, 2.171E-004);
        curveFitter.addObservedPoint(2.290E+000, 2.538E-004);
        curveFitter.addObservedPoint(2.310E+000, 3.138E-004);
        curveFitter.addObservedPoint(2.330E+000, 3.858E-004);
        curveFitter.addObservedPoint(2.350E+000, 4.591E-004);
        curveFitter.addObservedPoint(2.370E+000, 5.187E-004);
        curveFitter.addObservedPoint(2.390E+000, 5.605E-004);
        curveFitter.addObservedPoint(2.410E+000, 5.956E-004);
        curveFitter.addObservedPoint(2.430E+000, 6.259E-004);
        curveFitter.addObservedPoint(2.460E+000, 6.820E-004);
        curveFitter.addObservedPoint(2.500E+000, 7.530E-004);
        curveFitter.addObservedPoint(2.520E+000, 7.685E-004);
        curveFitter.addObservedPoint(2.550E+000, 7.647E-004);
        curveFitter.addObservedPoint(2.565E+000, 7.473E-004);
        curveFitter.addObservedPoint(2.580E+000, 7.392E-004);
        curveFitter.addObservedPoint(2.590E+000, 7.437E-004);
        curveFitter.addObservedPoint(2.600E+000, 7.543E-004);
        curveFitter.addObservedPoint(2.620E+000, 8.059E-004);
        curveFitter.addObservedPoint(2.675E+000, 1.367E-003);
        curveFitter.addObservedPoint(2.725E+000, 3.508E-003);
        curveFitter.addObservedPoint(2.778E+000, 1.346E-002);
        curveFitter.addObservedPoint(2.817E+000, 3.245E-002);
        curveFitter.addObservedPoint(2.833E+000, 4.572E-002);
        curveFitter.addObservedPoint(2.849E+000, 6.287E-002);
        curveFitter.addObservedPoint(2.865E+000, 8.548E-002);
        curveFitter.addObservedPoint(2.882E+000, 1.198E-001);
        curveFitter.addObservedPoint(2.899E+000, 1.690E-001);
        curveFitter.addObservedPoint(2.915E+000, 2.210E-001);
        curveFitter.addObservedPoint(2.933E+000, 2.760E-001);
        curveFitter.addObservedPoint(2.950E+000, 3.120E-001);
        curveFitter.addObservedPoint(2.967E+000, 3.470E-001);
        curveFitter.addObservedPoint(2.985E+000, 3.880E-001);
        curveFitter.addObservedPoint(3.003E+000, 4.380E-001);

        initialGuess = new double[]{1., 1., 1., 1., 1.};
        fit = curveFitter.fit(initialGuess);
        for (int i = 0; i < fit.length; i++) {
            System.out.printf("4th order fit 1020-3000nm: %d,%s%n", i, fit[i]);
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

        final double[] spectralSphericalAlbedos = new double[OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length];
        final Exp4ParamFunction exp4ParamFunction = new Exp4ParamFunction();
        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            spectralSphericalAlbedos[i] = exp4ParamFunction.value(wvl, fit);
        }

        for (int i = 0; i < 32; i++) {
            final double wvl = 0.4 + 0.02 * i;
            System.out.println("Exp4Param fit: " + wvl + ", " + exp4ParamFunction.value(wvl, fit));
        }

        System.out.println();
    }

    @Test
    public void testExp4ParamCurveFitting_aoki() throws Exception {

        // this test refers to AK TN 20171010 !!!

        Exp4ParamFitter curveFitter = new Exp4ParamFitter(new LevenbergMarquardtOptimizer());

        // put AOKI values in here and just do the fit...
        curveFitter.addObservedPoint(0.4, 0.975);
        curveFitter.addObservedPoint(0.49, 0.975);
//        curveFitter.addObservedPoint(0.753, 0.93);
        curveFitter.addObservedPoint(0.865, 0.89);
        curveFitter.addObservedPoint(1.02, 0.72);

        // fit does not always work, result strongly depends on initial guess!
        // see AK TN, 20171005
        // f(v; a, kappa1, L, b) := a *(exp(-(2*PI/v) * (kappa1*L + kappa2(v)*L)))^b

        final double kappa2_1020 = 1.E-6;
        final double r_1020 = 0.72;
        final double initialGuess_2 = 1.02 * Math.log(r_1020) * Math.log(r_1020) / (2.0 * Math.PI * kappa2_1020);

        double[] initialGuess = new double[]{0.97, 0., initialGuess_2, 0.7266642582580547};
        double[] fit = curveFitter.fit(initialGuess);

        for (int i = 0; i < fit.length; i++) {
            System.out.printf("Exp4Param 4 parameter fit for AOKI: %d,%s%n", i, fit[i]);
        }

        final Exp4ParamFunction exp4ParamFunction = new Exp4ParamFunction();
        for (int i = 0; i < 32; i++) {
            final double wvl = 0.4 + 0.02 * i;
            System.out.println("Exp4Param fit for AOKI: " + wvl + ", " + exp4ParamFunction.value(wvl, fit));
        }

        System.out.println();
    }

    @Test
    public void testComputeBroadbandAlbedo() throws Exception {
        final double grainDiamMicrons = 200.0;
        final double sza = 60.0;
        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final RefractiveIndexTable refractiveIndexTable = new RefractiveIndexTable();
        refractiveIndexTable.readTableFromFile();
        final SolarSpectrumTable solarSpectrumTable = new SolarSpectrumTable();
        solarSpectrumTable.readTableFromFile();

        RefractiveIndexTable refractiveIndexInterpolatedTable =
                SnowUtils.getRefractiveIndexInterpolated(refractiveIndexTable,
                                                         solarSpectrumTable);

        double[] planarBroadbandAlbedo = null;
        for (int k=0; k<10000; k++) {
            planarBroadbandAlbedo =
                    OlciSnowAlbedoAlgorithm.computeBroadbandAlbedo(mu_0, grainDiamMicrons,
                                                                   refractiveIndexInterpolatedTable, solarSpectrumTable);
            if (k % 1000 == 0) {
                System.out.println("k1 = " + k);
            }
        }
        for (double aPlanarBroadbandAlbedo : planarBroadbandAlbedo) {
            System.out.println("planarBroadbandAlbedo = " + aPlanarBroadbandAlbedo);
        }

        double[] sphericalBroadbandAlbedo = null;
        for (int k=0; k<10000; k++) {
            sphericalBroadbandAlbedo =
                    OlciSnowAlbedoAlgorithm.computeBroadbandAlbedo(1.0, grainDiamMicrons,
                                                                   refractiveIndexInterpolatedTable, solarSpectrumTable);
            if (k % 1000 == 0) {
                System.out.println("k2 = " + k);
            }
        }
        for (double aSphericalBroadbandAlbedo : sphericalBroadbandAlbedo) {
            System.out.println("sphericalBroadbandAlbedo = " + aSphericalBroadbandAlbedo);
        }

        System.out.println();
    }

    @Test
    public void testComputePollutedSnowParms() throws Exception {

        // test input provided by AK in 'input_soot.dat' (20180124):
        final double sza = 42.3621;
        final double vza = 2.561;
        final double raa = 42.088;
        final double brr400 = 0.5816;
        final double brr1200 = 0.3704;

        final double[] pollutedSnowParams =
                OlciSnowAlbedoAlgorithm.computePollutedSnowParams(brr400, brr1200, sza, vza, raa);
        final double grainDiam = pollutedSnowParams[0];
        final double soot = pollutedSnowParams[1];
        // test results provided by AK in 'snow_properties_CHECKOUT.txt' (20180124), first line:
        assertEquals(1.537, grainDiam, 1.E-3);
        assertEquals(1.519, soot, 1.E-3);
    }

    @Test
    public void testComputeSpectralAlbedosPolluted() throws Exception {

        // test input provided by AK in 'input_soot.dat' (20180124):
        final double sza = 42.3621;
        final double grainDiam = 1.537;
        final double soot = 1.519;
        final double[] pollutedSnowParams = new double[]{grainDiam, soot};

        final OlciSnowAlbedoAlgorithm.SpectralAlbedoResult spectralAlbedosPolluted =
                OlciSnowAlbedoAlgorithm.computeSpectralAlbedosPolluted(pollutedSnowParams, sza, Double.NaN, false);
        assertNotNull(spectralAlbedosPolluted);
        final double[] sphericalAlbedos = spectralAlbedosPolluted.getSpectralAlbedos()[0];
        assertNotNull(sphericalAlbedos);
        final double[] planarAlbedos = spectralAlbedosPolluted.getSpectralAlbedos()[1];
        assertNotNull(planarAlbedos);

        // test results provided by AK in 'snow_properties_CHECKOUT.txt' (20180124), first line:
        assertEquals(OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length, sphericalAlbedos.length);
        assertEquals(OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length, planarAlbedos.length);

        assertEquals(0.6377, sphericalAlbedos[0], 1.E-4);    // 400nm
        assertEquals(0.6658, sphericalAlbedos[3], 1.E-4);    // 490nm
        assertEquals(0.6710, sphericalAlbedos[4], 1.E-4);    // 510nm
        assertEquals(0.6826, sphericalAlbedos[5], 1.E-4);    // 560nm
        assertEquals(0.6934, sphericalAlbedos[6], 1.E-4);    // 620nm
        assertEquals(0.6988, sphericalAlbedos[7], 1.E-4);    // 665nm  // interpolated
        assertEquals(0.6676, sphericalAlbedos[16], 1.E-4);   // 865nm   // interpolated
        assertEquals(0.6440, sphericalAlbedos[17], 1.E-4);   // 885nm   // interpolated
        assertEquals(0.6344, sphericalAlbedos[18], 1.E-4);   // 900nm
        assertEquals(0.6166, sphericalAlbedos[19], 1.E-4);   // 940nm
        assertEquals(0.4512, sphericalAlbedos[20], 1.E-4);   // 1020nm

        assertEquals(0.6202, planarAlbedos[0], 1.E-4);      // 400nm
        assertEquals(0.6493, planarAlbedos[3], 1.E-4);      // 490nm
        assertEquals(0.6547, planarAlbedos[4], 1.E-4);      // 510nm
        assertEquals(0.6666, planarAlbedos[5], 1.E-4);      // 560nm
        assertEquals(0.6779, planarAlbedos[6], 1.E-4);      // 620nm
        assertEquals(0.6835, planarAlbedos[7], 1.E-4);      // 665nm // interpolated
        assertEquals(0.6511, planarAlbedos[16], 1.E-4);     // 865nm  // interpolated
        assertEquals(0.6267, planarAlbedos[17], 1.E-4);     // 885nm  // interpolated
        assertEquals(0.6167, planarAlbedos[18], 1.E-4);     // 900nm
        assertEquals(0.5984, planarAlbedos[19], 1.E-4);     // 940nm
        assertEquals(0.4295, planarAlbedos[20], 1.E-4);     // 1020nm

    }
}

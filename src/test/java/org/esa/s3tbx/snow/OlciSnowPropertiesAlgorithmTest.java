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


public class OlciSnowPropertiesAlgorithmTest {

    @Test
    public void testPolynominalCurveFittingWarrenBrandt() {
        double[] initialGuess = {0., 0., 0., 0., 0., 0.};

        PolynomialFitter curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());

        curveFitter.addObservedPoint(1.02, 2.25E-6);
        curveFitter.addObservedPoint(1.03, 2.33E-6);
        curveFitter.addObservedPoint(1.04, 2.33E-6);
        curveFitter.addObservedPoint(1.05, 2.17E-6);
        curveFitter.addObservedPoint(1.06, 1.96E-6);
        curveFitter.addObservedPoint(1.07, 1.81E-6);
        curveFitter.addObservedPoint(1.08, 1.74E-6);
        curveFitter.addObservedPoint(1.09, 1.73E-6);
        curveFitter.addObservedPoint(1.1, 1.7E-6);
        curveFitter.addObservedPoint(1.11, 1.76E-6);
        curveFitter.addObservedPoint(1.12, 1.82E-6);
        curveFitter.addObservedPoint(1.13, 2.04E-6);
        curveFitter.addObservedPoint(1.14, 2.25E-6);
        curveFitter.addObservedPoint(1.15, 2.29E-6);
        curveFitter.addObservedPoint(1.16, 3.04E-6);
        curveFitter.addObservedPoint(1.17, 3.84E-6);
        curveFitter.addObservedPoint(1.18, 4.77E-6);
        curveFitter.addObservedPoint(1.19, 5.76E-6);
        curveFitter.addObservedPoint(1.2, 6.71E-6);
        curveFitter.addObservedPoint(1.21, 8.66E-6);
        curveFitter.addObservedPoint(1.22, 1.02E-5);
        curveFitter.addObservedPoint(1.23, 1.13E-5);
        curveFitter.addObservedPoint(1.24, 1.22E-5);

        final double[] fit = curveFitter.fit(initialGuess);
        for (int i = 0; i < fit.length; i++) {
            System.out.printf("5th order fit 1.02-1.24: %d,%s%n", i, fit[i]);
        }

    }

    @Test
    public void testPolynominalCurveFitting() {
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
    public void testSigmoidalCurveFitting() {

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
    public void testExp4ParamCurveFitting() {

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

        final double[] spectralSphericalAlbedos = new double[OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length];
        final Exp4ParamFunction exp4ParamFunction = new Exp4ParamFunction();
        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];
            spectralSphericalAlbedos[i] = exp4ParamFunction.value(wvl, fit);
        }

        for (int i = 0; i < 32; i++) {
            final double wvl = 0.4 + 0.02 * i;
            System.out.println("Exp4Param fit: " + wvl + ", " + exp4ParamFunction.value(wvl, fit));
        }

        System.out.println();
    }

    @Test
    public void testExp4ParamCurveFitting_aoki() {

        // this test refers to AK TN 20171010 !!!

        Exp4ParamFitter curveFitter = new Exp4ParamFitter(new LevenbergMarquardtOptimizer());

        // put AOKI values in here and just do the fit...
        curveFitter.addObservedPoint(0.4, 0.975);
        curveFitter.addObservedPoint(0.49, 0.975);
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
    public void testComputeBroadbandAlbedo() {

        final double sza = 36.9;
        final double vza = 3.08;

        double brr400 = 0.68;
        double brr560 = 0.8623;
        double brr865 = 0.7378;
        double brr1020 = 0.4087;

        brr400 *= 0.9798;
        brr560 *= 0.9892;
        brr1020 *= 0.914;

        final double[] brr = new double[]{brr400, brr560, brr865, brr1020};

        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final RefractiveIndexTable refractiveIndexTable = new RefractiveIndexTable();
        refractiveIndexTable.readTableFromFile();
        final SolarSpectrumExtendedTable solarSpectrumExtendedTable = new SolarSpectrumExtendedTable();
        solarSpectrumExtendedTable.readTableFromFile();

        RefractiveIndexTable refractiveIndexInterpolatedTable =
                SnowUtils.getRefractiveIndexInterpolated(refractiveIndexTable,
                                                         solarSpectrumExtendedTable);

        double[][] broadbandAlbedo_simpson =
                OlciSnowPropertiesAlgorithm.computeBroadbandAlbedo(brr, false,
                                                                   refractiveIndexInterpolatedTable,
                                                                   solarSpectrumExtendedTable,
                                                                   sza, vza, 0.0);
        for (double sphericalBroadbandAlbedo : broadbandAlbedo_simpson[0]) {
            System.out.println("spherical broadbandAlbedo_simpson = " + sphericalBroadbandAlbedo);
        }
        for (double planarBroadbandAlbedo : broadbandAlbedo_simpson[1]) {
            System.out.println("planar broadbandAlbedo_simpson = " + planarBroadbandAlbedo);
        }

        System.out.println();
    }

    @Test
    public void testComputeBroadbandAlbedoPolluted() {

        // from Lautaret_pixel1.xls (ML, 20181123)

        final double sza = 36.9;
        final double vza = 3.08;

        double brr400 = 0.72778;
        double brr560 = 0.88665;
        double brr865 = 0.74993;
        double brr1020 = 0.41449;

        brr400 *= 0.9798;
        brr560 *= 0.9892;
        brr1020 *= 0.914;

        final double[] brr = new double[]{brr400, brr560, brr865, brr1020};

        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final RefractiveIndexTable refractiveIndexTable = new RefractiveIndexTable();
        refractiveIndexTable.readTableFromFile();
        final SolarSpectrumExtendedTable solarSpectrumExtendedTable = new SolarSpectrumExtendedTable();
        solarSpectrumExtendedTable.readTableFromFile();

        RefractiveIndexTable refractiveIndexInterpolatedTable =
                SnowUtils.getRefractiveIndexInterpolated(refractiveIndexTable,
                                                         solarSpectrumExtendedTable);

        double[][] broadbandAlbedo_simpson =
                OlciSnowPropertiesAlgorithm.computeBroadbandAlbedo(brr, true,
                                                                   refractiveIndexInterpolatedTable,
                                                                   solarSpectrumExtendedTable,
                                                                   sza, vza, 0.0);
        for (double sphericalBroadbandAlbedo : broadbandAlbedo_simpson[0]) {
            System.out.println("NOV 2018 sphericalBroadbandAlbedo_simpson = " + sphericalBroadbandAlbedo);
        }
        for (double planarBroadbandAlbedo : broadbandAlbedo_simpson[1]) {
            System.out.println("NOV 2018 planarBroadbandAlbedo_simpson = " + planarBroadbandAlbedo);
        }

        System.out.println();
    }

    @Test
    public void testComputePlanarSpectralAlbedoPolluted_nov2018() {

        // from Lautaret_pixel1.xls (ML, 20181123)

        final double sza = 36.9;
        final double vza = 3.08;
        final double saa = 143.8;
        final double vaa = 104.593;

        double brr400 = 0.72778;
        double brr560 = 0.88665;
        double brr865 = 0.74993;
        double brr1020 = 0.41449;

//        brr400 *= 0.9798;
//        brr560 *= 0.9892;
//        brr1020 *= 0.914;

        final double[] brr = new double[]{brr400, brr560, brr865, brr1020};

        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final RefractiveIndexTable refractiveIndexTable = new RefractiveIndexTable();
        refractiveIndexTable.readTableFromFile();
        final SolarSpectrumExtendedTable solarSpectrumExtendedTable = new SolarSpectrumExtendedTable();
        solarSpectrumExtendedTable.readTableFromFile();

        RefractiveIndexTable refractiveIndexInterpolatedTable =
                SnowUtils.getRefractiveIndexInterpolated(refractiveIndexTable,
                                                         solarSpectrumExtendedTable);
        final double raa = SnowUtils.getRelAzi(saa, vaa);
        final double r0Thresh =
                OlciSnowPropertiesAlgorithm.computeR0ReflectancePollutionThresh(sza, vza, raa);

        SpectralAlbedoResult spectralAlbedoPolluted_1 =
                OlciSnowPropertiesAlgorithm.computeFineGridSpectralAlbedo(brr,
                                                                          refractiveIndexInterpolatedTable,
                                                                          solarSpectrumExtendedTable,
                                                                          sza, vza, true, r0Thresh);
        for (double sphericalSpectralAlbedo : spectralAlbedoPolluted_1.getSpectralAlbedos()[0]) {
            System.out.println("NOV 2018 sphericalSpectralAlbedo = " + sphericalSpectralAlbedo);
        }
        for (double planarSpectralAlbedo : spectralAlbedoPolluted_1.getSpectralAlbedos()[1]) {
            System.out.println("NOV 2018 planarSpectralAlbedo = " + planarSpectralAlbedo);
        }

        System.out.println();
    }


    @Test
    public void testNewAlgorithmOct2018() {
        final double sza = 36.9;
        final double vza = 3.08;

        double brr400 = 0.68;
        double brr560 = 0.8623;
        double brr865 = 0.7378;
        double brr1020 = 0.4087;

        brr400 *= 0.9798;
        brr560 *= 0.9892;
        brr1020 *= 0.914;

        final double[] brr = new double[]{brr400, brr560, brr865, brr1020};

        final double amu1 = Math.cos(sza * MathUtils.DTOR);
        final double amu2 = Math.cos(vza * MathUtils.DTOR);

        final double[] wvl = new double[]{400., 560., 865., 1020.};
        final double[] akappa = new double[]{2.365e-11, 2.839e-9, 2.3877e-7, 2.25e-6};
        double[] alpha = new double[4];
        for (int i = 0; i < alpha.length; i++) {
            alpha[i] = 4.0 * Math.PI * akappa[i] / wvl[i];
        }

        // KOKHANOVSKY et al. (2018) paper:

        final double consb = 0.3537;
        final double eps1 = 1. / (1. - consb);
        final double eps2 = 1. - eps1;

        // R0=RR00
        final double rr00 = Math.pow(brr865, eps1) * Math.pow(brr1020, eps2);

        final double p1 = Math.log(brr[0] / rr00) * Math.log(brr[0] / rr00);
        final double p2 = Math.log(brr[1] / rr00) * Math.log(brr[1] / rr00);
        final double am = Math.log(p1 / p2) / Math.log(wvl[1] / wvl[0]);

        final double u1 = SnowUtils.computeU(amu1);
        final double u2 = SnowUtils.computeU(amu2);

        final double x = u1 * u1 * u2 * u2 / (rr00 * rr00);

        final double dlina = Math.log(brr[3] / rr00) * Math.log(brr[3] / rr00) / (x * x * alpha[3]);

        // dlina in mm:
        final double AL = 1.E-6 * dlina;
        final double aksi = 16.0 * 1.6 / 2.25;

        // diameter of grains (in mm):
        final double grainDiam = AL / aksi;

        // f (in 1/mm):
        final double SK = wvl[0] / wvl[3];
        final double f = p1 * Math.pow(SK, am) / (x * x * AL);    // todo: AL = al ??

        // results in ice_refl_output_par.dat:
        // write(21,12)  rr00,al,am,f,diamet:
        // 0.1071E+01  0.1255E+02  0.4371E+01  0.9399E-04  0.1103E+01
        assertEquals(1.07, rr00, 1.E-2);
        assertEquals(12.55, AL, 1.E-2);
        assertEquals(4.371, am, 1.E-3);
        assertEquals(0.9399E-04, f, 1.E-6);
        assertEquals(1.103, grainDiam, 1.E-3);

        // reflectance calculation at 4 points:
        double[] ar = new double[4];
        double[] ad = new double[4];
        for (int i = 0; i < ar.length; i++) {
            final double t = alpha[i] * 1.E6 + f * Math.pow(wvl[i] / wvl[3], -am);
            ar[i] = rr00 * Math.exp(-x * Math.sqrt(t * AL));
            ad[i] = 100.0 * (ar[i] - brr[i]) / brr[i];
        }

        final double rssk = 0.4087;
        // ddd=alog(rssk)*alog(rssk)/3.62**2./alfa(4)  /1.e+6:
        final double ddd = Math.log(rssk) * Math.log(rssk) / (3.62 * 3.62 * alpha[3] * 1.E6);

        final double cvf = rr00 / u2;
        final double rplane = Math.pow(brr[3] / rr00, cvf);

        // results in ice_refl_output.dat:
        // write(22,12)ar(1),refk1,ar(2),refk2,ar(3), refk3,ar(4),refk4
        // 0.6662E+00  0.6663E+00  0.8483E+00  0.8530E+00  0.7303E+00  0.7378E+00  0.3729E+00  0.3736E+00
        assertEquals(0.6662, ar[0], 1.E-4);
        assertEquals(0.6663, brr[0], 1.E-4);
        assertEquals(0.8483, ar[1], 1.E-4);
        assertEquals(0.853, brr[1], 1.E-3);
        assertEquals(0.7303, ar[2], 1.E-4);
        assertEquals(0.7378, brr[2], 1.E-4);
        assertEquals(0.3729, ar[3], 1.E-4);
        assertEquals(0.3736, brr[3], 1.E-4);

        // write(22,12) ad1,ad2,ad3,ad4
        // -0.3131E-02 -0.5521E+00 -0.1023E+01 -0.1782E+00
        assertEquals(-0.3131E-02, ad[0], 1.E-4);
        assertEquals(-0.5521, ad[1], 1.E-4);
        assertEquals(-0.1023E+01, ad[2], 1.E-3);
        assertEquals(-0.1782, ad[3], 1.E-4);

        // write(22,12) ddd
        // 0.2204E+01
        assertEquals(0.2204E+01, ddd, 1.E-4);
        // write(22,12) rplane,cvf,rr00
        // 0.4157E+00  0.8336E+00  0.1071E+01
        assertEquals(0.4157, rplane, 1.E-4);
        assertEquals(0.8336, cvf, 1.E-4);

    }

    @Test
    public void testComputeSpectralAlbedosOct2018() {
        // from 'ice_refl_input.dat'  (Lautaret example pixel)
        final double sza = 36.9;
        final double vza = 3.08;
        final double saa = 143.8;
        final double vaa = 104.593;

        double brr400 = 0.68;
        double brr560 = 0.8623;
        double brr865 = 0.7378;
        double brr1020 = 0.4087;

        brr400 *= 0.9798;
        brr560 *= 0.9892;
        brr1020 *= 0.914;

        final double[] brr = new double[]{brr400, brr560, brr865, brr1020};

        // clean snow
        SpectralAlbedoResult spectralAlbedoResult =
                OlciSnowPropertiesAlgorithm.computeCoarseGridSpectralAlbedoWithErrorEstimates(brr, 0.01, sza, vza, 0.0, false);

        assertNotNull(spectralAlbedoResult.getSpectralAlbedos());

        // polluted snow
        final double raa = SnowUtils.getRelAzi(saa, vaa);
        final double r0Thresh =
                OlciSnowPropertiesAlgorithm.computeR0ReflectancePollutionThresh(sza, vza, raa);
        spectralAlbedoResult =
                OlciSnowPropertiesAlgorithm.computeCoarseGridSpectralAlbedoWithErrorEstimates(brr, 0.01, sza, vza, r0Thresh, true);

        assertNotNull(spectralAlbedoResult.getSpectralAlbedos());

    }

    @Test
    public void testComputeSpectralAlbedosOct2018_2() {
        // from 'ice_refl_input_20181022.dat'   (Dome-C example pixel)
        final double sza = 66.37;
        final double vza = 9.7663;
        final double saa = 143.8;
        final double vaa = 104.593;

        double brr400 = 1.0496;
        double brr560 = 0.9439;
        double brr865 = 0.8263;
        double brr1020 = 0.6728;

        brr400 *= 0.9798;
        brr560 *= 0.9892;
        brr1020 *= 0.914;

        final double[] brr = new double[]{brr400, brr560, brr865, brr1020};

        // clean snow
        SpectralAlbedoResult spectralAlbedoResult =
                OlciSnowPropertiesAlgorithm.computeCoarseGridSpectralAlbedoWithErrorEstimates(brr, 0.01, sza, vza, 0.0, false);

        assertNotNull(spectralAlbedoResult.getSpectralAlbedos());

        // polluted snow
        final double raa = SnowUtils.getRelAzi(saa, vaa);
        final double r0Thresh =
                OlciSnowPropertiesAlgorithm.computeR0ReflectancePollutionThresh(sza, vza, raa);
        spectralAlbedoResult =
                OlciSnowPropertiesAlgorithm.computeCoarseGridSpectralAlbedoWithErrorEstimates(brr, 0.01, sza, vza, r0Thresh, true);

        assertNotNull(spectralAlbedoResult.getSpectralAlbedos());

    }
}

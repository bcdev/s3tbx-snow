package org.esa.s3tbx.snow;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.fitting.PolynomialFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.esa.s3tbx.snow.math.SigmoidalFitter;
import org.esa.s3tbx.snow.math.SigmoidalFunction;
import org.esa.s3tbx.snow.math.lma.LMA;
import org.esa.s3tbx.snow.math.lma.implementations.Polynomial;
import org.esa.snap.core.util.math.MathUtils;

/**
 * Snow Albedo Algorithm for OLCI following A. Kokhanovsky (EUMETSAT).
 * <p/>
 * References:
 * - [TN1] Snow planar broadband albedo determination from spectral OLCI snow spherical albedo measurements. (2017)
 * - [TN2] Snow spectral albedo determination using OLCI measurement. (2017)
 * <p/>
 * todo: extract magic numbers as constants
 *
 * @author olafd
 */
public class OlciSnowAlbedoAlgorithm {

    /**
     * Computes snow spectral albedo for considered wavelengths from input reflectances (after AC) and given geometry.
     * Initially follows: 'snow_albedo_algorithm_1.docx' (AK, 20170519)
     * <p/>
     * The approach requires OLCI BRR bands 1, 2, 3, 4, 5, 12, 17, 21
     * <p/>
     * Update 20170922: 'An update on the algorithm to retrieve snow spectral albedo using OLCI measurements
     * over fresh snow layers (no pollution)', in: 'alex_sept22_2017.pdf'
     * Update 20170929: 'An update on the algorithm to retrieve snow grain size and snow spectral  albedo using OLCI
     * measurements over fresh snow layers (no pollution)', in: 'spectral_albedo_exp_eq.docx'
     *
     * @param brr - Rayleigh corrected reflectance at considered wavelengths (see {@link OlciSnowAlbedoConstants })
     * @param sza - sun zenith angle
     * @param vza - view zenith angle
     * @return array of snow spectral albedos
     */
    static double[] computeSpectralSphericalAlbedos(double[] brr, double sza, double vza,
                                                    SpectralAlbedoMode mode) {
//        double[] spectralSphericalAlbedos = new double[brr.length];
        double[] spectralSphericalAlbedos = new double[OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length];

        // ****************************************************************************************
        // visible subrange (400-510nm, 5 bands):
        final double[] brrVis = new double[5];
        System.arraycopy(brr, 0, brrVis, 0, 5);
        double[] spectralSphericalAlbedosFromBrrVis = new double[4];

        // step 1): Use Eq. (2) at channels: 1( 400nm), 12 (753.75nm), 17(865nm), and 21 (1020nm)
        spectralSphericalAlbedosFromBrrVis[0] = computeSpectralAlbedoFromBrrVis(brr[0], brrVis, sza, vza);
        spectralSphericalAlbedosFromBrrVis[1] = computeSpectralAlbedoFromBrrVis(brr[5], brrVis, sza, vza);
        spectralSphericalAlbedosFromBrrVis[2] = computeSpectralAlbedoFromBrrVis(brr[6], brrVis, sza, vza);
        final double brr21 = brr[Sensor.OLCI.getRequiredBrrBandNames().length - 1];
        spectralSphericalAlbedosFromBrrVis[3] = computeSpectralAlbedoFromBrrVis(brr21, brrVis, sza, vza);

        // AK 20170922:
        // "An update on the algorithm to retrieve snow spectral albedo1020 using OLCI measurements
        // over fresh snow layers (no pollution)":
        if (mode == SpectralAlbedoMode.POLYNOMINAL_FIT) {

            // OD proposed also to provide polynominal fit
            double[] initialGuess = {0., 0., 0., 0., 0., 0., 0., 0.};
            PolynomialFitter curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());

            curveFitter.clearObservations();
            curveFitter.addObservedPoint(0.4, spectralSphericalAlbedosFromBrrVis[0]);
            curveFitter.addObservedPoint(0.753, spectralSphericalAlbedosFromBrrVis[1]);
            curveFitter.addObservedPoint(0.865, spectralSphericalAlbedosFromBrrVis[2]);
            curveFitter.addObservedPoint(1.02, spectralSphericalAlbedosFromBrrVis[3]);
            double[] fit = ((PolynomialFitter) curveFitter).fit(initialGuess);

            // use eq. (3) simplified to 2 parameters:
            for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
                final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
                spectralSphericalAlbedos[i] = new PolynomialFunction(fit).value(wvl);
            }

            // ****************************************************************************************
        } else if (mode == SpectralAlbedoMode.SIGMOIDAL_FIT) {
            // step 2): Fit the derived values of albedo1020
            SigmoidalFitter curveFitter = new SigmoidalFitter(new LevenbergMarquardtOptimizer());

            double[] initialGuess = {1., 1.};
            final SigmoidalFunction sigmoidalFunction2 = new SigmoidalFunction(2);
//            double[] initialGuess = {1., 1., 1., 1.};
//            double[] fit = curveFitter.fit(initialGuess, 4);
            curveFitter.clearObservations();
            curveFitter.addObservedPoint(0.4, spectralSphericalAlbedosFromBrrVis[0]);
            curveFitter.addObservedPoint(0.753, spectralSphericalAlbedosFromBrrVis[1]);
            curveFitter.addObservedPoint(0.865, spectralSphericalAlbedosFromBrrVis[2]);
            curveFitter.addObservedPoint(1.02, spectralSphericalAlbedosFromBrrVis[3]);
            double[] fit = curveFitter.fit(initialGuess, 2);    // todo: this is the performance killer!!
            // use eq. (3) simplified to 2 parameters:
            for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
                final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
                spectralSphericalAlbedos[i] = sigmoidalFunction2.value(wvl, fit);
            }
            // ****************************************************************************************
        } else if (mode == SpectralAlbedoMode.EXPONENTIAL_SQRT_FIT) {

            // latest approach, AK 20170929: "spectral_albedo_exp_eq.doc":
            final double albedo1020 = spectralSphericalAlbedosFromBrrVis[3];
            double grainDiameter = computeGrainDiameter(albedo1020);
            final double b = 3.62;
            final double[] a = new double[spectralSphericalAlbedos.length];
            for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
                final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
                a[i] = computeA(wvl);
                spectralSphericalAlbedos[i] = Math.exp(-b * Math.sqrt(a[i] * grainDiameter));
            }
        } else if (mode == SpectralAlbedoMode.EXPONENTIAL_3PARAM_FIT) {
            // todo

        } else {
            throw new IllegalArgumentException("spectral albedo algoritm mode " + mode.getName() + " not supported");
        }

        // *****************************************************************************************

        return spectralSphericalAlbedos;
    }

    static double[] computeSpectralSphericalAlbedos_test(double[] brr, double sza, double vza) {
//        double[] spectralSphericalAlbedos = new double[brr.length];
        double[] spectralSphericalAlbedos = new double[OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length];

        // ****************************************************************************************
        // visible subrange (400-510nm, 5 bands):
        final double[] brrVis = new double[5];
        System.arraycopy(brr, 0, brrVis, 0, 5);
        double[] spectralSphericalAlbedosFromBrrVis = new double[4];

        // step 1): Use Eq. (2) at channels: 1( 400nm), 12 (753.75nm), 17(865nm), and 21 (1020nm)
        spectralSphericalAlbedosFromBrrVis[0] = computeSpectralAlbedoFromBrrVis(brr[0], brrVis, sza, vza);
        spectralSphericalAlbedosFromBrrVis[1] = computeSpectralAlbedoFromBrrVis(brr[5], brrVis, sza, vza);
        spectralSphericalAlbedosFromBrrVis[2] = computeSpectralAlbedoFromBrrVis(brr[6], brrVis, sza, vza);
        final double brr21 = brr[Sensor.OLCI.getRequiredBrrBandNames().length - 1];
        spectralSphericalAlbedosFromBrrVis[3] = computeSpectralAlbedoFromBrrVis(brr21, brrVis, sza, vza);

        SigmoidalFitter curveFitter = new SigmoidalFitter(new LevenbergMarquardtOptimizer());

        double[] initialGuess = {1., 1.};
        final SigmoidalFunction sigmoidalFunction_2 = new SigmoidalFunction(2);
        curveFitter.clearObservations();
        curveFitter.addObservedPoint(0.4, spectralSphericalAlbedosFromBrrVis[0]);
        curveFitter.addObservedPoint(0.753, spectralSphericalAlbedosFromBrrVis[1]);
        curveFitter.addObservedPoint(0.865, spectralSphericalAlbedosFromBrrVis[2]);
        curveFitter.addObservedPoint(1.02, spectralSphericalAlbedosFromBrrVis[3]);
        double[] fit = curveFitter.fit(initialGuess, 2);    // todo: this is the performance killer!!
//        double[] fit = initialGuess;
        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            spectralSphericalAlbedos[i] = sigmoidalFunction_2.value(wvl, fit);
        }

        return spectralSphericalAlbedos;
    }

    static double[] computeSpectralSphericalAlbedos_test_lma(double[] brr, double sza, double vza) {
        double[] spectralSphericalAlbedos = new double[OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length];

        // ****************************************************************************************
        // visible subrange (400-510nm, 5 bands):
        final double[] brrVis = new double[5];
        System.arraycopy(brr, 0, brrVis, 0, 5);
        double[] spectralSphericalAlbedosFromBrrVis = new double[4];

        // step 1): Use Eq. (2) at channels: 1( 400nm), 12 (753.75nm), 17(865nm), and 21 (1020nm)
        spectralSphericalAlbedosFromBrrVis[0] = computeSpectralAlbedoFromBrrVis(brr[0], brrVis, sza, vza);
        spectralSphericalAlbedosFromBrrVis[1] = computeSpectralAlbedoFromBrrVis(brr[5], brrVis, sza, vza);
        spectralSphericalAlbedosFromBrrVis[2] = computeSpectralAlbedoFromBrrVis(brr[6], brrVis, sza, vza);
        final double brr21 = brr[Sensor.OLCI.getRequiredBrrBandNames().length - 1];
        spectralSphericalAlbedosFromBrrVis[3] = computeSpectralAlbedoFromBrrVis(brr21, brrVis, sza, vza);

        LMA lma = new LMA(
                new Polynomial(),
//                new double[]{1, 1, 1., 1., 1., 1., 1., 1.},
                new double[]{0., 0., 0., 0., 0., 0., 0., 0.},
                new double[][]{
                        {0.4, 0.753, 0.865, 1.02},
                        {spectralSphericalAlbedosFromBrrVis[0],
                                spectralSphericalAlbedosFromBrrVis[1],
                                spectralSphericalAlbedosFromBrrVis[2],
                                spectralSphericalAlbedosFromBrrVis[3]}}
        );
        lma.fit();
        double[] fit = lma.parameters;

        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            spectralSphericalAlbedos[i] = new Polynomial().getY(wvl, fit);
        }

        return spectralSphericalAlbedos;
    }


    /**
     * New algorithm from 'alex_sept_22_2017.pdf':
     *
     * @param brr    - brr value
     * @param brrVis - array of BRR values in VIS range (400-510nm)
     * @param sza    - SZA
     * @param vza    - VZA
     * @return spectral albedo
     */
    static double computeSpectralAlbedoFromBrrVis(double brr, double[] brrVis, double sza, double vza) {
        // eq. (1):
        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);

        final double k_mu_0 = computeK(mu_0);
        final double k_mu = computeK(mu);

        // get R_v for eq. (2):
        double brrMax = Double.MIN_VALUE;
        for (double vis : brrVis) {
            if (vis > brrMax) {
                brrMax = vis;
            }
        }

        // eq. (2):
        final double x = k_mu_0 * k_mu / brrMax;

        return Math.pow(brr / brrMax, 1.0 / x);
    }

    /**
     * Computes planar from spherical albedos at considered wavelengths.
     * Follows 'snow_albedo_algorithm_1.docx' (AK, 20170519).
     *
     * @param sphericalAlbedos - the spherical albedos at considered wavelengths
     * @param sza              - sun zenith angle
     * @return array of planar albedos
     */
    static double[] computePlanarFromSphericalAlbedos(double[] sphericalAlbedos, double sza) {
        double[] planarAlbedos = new double[sphericalAlbedos.length];
        for (int i = 0; i < planarAlbedos.length; i++) {
            planarAlbedos[i] = computePlanarFromSphericalAlbedo(sphericalAlbedos[i], sza);
        }
        return planarAlbedos;
    }


    /**
     * Computes planar from spherical albedos at considered wavelengths.
     * Follows 'snow_albedo_algorithm_1.docx' (AK, 20170519)
     *
     * @param sphericaAlbedo - the spherical albedo value
     * @param sza            - sun zenith angle
     * @return theplanar albedo value
     */
    static double computePlanarFromSphericalAlbedo(double sphericaAlbedo,
                                                   double sza) {
        // eq. (5):
        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        return Math.pow(sphericaAlbedo, computeK(mu_0));
    }


    /**
     * Computes spherical broadband albedo terms (b_1, b_2 integrals, grain diameter).
     * Follows 'algorithm__BROADBAND_SPHERICAL_ALBEDO.docx' (AK, 20170530)
     *
     * @param sphericalSpectralAlbedos - the spherical spectral albedos at considered wavelengths.
     * @return SphericalBroadbandAlbedo object holding b_1, b_2 integrals and grain diameter.
     */
    static SphericalBroadbandAlbedo computeSphericalBroadbandAlbedoTerms(double[] sphericalSpectralAlbedos) {
        SphericalBroadbandAlbedo sbba = new SphericalBroadbandAlbedo();
        sbba.setR_b1(integrateR_b1(sphericalSpectralAlbedos));
        final double albedo1020 = sphericalSpectralAlbedos[sphericalSpectralAlbedos.length - 1];
        final double grainDiameter = computeGrainDiameter(albedo1020);
        sbba.setGrainDiameter(grainDiameter);
        sbba.setR_b2(integrateR_b2(grainDiameter));

        return sbba;
    }

    /**
     * Computes the snow grain diameter for given Rayleigh corrected reflectance at band 21 (1020nm).
     * Follows 'algorithm__BROADBAND_SPHERICAL_ALBEDO.docx' (AK, 20170530)
     *
     * @param albedo1020 - Spectral albedo at band 21 (1020nm)
     * @return the snow grain diameter in microns
     */
    static double computeGrainDiameter(double albedo1020) {
        // eq. (5):
        final double b = 3.62;
        final double a_21 = computeA(OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[20]);
        return Math.log(albedo1020) * Math.log(albedo1020) / (b * b * a_21);    // this is the grain diameter in microns!!
    }

    private static double computeA(double lambda) {
        final double chi = 2.44E-13 * Math.exp(lambda / 0.06367);
        return 4.0 * Math.PI * chi / lambda;
    }

    /**
     * Computes the 'r_b1' term for broadband snow albedo
     * Follows 'algorithm__BROADBAND_SPHERICAL_ALBEDO.docx' (AK, 20170530)
     *
     * @param spectralAlbedos - the spectral albedos at considered wavelengths.
     * @return r_b1
     */
    static double integrateR_b1(double[] spectralAlbedos) {
        double r_b1 = 0.0;
        // interpolate input spectralAlbedos (21 OLCI wavelengths) to full grid 300-1020nm (53 wavelengths)
        final double[] wvlFull = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI_EXTENDED;
        final double[] wvlOlci = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI;
        // we need to extrapolate array to 300nm on the lower side:
        double[] wvl22 = new double[wvlOlci.length + 1];
        wvl22[0] = 0.299;
        System.arraycopy(wvlOlci, 0, wvl22, 1, wvlOlci.length);
        double[] spectralAlbedos22 = new double[spectralAlbedos.length + 1];
        spectralAlbedos22[0] = spectralAlbedos[0];
        System.arraycopy(spectralAlbedos, 0, spectralAlbedos22, 1, spectralAlbedos.length);
        final double[] spectralAlbedosInterpolated = interpolateSpectralAlbedos(wvl22, spectralAlbedos22, wvlFull);

        // integration: eq. (A.1) with f_lambda from Table (A.2)
        for (int i = 0; i < OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI_EXTENDED.length; i++) {
            r_b1 += spectralAlbedosInterpolated[i] * OlciSnowAlbedoConstants.F_LAMBDA_EXTENDED[i];
        }
        return r_b1;
    }

    /**
     * Interpolates array of spectral albedos (here: 21 OLCI) to coarser grid (here 53 wavelengths 300-1020nm).
     * Follows 'algorithm__BROADBAND_SPHERICAL_ALBEDO.docx' (AK, 20170530)
     *
     * @param x  - array of x values (OLCI wavelengths)
     * @param y  - array of spectral albedos at OLCI wavelengths
     * @param xi - array of 53 wavelengths 300-1200nm
     * @return array of spectral albedos at 53 wavelengths 300-1200nm
     */
    static double[] interpolateSpectralAlbedos(double[] x, double[] y, double[] xi) {
        final LinearInterpolator linearInterpolator = new LinearInterpolator();
        PolynomialSplineFunction psf = linearInterpolator.interpolate(x, y);

        double[] yi = new double[xi.length];
        for (int i = 0; i < xi.length; i++) {
            yi[i] = psf.value(xi[i]);
        }
        return yi;
    }

    /**
     * Computes the 'r_b2' term for broadband snow albedo
     * Follows 'algorithm__BROADBAND_SPHERICAL_ALBEDO.docx' (AK, 20170530)
     *
     * @param d - the snow grain diameter.
     * @return r_b2
     */
    static double integrateR_b2(double d) {
        // eq. (4):
        final double q0 = 0.0947;
        final double q1 = 0.0569;
        final double d0 = 200.0;
        return q0 - q1 * Math.log10(d / d0);
    }

    private static double computeK(double mu) {
        return 3.0 * (1.0 + 2.0 * mu) / 7.0;
    }

    /**
     * Container object for spherical broadband albedo terms (b_1, b_2 integrals, grain diameter).
     */
    static class SphericalBroadbandAlbedo {
        double r_b1;
        double r_b2;
        double grainDiameter;

        public double getR_b1() {
            return r_b1;
        }

        public void setR_b1(double r_b1) {
            this.r_b1 = r_b1;
        }

        public double getR_b2() {
            return r_b2;
        }

        public void setR_b2(double r_b2) {
            this.r_b2 = r_b2;
        }

        public double getGrainDiameter() {
            return grainDiameter;
        }

        public void setGrainDiameter(double grainDiameter) {
            this.grainDiameter = grainDiameter;
        }
    }
}

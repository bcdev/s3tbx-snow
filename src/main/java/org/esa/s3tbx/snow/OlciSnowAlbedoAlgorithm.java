package org.esa.s3tbx.snow;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.esa.s3tbx.snow.math.Exp4ParamFitter;
import org.esa.s3tbx.snow.math.Exp4ParamFunction;
import org.esa.s3tbx.snow.math.SigmoidalFitter;
import org.esa.s3tbx.snow.math.SigmoidalFunction;
import org.esa.snap.core.util.math.MathUtils;

/**
 * Snow Albedo Algorithm for OLCI following A. Kokhanovsky (former EUMETSAT, now Vitrociset_Belgium).
 * <p>
 * Initial references, updated several times:
 * - [TN1] Snow planar broadband albedo determination from spectral OLCI snow spherical albedo measurements. (2017)
 * - [TN2] Snow spectral albedo determination using OLCI measurement. (2017)
 * <p>
 *
 * @author olafd
 */
class OlciSnowAlbedoAlgorithm {

    /**
     * Computes spectral spherical and planar albedos using given computation mode from the ones proposed by AK.
     * Currently we only use SIMPLE_APPROXIMATION (20171120).
     *
     * @param brr                 - subset of BRR spectrum
     * @param sza                 - sun zenith angle (deg)
     * @param vza                 - view zenith angle (deg)
     * @param referenceWavelength - OLCI reference wavelength (1020 or 865 nm)
     * @param mode                - computation mode  @return double[][]{spectralSphericalAlbedo, spectralPlanarAlbedo}
     */
    static double[][] computeSphericalAlbedos(double[] brr, double sza, double vza,
                                              double referenceWavelength,
                                              SpectralAlbedoMode mode) {
        final int numWvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length;
        double[][] sphericalAlbedos = new double[2][numWvl];

        if (mode == SpectralAlbedoMode.SIMPLE_APPROXIMATION) {
            computeSpectralSphericalAlbedoWithSimpleApproximation(brr, sza, vza, numWvl,
                                                                  referenceWavelength, sphericalAlbedos);
        } else {
            // we need the visible subrange (400-510nm, 5 bands):
            final double[] brrVis = new double[5];
            System.arraycopy(brr, 0, brrVis, 0, 5);
            double[] spectralSphericalAlbedosFromBrrVis = new double[4];

            // step 1): Use Eq. (2) at channels: 1( 400nm), 12 (753.75nm), 17(865nm), and 21 (1020nm)
            spectralSphericalAlbedosFromBrrVis[0] = computeSpectralAlbedoFromBrrVis(brr[0], brrVis, sza, vza);
            spectralSphericalAlbedosFromBrrVis[1] = computeSpectralAlbedoFromBrrVis(brr[5], brrVis, sza, vza);
            spectralSphericalAlbedosFromBrrVis[2] = computeSpectralAlbedoFromBrrVis(brr[6], brrVis, sza, vza);
            final double brr21 = brr[brr.length - 1];
            spectralSphericalAlbedosFromBrrVis[3] = computeSpectralAlbedoFromBrrVis(brr21, brrVis, sza, vza);

            if (mode == SpectralAlbedoMode.POLYNOMINAL_FIT) {
                computeSpectralSphericalAlbedoWithPolynominalFit(sphericalAlbedos[0], spectralSphericalAlbedosFromBrrVis);
            } else if (mode == SpectralAlbedoMode.SIGMOIDAL_FIT) {
                computeSpectralSphericalAlbedoWithSigmoidalFit(sphericalAlbedos[0], spectralSphericalAlbedosFromBrrVis);
                // ****************************************************************************************
            } else if (mode == SpectralAlbedoMode.EXPONENTIAL_SQRT_FIT) {
                computeSpectralSphericalAlbedoWithExponentialSqrtFit(sphericalAlbedos[0], spectralSphericalAlbedosFromBrrVis[3]);

            } else if (mode == SpectralAlbedoMode.EXPONENTIAL_4PARAM_FIT) {
                computeSpectralSphericalAlbedoWithExponential4ParamFit(brr, sza, vza, sphericalAlbedos[0],
                                                                       spectralSphericalAlbedosFromBrrVis);
            } else {
                throw new IllegalArgumentException("spectral albedo algoritm mode " + mode.getName() + " not supported");
            }
        }
        sphericalAlbedos[1] = computePlanarFromSphericalAlbedos(sphericalAlbedos[0], sza);

        return sphericalAlbedos;
    }

    /**
     * Computes broadband albedos following AK algorithm given in 'Technical note_BBA_DECEMBER_2017.doc' (20171204).
     *
     * @param mu_0 - mu0
     * @param d - snow grain diameter
     * @param refractiveIndexTable - table with refractive indices
     * @param solarSpectrumTable - table with solar spectrum
     *
     * @return double[]{pbbaVis, pbbaNir, pbbaSw};
     */
    static double[] computeBroadbandAlbedo(double mu_0, double d,
                                                  RefractiveIndexTable refractiveIndexTable,
                                                  SolarSpectrumTable solarSpectrumTable) {
        double pbbaVis;
        double pbbaNir;
        double pbbaSw;
        double numeratorVis = 0.0;
        double denominatorVis = 0.0;
        double numeratorNir = 0.0;
        double denominatorNir = 0.0;
        double numeratorSw = 0.0;
        double denominatorSw = 0.0;

        double[] wvlsFull = solarSpectrumTable.getWvl();
        final double[] f_lambda = SnowUtils.getFLambda(solarSpectrumTable);

        for (int i = 0; i < wvlsFull.length - 1; i++) {
            final double wvl = wvlsFull[i];
            double u;
            if (mu_0 == 1.0) {
                u = 1.0;
            } else {
                u = SnowUtils.computeU(mu_0);
            }
            final double chi = refractiveIndexTable.getRefractiveIndexImag(i);
            final double k = 4.0 * Math.PI * chi / wvl;
            final double planarSpectralAlbedo = Math.exp(-3.6 * u * Math.sqrt(k * d));

            double dx = wvlsFull[i + 1] - wvl;
            // VIS 0.3-0.7
            if (wvl > OlciSnowAlbedoConstants.BB_WVL_1 && wvl < OlciSnowAlbedoConstants.BB_WVL_2) {
                numeratorVis += planarSpectralAlbedo * f_lambda[i] * dx;
                denominatorVis += f_lambda[i] * dx;
            }
            // NIR 0.7-2.4
            if (wvl > OlciSnowAlbedoConstants.BB_WVL_2 && wvl < OlciSnowAlbedoConstants.BB_WVL_3) {
                numeratorNir += planarSpectralAlbedo * f_lambda[i] * dx;
                denominatorNir += f_lambda[i] * dx;
            }
            // SW 0.3-2.4
            if (wvl > OlciSnowAlbedoConstants.BB_WVL_1 && wvl < OlciSnowAlbedoConstants.BB_WVL_3) {
                numeratorSw += planarSpectralAlbedo * f_lambda[i] * dx;
                denominatorSw += f_lambda[i] * dx;
            }
        }

        pbbaVis = numeratorVis / denominatorVis;
        pbbaNir = numeratorNir / denominatorNir;
        pbbaSw = numeratorSw / denominatorSw;

        return new double[]{pbbaVis, pbbaNir, pbbaSw};
    }


    static double[] computeSpectralPpa(double[] brr, double sza, double vza) {
        double[] ppa = new double[brr.length];

        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);
        final double k_mu_0 = SnowUtils.computeU(mu_0);
        final double k_mu = SnowUtils.computeU(mu);
        final double brr_400 = Math.min(1.0, brr[0]);
        final double m = k_mu_0 * k_mu/brr_400;

        for (int i = 0; i < brr.length; i++) {
            final double y = Math.log(brr_400/brr[i]) / m;
            ppa[i] = 3.0*y*y/64.0;
        }

        return ppa;
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
        return Math.pow(sphericaAlbedo, SnowUtils.computeU(mu_0));
    }

    /**
     * Computes the snow grain diameter for given Rayleigh corrected reflectance at band 21 (1020nm).
     * Follows 'sgs_nov_20_865nm.doc' (AK, 20171120)
     *
     * @param refAlbedo - Spectral albedo at band 18 or 21 (865 or 1020nm)
     * @param refWvl - Reference wavelength (865 or 1020nm)
     * @return the snow grain diameter in microns
     */
    static double computeGrainDiameter(double refAlbedo, double refWvl) {
        final double b = 3.62;
        final int numWvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length;
        final double chiRef = refWvl == 1020.0 ? OlciSnowAlbedoConstants.ICE_REFR_INDEX[numWvl - 1] :
                OlciSnowAlbedoConstants.ICE_REFR_INDEX[numWvl - 5];
        final double kappaRef = 4.0 * Math.PI * chiRef / (refWvl/1000.0);  // eq. (4)
        return Math.log(refAlbedo) * Math.log(refAlbedo) / (b * b * kappaRef);
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

    /**
     * Computation of spectral spherical albedo following AK new approach, 20171120.
     * Currently used as default.
     * <p>
     * References:
     * [1]: The simple approximation  for the spectral planar albedo. TN AK, 20171120. File: nov_20.doc.
     * [2]: The snow grain size determination. TN AK, 20171120. File: sgs_nov_20.doc.
     */
    private static void computeSpectralSphericalAlbedoWithSimpleApproximation(double[] brr,
                                                                              double sza, double vza, int numWvl,
                                                                              double refWvl,
                                                                              double[][] sphericalAlbedos) {
        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);
        final double k_mu_0 = SnowUtils.computeU(mu_0);
        final double k_mu = SnowUtils.computeU(mu);
        final double brr_400 = Math.min(1.0, brr[0]);
        final double xi = brr_400 / (k_mu_0 * k_mu);   // [1], eq. (5)

        final double brr_ref = brr[brr.length - 1];
        final double r_ref = Math.pow(brr_ref / brr_400, xi); // [1], eq. (5)
        final double wvl_ref = refWvl/1000.0;  // in microns!
        final double chi_ref = refWvl == 1020.0 ? OlciSnowAlbedoConstants.ICE_REFR_INDEX[numWvl - 1] :
                OlciSnowAlbedoConstants.ICE_REFR_INDEX[numWvl - 5];
        final double kappa_ref = 4.0 * Math.PI * chi_ref / wvl_ref;  // [2], eqs. (1), (2)
        final double B = Math.log(r_ref) * Math.log(r_ref) / kappa_ref;  // [2], eq. (2) transformed

        for (int i = 0; i < numWvl; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            final double chi = OlciSnowAlbedoConstants.ICE_REFR_INDEX[i];
            final double kappa = 4.0 * Math.PI * chi / wvl;  // [2], eqs. (1), (2)

            sphericalAlbedos[0][i] = Math.exp(-Math.sqrt(B * kappa));  // spectral spherical abledo
        }
    }

    private static void computeSpectralSphericalAlbedoWithExponential4ParamFit(double[] brr, double sza, double vza,
                                                                               double[] spectralSphericalAlbedos,
                                                                               double[] spectralSphericalAlbedosFromBrrVis) {
        Exp4ParamFitter curveFitter = new Exp4ParamFitter(new LevenbergMarquardtOptimizer());

        double[] initialGuess = new double[4]; // a, kappa_1, L, b, from AK TN, 20171005
        initialGuess[0] = brr[0]; // ~1.0?!
        initialGuess[1] = 0;
        final double kappa2_1020 = 1.E-6;
        final double r_1020 = spectralSphericalAlbedosFromBrrVis[3];
        initialGuess[2] = 1.02 * Math.log(r_1020) * Math.log(r_1020) / (2.0 * Math.PI * kappa2_1020);
        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);
        final double k_mu_0 = SnowUtils.computeU(mu_0);
        final double k_mu = SnowUtils.computeU(mu);
        initialGuess[3] = k_mu * k_mu_0 / brr[0];
        final Exp4ParamFunction exp4ParamFunction = new Exp4ParamFunction();
        curveFitter.clearObservations();
        curveFitter.addObservedPoint(0.4, spectralSphericalAlbedosFromBrrVis[0]);
        curveFitter.addObservedPoint(0.753, spectralSphericalAlbedosFromBrrVis[1]);
        curveFitter.addObservedPoint(0.865, spectralSphericalAlbedosFromBrrVis[2]);
        curveFitter.addObservedPoint(1.02, spectralSphericalAlbedosFromBrrVis[3]);
        double[] fit = curveFitter.fit(initialGuess);
        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            spectralSphericalAlbedos[i] = exp4ParamFunction.value(wvl, fit);
        }
    }

    private static void computeSpectralSphericalAlbedoWithExponentialSqrtFit(double[] spectralSphericalAlbedos,
                                                                             double spectralSphericalAlbedosFromBrrVis) {
        // latest approach, AK 20170929: "spectral_albedo_exp_eq.doc":
        double grainDiameter = computeGrainDiameter(spectralSphericalAlbedosFromBrrVis, 1020.0);
        final double b = 3.62;
        final double[] a = new double[spectralSphericalAlbedos.length];
        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            final double chi = 2.44E-13 * Math.exp(wvl / 0.06367);
            a[i] = 4.0 * Math.PI * chi / wvl;
            spectralSphericalAlbedos[i] = Math.exp(-b * Math.sqrt(a[i] * grainDiameter));
        }
    }

    private static void computeSpectralSphericalAlbedoWithSigmoidalFit(double[] spectralSphericalAlbedos,
                                                                       double[] spectralSphericalAlbedosFromBrrVis) {
        SigmoidalFitter curveFitter = new SigmoidalFitter(new LevenbergMarquardtOptimizer());

        double[] initialGuess = {1., 1.};
        final SigmoidalFunction sigmoidalFunction2 = new SigmoidalFunction(2);
        curveFitter.clearObservations();
        curveFitter.addObservedPoint(0.4, spectralSphericalAlbedosFromBrrVis[0]);
        curveFitter.addObservedPoint(0.753, spectralSphericalAlbedosFromBrrVis[1]);
        curveFitter.addObservedPoint(0.865, spectralSphericalAlbedosFromBrrVis[2]);
        curveFitter.addObservedPoint(1.02, spectralSphericalAlbedosFromBrrVis[3]);
        double[] fit = curveFitter.fit(initialGuess, 2); // CAREFUL: this is a performance killer in commons-math if observed points contain NaNs!!
        // use eq. (3) simplified to 2 parameters:
        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            spectralSphericalAlbedos[i] = sigmoidalFunction2.value(wvl, fit);
        }
    }

    private static void computeSpectralSphericalAlbedoWithPolynominalFit(double[] spectralSphericalAlbedos,
                                                                         double[] spectralSphericalAlbedosFromBrrVis) {
        // OD proposed also to provide polynominal fit
        double[] initialGuess = {0., 0., 0., 0., 0., 0., 0., 0.};
        PolynomialFitter curveFitter = new PolynomialFitter(new LevenbergMarquardtOptimizer());

        curveFitter.clearObservations();
        curveFitter.addObservedPoint(0.4, spectralSphericalAlbedosFromBrrVis[0]);
        curveFitter.addObservedPoint(0.753, spectralSphericalAlbedosFromBrrVis[1]);
        curveFitter.addObservedPoint(0.865, spectralSphericalAlbedosFromBrrVis[2]);
        curveFitter.addObservedPoint(1.02, spectralSphericalAlbedosFromBrrVis[3]);
        double[] fit = curveFitter.fit(initialGuess);

        // use eq. (3) simplified to 2 parameters:
        for (int i = 0; i < spectralSphericalAlbedos.length; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            spectralSphericalAlbedos[i] = new PolynomialFunction(fit).value(wvl);
        }
    }

    /**
     * New algorithm from 'alex_sept_22_2017.pdf':
     * Computes a spectral albedo in visible range for a given BRR from given BRR visible sub-spectrum.
     *
     * @param brr    - brr value
     * @param brrVis - array of BRR values in VIS range (400-510nm)
     * @param sza    - SZA
     * @param vza    - VZA
     * @return spectral albedo
     */
    private static double computeSpectralAlbedoFromBrrVis(double brr, double[] brrVis, double sza, double vza) {
        // eq. (1):
        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);

        final double k_mu_0 = SnowUtils.computeU(mu_0);
        final double k_mu = SnowUtils.computeU(mu);

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
    private static double[] computePlanarFromSphericalAlbedos(double[] sphericalAlbedos, double sza) {
        double[] planarAlbedos = new double[sphericalAlbedos.length];
        for (int i = 0; i < planarAlbedos.length; i++) {
            planarAlbedos[i] = computePlanarFromSphericalAlbedo(sphericalAlbedos[i], sza);
        }
        return planarAlbedos;
    }

}

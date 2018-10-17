package org.esa.s3tbx.snow;

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
class OlciSnowPropertiesAlgorithm {

    /**
     * Computes spectral spherical and planar albedos using given computation mode from the ones proposed by AK.
     * Currently we only use SIMPLE_APPROXIMATION (20171120).
     *
     * @param brr                 - subset of BRR spectrum  (BRR_01, BRR_05, BRR_21, BRR_17)
     * @param sza                 - sun zenith angle (deg)
     * @param vza                 - view zenith angle (deg)
     * @param referenceWavelength - OLCI reference wavelength (1020 or 865 nm)
     * @param mode                - computation mode  @return double[][]{spectralSphericalAlbedo, spectralPlanarAlbedo}
     * @param useAlgoApril2018 - true if new AK algorithm from April 2018 is used
     * @return double[][] spectralAlbedos : spherical and planar at OLCI wavelengths
     */
    static double[][] computeSpectralAlbedos(double[] brr, double sza, double vza,
                                             double referenceWavelength,
                                             SpectralAlbedoMode mode, boolean useAlgoApril2018) {
        final int numWvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length;
        double[][] spectralAlbedos = new double[2][numWvl];

        if (mode == SpectralAlbedoMode.SIMPLE_APPROXIMATION) {
            if (useAlgoApril2018) {
                computeSpectralAlbedoFromTwoWavelengths(brr, sza, vza, numWvl, spectralAlbedos);
            } else {
                computeSpectralAlbedoWithSimpleApproximation(brr, sza, vza, numWvl,
                                                             referenceWavelength, spectralAlbedos);
            }
        } else {
            // we no longer support this
            throw new IllegalArgumentException("spectral albedo algoritm mode " + mode.getName() + " not supported");
        }

        return spectralAlbedos;
    }

    /**
     * Computes spectral spherical and planar albedos in case of polluted snow.
     * Algorithm Tech Note: 'Manual_31_01_2018.docx', AK 20180131
     *
     * @param brr - subset of BRR spectrum  (BRR_01, BRR_05, BRR_21, BRR_17)
     * @param pollutedSnowParams - parameters for polluted snow,
     *                           see {@link OlciSnowPropertiesAlgorithm#computePollutedSnowParams
     *                           (double, double, double, double, double)}
     * @param sza - sun zenith angle (deg)
     * @param vza - view zenith angle (deg)
     * @param useAlgoApril2018 - true if new AK algorithm from April 2018 is used ('Manual_04_04_2018.docx')
     *
     * @return double[][] spectralAlbedos : spherical and planar at OLCI wavelengths
     */
    static SpectralAlbedoResult computeSpectralAlbedosPolluted(double[] brr,
                                                               double[] pollutedSnowParams,
                                                               double sza,
                                                               double vza,
                                                               double deltaBrr,
                                                               boolean useAlgoApril2018) {
        if (useAlgoApril2018) {
            return computeSpectralAlbedosPollutedFromFourWavelengths(brr, deltaBrr, sza, vza);
        } else {
            return computeSpectralAlbedosPolluted(pollutedSnowParams, sza);
        }
    }

    private static SpectralAlbedoResult computeSpectralAlbedosPolluted(double[] pollutedSnowParams, double sza) {
        final int numWvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length;
        double[][] spectralAlbedos = new double[2][numWvl];

        final double grainDiam = pollutedSnowParams[0] * 1000.0;  // in microns here
        final double soot = pollutedSnowParams[1] * 1.E-6; // dimensionless here

        final double akas = 0.46;
        final double am0 = Math.cos(Math.PI * sza / 180.);
        final double ak0 = (3. / 7.) * (1. + 2. * am0);

        for (int i = 0; i < numWvl; i++) {
            final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];
            final double p = 6. * 4. * Math.PI * akas * soot * grainDiam / wvl;
            final double chi = OlciSnowPropertiesConstants.ICE_REFR_INDEX[i];   // chi = ak(i) in 'soot.f' breadboard
            final double x = Math.sqrt(13. * 4. * Math.PI * grainDiam * chi / wvl + p);
            spectralAlbedos[0][i] = Math.exp(-x);  // spherical albedo
            spectralAlbedos[1][i] = Math.pow(spectralAlbedos[0][i], ak0);  // planar albedo
        }

        return new SpectralAlbedoResult(spectralAlbedos, 0, 0, 0, 0, 0, 0, 0, 0);
    }


    /**
     * Computes broadband albedos following AK algorithm given in 'Technical note_BBA_DECEMBER_2017.doc' (20171204).
     *
     * @param mu_0                 - mu0
     * @param d                    - snow grain diameter
     * @param refractiveIndexTable - table with refractive indices
     * @param solarSpectrumExtendedTable   - table with extended solar spectrum for sza = 0, 15, 30, 45, 60, 75 deg
     * @return double[]{pbbaVis, pbbaNir, pbbaSw};
     */
    static double[] computeBroadbandAlbedo(double mu_0, double d,
                                           RefractiveIndexTable refractiveIndexTable,
//                                           SolarSpectrumTable solarSpectrumTable,
                                           SolarSpectrumExtendedTable solarSpectrumExtendedTable,
                                           double sza) {

        if (Double.isNaN(d)) {
            return new double[]{Double.NaN, Double.NaN, Double.NaN};
        }

        double pbbaVis;
        double pbbaNir;
        double pbbaSw;
        double numeratorVis = 0.0;
        double denominatorVis = 0.0;
        double numeratorNir = 0.0;
        double denominatorNir = 0.0;
        double numeratorSw = 0.0;
        double denominatorSw = 0.0;

        double[] wvlsFull = solarSpectrumExtendedTable.getWvl();
//        final double[] f_lambda = SnowUtils.getFLambda(solarSpectrumTable);
        final double[] f_lambda = SnowUtils.getFLambda(solarSpectrumExtendedTable, sza);

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
            if (wvl > OlciSnowPropertiesConstants.BB_WVL_1 && wvl < OlciSnowPropertiesConstants.BB_WVL_2) {
                numeratorVis += planarSpectralAlbedo * f_lambda[i] * dx;
                denominatorVis += f_lambda[i] * dx;
            }
            // NIR 0.7-2.4
            if (wvl > OlciSnowPropertiesConstants.BB_WVL_2 && wvl < OlciSnowPropertiesConstants.BB_WVL_3) {
                numeratorNir += planarSpectralAlbedo * f_lambda[i] * dx;
                denominatorNir += f_lambda[i] * dx;
            }
            // SW 0.3-2.4
            if (wvl > OlciSnowPropertiesConstants.BB_WVL_1 && wvl < OlciSnowPropertiesConstants.BB_WVL_3) {
                numeratorSw += planarSpectralAlbedo * f_lambda[i] * dx;
                denominatorSw += f_lambda[i] * dx;
            }
        }

        pbbaVis = numeratorVis / denominatorVis;
        pbbaNir = numeratorNir / denominatorNir;
        pbbaSw = numeratorSw / denominatorSw;

        return new double[]{pbbaVis, pbbaNir, pbbaSw};
    }

    /**
     * Computes probability of photon absorption (PPA) at considered wavelengths.
     * Follows 'ppa_20171207.doc' (AK, 20171207)
     *
     * @param brr - Rayleigh reflectance
     * @param sza - sun zenith angle
     * @param vza - view zenith angle
     * @return double [] ppa
     */
    static double[] computeSpectralPPA(double[] brr, double sza, double vza) {
        double[] ppa = new double[brr.length];

        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);
        final double k_mu_0 = SnowUtils.computeU(mu_0);
        final double k_mu = SnowUtils.computeU(mu);
        final double brr_400 = Math.min(1.0, brr[0]);
        final double m = k_mu_0 * k_mu / brr_400;

        for (int i = 0; i < brr.length; i++) {
            final double y = Math.log(brr_400 / brr[i]) / m;
            ppa[i] = 3.0 * y * y / 64.0;
        }

        return ppa;
    }

    /**
     * Computes the snow grain diameter for given Rayleigh corrected reflectance at band 21 (1020nm).
     * Follows 'sgs_nov_20_865nm.doc' (AK, 20171120)
     *
     * @param refAlbedo - Spectral albedo at band 18 or 21 (865 or 1020nm)
     * @param refWvl    - Reference wavelength (865 or 1020nm)
     * @return the snow grain diameter in microns
     */
    static double computeGrainDiameter(double refAlbedo, double refWvl) {
        final double b = 3.62;
        final int numWvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length;
        final double chiRef = refWvl == 1020.0 ? OlciSnowPropertiesConstants.ICE_REFR_INDEX[numWvl - 1] :
                OlciSnowPropertiesConstants.ICE_REFR_INDEX[numWvl - 5];
        final double kappaRef = 4.0 * Math.PI * chiRef / (refWvl / 1000.0);  // eq. (4)
        return Math.log(refAlbedo) * Math.log(refAlbedo) / (b * b * kappaRef);
    }

    /**
     * Computes the snow grain diameter and soot concentration in case of polluted snow.
     * Algorithm Tech Note: 'Manual_31_01_2018.docx', AK 20180131
     *
     * @param brr400  - reflectance at 400nm
     * @param brr1020 - reflectance at 1200nm
     * @param sza     - SZA
     * @param vza     - VZA
     * @param raa     - rel. azimuth
     * @return double[] {grainDiam, soot} in mm and PPM
     */
    static double[] computePollutedSnowParams(double brr400, double brr1020,
                                              double sza, double vza, double raa) {

        // Wavelengths
        final double alam1 = 0.4;
        final double alam2 = 1.02;

        // Im(m)
        final double akas = 0.46;
        final double akai1 = 2.365e-11;
        final double akai2 = 2.25e-6;

        final double gi1 = 4. * Math.PI * akai1 / alam1;
        final double gi2 = 4. * Math.PI * akai2 / alam2;
        final double gs1 = 4. * Math.PI * akas / alam1;
        final double gs2 = 4. * Math.PI * akas / alam2;

        // Cosines
        final double raa1 = Math.abs(180. - raa);
        final double am0 = Math.cos(Math.PI * sza / 180.);
        final double am = Math.cos(Math.PI * vza / 180.);
        final double sam0 = Math.sin(Math.PI * sza / 180.);
        final double sam = Math.sin(Math.PI * vza / 180.);
        final double cosazi = Math.cos(Math.PI * raa1 / 180.);
        final double t = -am * am0 + sam0 * sam * cosazi;
        final double scat = Math.acos(t);
        final double theta = scat * 180. / Math.PI;

        // K(mu)
        final double ak0 = (3. / 7.) * (1. + 2. * am0);
        final double ak1 = (3. / 7.) * (1. + 2. * am);

        final double ax = 1.247;
        final double bx = 1.186;
        final double cx = 5.157;
        final double s1 = ax + bx * (am + am0) + cx * am * am0;
        final double s2 = 4. * (am + am0);
        final double px = 11.1 * Math.exp(-0.087 * theta) + 1.1 * Math.exp(-0.014 * theta);
        final double r0 = (s1 + px) / s2;

        // Calculations
        final double xx = r0 / ak0 / ak1;
        final double yy1 = Math.log(brr400 / r0);
        final double yy2 = Math.log(brr1020 / r0);
        final double y1 = xx * xx * yy1 * yy1;
        final double y2 = xx * xx * yy2 * yy2;
        final double v = y1 / y2;

        // Impurity content (soot)
        final double soot = (13. / 6.) * (gi2 * v - gi1) / (gs1 - gs2 * v);

        // Grain diameter
        final double grainDiam = y2 / (13. * gi2 + 6. * gs2 * soot);

        return new double[]{grainDiam / 1000., soot * 1.E6};  // in mm and PPM
    }

    static double computeR0PollutionThresh(double sza, double vza, double raa) {

        final double A = 1.247;
        final double B = 1.186;
        final double C = 5.157;

        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);
        final double cosTheta = SnowUtils.calcScatteringCos(sza, vza, raa);
        final double theta = Math.acos(cosTheta) * MathUtils.RTOD;
        final double p = 11.1 * Math.exp(-0.087 * theta) + 1.1 * Math.exp(-0.014 * theta);

        return (A + B * (mu_0 + mu) + C * mu_0 * mu + p) / (4.0 * (mu_0 + mu));
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
    private static void computeSpectralAlbedoWithSimpleApproximation(double[] brr,
                                                                     double sza, double vza, int numWvl,
                                                                     double refWvl,
                                                                     double[][] sphericalAlbedos) {
        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);
        final double k_mu_0 = SnowUtils.computeU(mu_0);
        final double k_mu = SnowUtils.computeU(mu);
        final double brr_400 = Math.min(1.0, brr[0]);
        final double xi = brr_400 / (k_mu_0 * k_mu);   // [1], eq. (5)

        final double brr_ref = refWvl == 1020.0 ? brr[3] : brr[2];
        final double r_ref = Math.pow(brr_ref / brr_400, xi); // [1], eq. (5)
        final double wvl_ref = refWvl / 1000.0;  // in microns!
        final double chi_ref = refWvl == 1020.0 ? OlciSnowPropertiesConstants.ICE_REFR_INDEX[numWvl - 1] :
                OlciSnowPropertiesConstants.ICE_REFR_INDEX[numWvl - 5];
        final double kappa_ref = 4.0 * Math.PI * chi_ref / wvl_ref;  // [2], eqs. (1), (2)
        final double B = Math.log(r_ref) * Math.log(r_ref) / kappa_ref;  // [2], eq. (2) transformed

        for (int i = 0; i < numWvl; i++) {
            final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];
            final double chi = OlciSnowPropertiesConstants.ICE_REFR_INDEX[i];
            final double kappa = 4.0 * Math.PI * chi / wvl;  // [2], eqs. (1), (2)

            sphericalAlbedos[0][i] = Math.exp(-Math.sqrt(B * kappa));  // spectral spherical abledo
        }
        sphericalAlbedos[1] = computePlanarFromSphericalAlbedos(sphericalAlbedos[0], sza);    // spectral planar abledo
    }

    /**
     * Computation of spectral spherical albedo following AK latest approach, 20180404.
     * Currently used as default.
     * <p>
     * References:
     * [1]: The determination of snow parameters using OLCI observations. TN AK, 20180404. File: Manual_04_04_2018-1.docx.
     */
    static void computeSpectralAlbedoFromTwoWavelengths(double[] brr,
                                                                double sza, double vza, int numWvl,
                                                                double[][] spectralAlbedos) {
        final double l = computeLFromTwoWavelengths(brr, sza, vza);
        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double k_mu_0 = SnowUtils.computeU(mu_0);

        for (int i = 0; i < numWvl; i++) {
            final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];
            final double chi = OlciSnowPropertiesConstants.ICE_REFR_INDEX[i];
            final double alpha_ice = 4.0 * Math.PI * chi / wvl;  // [2], eqs. (1), (2)

            spectralAlbedos[0][i] = Math.exp(-Math.sqrt(l * alpha_ice));  // spectral spherical abledo
            spectralAlbedos[1][i] = Math.exp(-k_mu_0 * Math.sqrt(l * alpha_ice));  // spectral planar abledo
        }
    }


    static void computeSpectralAlbedoFromTwoWavelengths_Oct2018(double[] brr,
                                                                double sza, double vza, int numWvl,
                                                                boolean polluted,
                                                                double[][] spectralAlbedos) {

        final double[] refWvl = new double[]{400., 560., 865., 1020.};
        final double[] akappa = new double[]{2.365e-11, 2.839e-9, 2.3877e-7, 2.25e-6};
        double[] alpha = new double[4];
        for (int i = 0; i < alpha.length; i++) {
            alpha[i] = 4.0*Math.PI*akappa[i]/refWvl[i];
        }

        final double consb = 0.3537;
        final double eps1 = 1./(1. - consb);
        final double eps2 = 1. - eps1;

        final double rr00 = Math.pow(brr[2], eps1) * Math.pow(brr[3], eps2);

        final double p1 = Math.log(brr[0]/rr00) * Math.log(brr[0]/rr00);
        final double p2 = Math.log(brr[1]/rr00) * Math.log(brr[1]/rr00);
        final double am = Math.log(p1/p2) / Math.log(refWvl[1]/refWvl[0]);

        final double amu1 = Math.cos(sza*MathUtils.DTOR);
        final double amu2 = Math.cos(vza*MathUtils.DTOR);

        final double u1 = SnowUtils.computeU(amu1);
        final double u2 = SnowUtils.computeU(amu2);

        final double x = u1 * u1 * u2 * u2 / (rr00*rr00);

        final double dlina = Math.log(brr[3]/rr00) * Math.log(brr[3]/rr00) / (x * x * alpha[3]);

        // dlina in mm:
        final double AL = 1.E-6 * dlina;
        final double aksi = 16.0 * 1.6 / 2.25;

        // diameter of grains (in mm):
        final double grainDiam = AL/aksi;

        // f (in 1/mm):
        final double SK = refWvl[0] / refWvl[3];
        final double f = p1 * Math.pow(SK, am) / (x * x * AL);


        for (int i = 0; i < numWvl; i++) {
            final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];
            final double chi = OlciSnowPropertiesConstants.ICE_REFR_INDEX[i];
            final double wvlNm = 1000. * wvl;
            final double alka = 4.0 * Math.PI * chi / wvlNm;
            if (polluted) {
                // sd9=f*(wsk/ws(4))**(-am)
                // TT(j)=alka(j)*1.e+6+sd9
                // arr(j)=rr00*exp(-x*sqrt(TT(j)*AL))
                // plane1= (arr(j)/rr00)**(rr00/u2)
                // spher1=(arr(j)/rr00)**(rr00/u2/u1)
                final double sd9 = f * Math.pow(wvlNm/refWvl[3], -am);
                final double tt = alka*1.E6 + sd9;
                final double arr = rr00 * Math.exp(-x*Math.sqrt(tt*AL));
                spectralAlbedos[0][i] = Math.pow(arr/rr00, rr00/(u1*u2));   // spectral spherical albedo
                spectralAlbedos[1][i] = Math.pow(arr/rr00, rr00/u2);        // spectral planar albedo
            } else {
                // TT(j)=alka(j)*1.e+6
                // arr(j)=rr00*exp(-x*sqrt(TT(j)*AL))
                // plane=exp(-u1*sqrt(TT(j)*AL))
                // spher=exp(-sqrt(TT(j)*AL))
                final double tt = alka*1.E6;
                spectralAlbedos[0][i] = Math.exp(-Math.sqrt(tt*AL));        // spectral spherical albedo
                spectralAlbedos[1][i] = Math.exp(-u1*Math.sqrt(tt*AL));     // spectral planar albedo
            }
        }
    }

    public static double computeLFromTwoWavelengths(double[] brr, double sza, double vza) {
        final int numWvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length;

        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);
        final double k_mu_0 = SnowUtils.computeU(mu_0);
        final double k_mu = SnowUtils.computeU(mu);

        final double chi_1 = OlciSnowPropertiesConstants.ICE_REFR_INDEX[numWvl - 5];
        final double chi_2 = OlciSnowPropertiesConstants.ICE_REFR_INDEX[numWvl - 1];
        final double alpha_1 = 4.0 * Math.PI * chi_1 / 0.865;
        final double alpha_2 = 4.0 * Math.PI * chi_2 / 1.02;
        final double b = Math.sqrt(alpha_1 / alpha_2);
        final double eps_1 = 1.0 / (1.0 - b);
        final double eps_2 = 1.0 / (1.0 - 1.0 / b);
        final double r_0 = Math.pow(brr[2], eps_1) * Math.pow(brr[3], eps_2);
        final double x = (k_mu_0 * k_mu) / r_0;   // [1], eq. (5)

        return Math.log(brr[3] / brr[2]) * Math.log(brr[3] / brr[2]) / (x * x * alpha_2);
    }

    static SpectralAlbedoResult computeSpectralAlbedosPollutedFromFourWavelengths(double[] brr,
                                                                                          double deltaBrr,
                                                                                          double sza,
                                                                                          double vza) {
        final int numWvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length;
        double[][] spectralAlbedos = new double[2][numWvl];

        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);
        final double k_mu_0 = SnowUtils.computeU(mu_0);
        final double k_mu = SnowUtils.computeU(mu);

        final double chi_3 = OlciSnowPropertiesConstants.ICE_REFR_INDEX[numWvl - 5];
        final double chi_4 = OlciSnowPropertiesConstants.ICE_REFR_INDEX[numWvl - 1];
//        final double alpha_3 = 4.0 * Math.PI * chi_3 / 0.865;
//        final double alpha_4 = 4.0 * Math.PI * chi_4 / 1.02;
        final double alpha_3 = 1.E-3 * 4.0 * Math.PI * chi_3 / 0.865;  // in mm!   (hint of AK, 20180705)
        final double alpha_4 = 1.E-3 * 4.0 * Math.PI * chi_4 / 1.02;
        final double b = Math.sqrt(alpha_3 / alpha_4);
        final double eps_1 = 1.0 / (1.0 - b);
        final double eps_2 = 1.0 / (1.0 - 1.0 / b);
        final double r0 = Math.pow(brr[2], eps_1) * Math.pow(brr[3], eps_2);
        final double x = (k_mu_0 * k_mu) / r0;   // [1], eq. (5)
//        final double l = Math.log(brr[3] / brr[0]) * Math.log(brr[3] / brr[0]) / (x * x * alpha_4);
        final double l = 1.E-6 * Math.log(brr[3] / r0) * Math.log(brr[3] / r0) / (x * x * alpha_4); // in mm!

        final double p_1 = Math.log(brr[0] / r0) * Math.log(brr[0] / r0);
        final double p_2 = Math.log(brr[1] / r0) * Math.log(brr[1] / r0);
        final double wvl_1 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[0];
        final double wvl_2 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[5];
//        final double m = Math.log(p_1 / p_2) / Math.log(wvl_1 / wvl_2);
        final double m = Math.log(p_1 / p_2) / Math.log(wvl_2 / wvl_1);  // wrong sign!  (hint of AK, 20180705)
//        final double f = p_1 * Math.pow(wvl_1, m) / (x * x * l);
        final double f = p_1 * Math.pow(wvl_1, m) / (x * x * l);

        for (int i = 0; i < numWvl; i++) {
            final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];
            final double chi = OlciSnowPropertiesConstants.ICE_REFR_INDEX[i];
//            final double alpha_ice = 4.0 * Math.PI * chi / wvl;  // [2], eqs. (1), (2)
            final double alpha_ice = 1.E3 * 4.0 * Math.PI * chi / wvl;  // [2], eqs. (1), (2) , wvl in mm !!

//            spectralAlbedos[0][i] = Math.exp(-Math.sqrt(alpha_ice + f * Math.pow(wvl, -m)) * l);  // spectral spherical albedo
//            spectralAlbedos[1][i] = Math.exp(-k_mu_0 * Math.sqrt(alpha_ice + f * Math.pow(wvl, -m)) * l);  // spectral planar albedo

            // l is under sqrt, found 20180706:
//            spectralAlbedos[0][i] = Math.exp(-Math.sqrt((alpha_ice + f * Math.pow(wvl, -m)) * l));  // spectral spherical albedo
//            spectralAlbedos[1][i] = Math.exp(-k_mu_0 * Math.sqrt((alpha_ice + f * Math.pow(wvl, -m)) * l));  // spectral planar albedo

            // AK code differs from tech note last eqs. on p.3:
//            final double beta = (alpha_ice + f * Math.pow(wvl, 1.0) * m) * l;
            final double beta = (alpha_ice + f * Math.pow(wvl, -m)) * l;     // switched again, AK 20180709
            spectralAlbedos[0][i] = Math.exp(-Math.sqrt(beta));  // spectral spherical albedo
            spectralAlbedos[1][i] = Math.exp(-k_mu_0 * Math.sqrt(beta));  // spectral planar albedo
        }

        // relative error estimation...
        double[] z = new double[4];
        for (int i = 0; i < z.length; i++) {
            z[i] = 1.0 / Math.log(brr[i] / r0);
        }

        final double s = 0.0; // todo: unreadable, check with AK, set to 0.0 in the meantime

        final double nu_1 = 2.0*eps_1 - 2.0 * eps_1 * z[3];
        final double nu_2 = 2.0*eps_2 + 2.0*eps_1*1.0/Math.log(brr[3] / r0);

        double[] w = new double[4];
        w[0] = 2.0*z[0] / (m*Math.log(wvl_2/wvl_1));
        w[1] = -2.0*z[1] / (m*Math.log(wvl_2/wvl_1));
        w[2] = -(w[0] + w[1])*eps_1;
        w[3] = -(w[0] + w[1])*eps_2;

        double[] h = new double[4];
        h[0] = 2.0*(1.0 + s) *z[0];
        h[1] = -2.0*s*z[1];
        h[2] = -eps_1*(h[0] + h[1] + 2.0*z[3]);
        h[3] = 2.0*eps_2*(z[3] + s*z[1] - (1.0 + s*z[0])*z[0]) - 2.0*z[3];

        final double r0RelErr = computeR0RelErr(r0, brr, eps_1, eps_2, deltaBrr);
        final double lRelErr = computeLRelErr(l, brr, nu_1, nu_2, deltaBrr);
        final double mRelErr = computeMRelErr(m, brr, w, deltaBrr);
        final double fRelErr = computeFRelErr(f, brr, h, deltaBrr);

        return new SpectralAlbedoResult(spectralAlbedos, r0, f, l, m, r0RelErr, fRelErr, lRelErr, mRelErr);
    }

    /**
     * Computes planar from spherical albedos at considered wavelengths.
     * Follows 'nov_20.doc' (AK, 20171120)
     *
     * @param sphericaAlbedo - the spherical albedo value
     * @param sza            - sun zenith angle
     * @return theplanar albedo value
     */
    private static double computePlanarFromSphericalAlbedo(double sphericaAlbedo,
                                                           double sza) {
        // eq. (1):
        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        return Math.pow(sphericaAlbedo, SnowUtils.computeU(mu_0));
    }

    private static double[] computePlanarFromSphericalAlbedos(double[] sphericalAlbedos, double sza) {
        double[] planarAlbedos = new double[sphericalAlbedos.length];
        for (int i = 0; i < planarAlbedos.length; i++) {
            planarAlbedos[i] = computePlanarFromSphericalAlbedo(sphericalAlbedos[i], sza);
        }
        return planarAlbedos;
    }

    private static double computeR0RelErr(double r0, double[] brr, double eps_1, double eps_2, double deltaBrr) {
        // AK: 'technical_note_JUNE_20_2018.docx', eq. (4)
        double r0RelErr = Math.sqrt(eps_1 * eps_1 * deltaBrr * deltaBrr / (brr[2] * brr[2]) +
                                            eps_2 * eps_2 * deltaBrr * deltaBrr / (brr[3] * brr[3]));
        // make sure error is positive and does not exceed value itself
        return Math.min(Math.abs(r0), Math.abs(r0RelErr));
    }

    private static double computeLRelErr(double l, double[] brr, double nu_1, double nu_2, double deltaBrr) {
        // AK: 'technical_note_JUNE_20_2018.docx', eq. (4)
        double lRelErr = Math.sqrt(nu_1 * nu_1 * deltaBrr * deltaBrr / (brr[2] * brr[2]) +
                                           nu_2 * nu_2 * deltaBrr * deltaBrr / (brr[3] * brr[3]));
        // make sure error is positive and does not exceed value itself
        return Math.min(Math.abs(l), Math.abs(lRelErr));
    }

    private static double computeMRelErr(double m, double[] brr, double[] w, double deltaBrr) {
        // AK: 'technical_note_JUNE_20_2018.docx', eq. (4)
        double sum = 0.0;
        for (int i = 0; i < brr.length; i++) {
            sum += w[i]*w[i]*deltaBrr*deltaBrr/(brr[i]*brr[i]);
        }
        final double mRelErr = Math.sqrt(sum);
        // make sure error is positive and does not exceed value itself
        return Math.min(Math.abs(m), Math.abs(mRelErr));
    }

    private static double computeFRelErr(double f, double[] brr, double[] h, double deltaBrr) {
        // AK: 'technical_note_JUNE_20_2018.docx', eq. (4)
        double sum = 0.0;
        for (int i = 0; i < brr.length; i++) {
            sum += h[i]*h[i]*deltaBrr*deltaBrr/(brr[i]*brr[i]);
        }
        final double fRelErr = Math.sqrt(sum);
        return Math.min(Math.abs(f), Math.abs(fRelErr));
    }

    
    static class SpectralAlbedoResult {
        private final double[][] spectralAlbedos;
        private final double r0;
        private final double f;
        private final double l;
        private final double m;

        private final double r0RelErr;
        private final double fRelErr;
        private final double lRelErr;
        private final double mRelErr;

        SpectralAlbedoResult(double[][] spectralAlbedos, double r0, double f, double l, double m,
                             double r0RelErr, double fRelErr, double lRelErr, double mRelErr) {
            this.spectralAlbedos = spectralAlbedos;
            this.r0 = r0;
            this.f = f;
            this.l = l;
            this.m = m;
            this.r0RelErr = r0RelErr;
            this.fRelErr = fRelErr;
            this.lRelErr = lRelErr;
            this.mRelErr = mRelErr;
        }

        public double[][] getSpectralAlbedos() {
            return spectralAlbedos;
        }

        public double getR0() {
            return r0;
        }

        public double getF() {
            return f;
        }

        public double getL() {
            return l;
        }

        public double getM() {
            return m;
        }

        public double getR0RelErr() {
            return r0RelErr;
        }

        public double getfRelErr() {
            return fRelErr;
        }

        public double getlRelErr() {
            return lRelErr;
        }

        public double getmRelErr() {
            return mRelErr;
        }
    }
}

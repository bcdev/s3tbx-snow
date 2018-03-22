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
     *
     * @return double[][] spectralAlbedos : spherical and planar at OLCI wavelengths
     */
    static double[][] computeSpectralAlbedos(double[] brr, double sza, double vza,
                                             double referenceWavelength,
                                             SpectralAlbedoMode mode) {
        final int numWvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length;
        double[][] spectralAlbedos = new double[2][numWvl];

        if (mode == SpectralAlbedoMode.SIMPLE_APPROXIMATION) {
            computeSpectralSphericalAlbedoWithSimpleApproximation(brr, sza, vza, numWvl,
                                                                  referenceWavelength, spectralAlbedos);
        } else {
            // we no longer support this
            throw new IllegalArgumentException("spectral albedo algoritm mode " + mode.getName() + " not supported");
        }
        spectralAlbedos[1] = computePlanarFromSphericalAlbedos(spectralAlbedos[0], sza);

        return spectralAlbedos;
    }

    /**
     * Computes spectral spherical and planar albedos in case of polluted snow.
     * Algorithm Tech Note: 'Manual_31_01_2018.docx', AK 20180131
     *
     * @param pollutedSnowParams - double[]{grainDiam, soot}, the result from {@code computePollutedSnowParams}.
     * @param sza - SZA
     *
     * @return double[][] spectralAlbedos : spherical and planar at OLCI wavelengths
     *
     * @see OlciSnowAlbedoAlgorithm#computePollutedSnowParams(double, double, double, double, double)
     */
    static double[][] computeSpectralAlbedosPolluted(double[] pollutedSnowParams, double sza) {
        final int numWvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length;
        double[][] spectralAlbedos = new double[2][numWvl];

        final double grainDiam = pollutedSnowParams[0] * 1000.0;  // in microns here
        final double soot = pollutedSnowParams[1] * 1.E-6; // dimensionless here

        final double akas = 0.46;
        final double am0=Math.cos(Math.PI*sza/180.);
        final double ak0=(3./7.) *(1.+2.*am0);

        for (int i = 0; i < numWvl; i++) {
            final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[i];
            final double p = 6.*4.*Math.PI*akas*soot*grainDiam/wvl;
            final double chi = OlciSnowAlbedoConstants.ICE_REFR_INDEX[i];   // chi = ak(i) in 'soot.f' breadboard
            final double x = Math.sqrt (13.*4.*Math.PI*grainDiam*chi/wvl + p);
            spectralAlbedos[0][i] = Math.exp(-x);  // spherical albedo
            spectralAlbedos[1][i] = Math.pow(spectralAlbedos[0][i], ak0);  // planar albedo
        }

        return spectralAlbedos;
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

    /**
     * Computes probability of photon absorption (PPA) at considered wavelengths.
     * Follows 'ppa_20171207.doc' (AK, 20171207)
     *
     * @param brr - Rayleigh reflectance
     * @param sza - sun zenith angle
     * @param vza - view zenith angle
     *
     * @return double [] ppa
     */
    static double[] computeSpectralPPA(double[] brr, double sza, double vza) {
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
     * Computes the snow grain diameter and soot concentration in case of polluted snow.
     * Algorithm Tech Note: 'Manual_31_01_2018.docx', AK 20180131
     *
     * @param brr400 - reflectance at 400nm
     * @param brr1020 - reflectance at 1200nm
     * @param sza - SZA
     * @param vza - VZA
     * @param raa - rel. azimuth
     *
     * @return double[] {grainDiam, soot} in mm and PPM
     */
    static double[] computePollutedSnowParams(double brr400, double brr1020,
                                              double sza, double vza, double raa) {

        // Wavelengths
        final double alam1=0.4;
        final double alam2=1.02;

        // Im(m)
        final double akas=0.46;
        final double akai1=2.365e-11;
        final double akai2=2.25e-6;

        final double gi1=4.*Math.PI*akai1/alam1;
        final double gi2=4.*Math.PI*akai2/alam2;
        final double gs1=4.*Math.PI*akas/alam1;
        final double gs2=4.*Math.PI*akas/alam2;

        // Cosines
        final double raa1=Math.abs(180.-raa);
        final double am0=Math.cos(Math.PI*sza/180.);
        final double am=Math.cos(Math.PI*vza/180.);
        final double sam0=Math.sin(Math.PI*sza/180.);
        final double sam=Math.sin(Math.PI*vza/180.);
        final double cosazi=Math.cos(Math.PI*raa1/180.);
        final double t=-am*am0+sam0*sam*cosazi;
        final double scat=Math.acos(t);
        final double theta=scat*180./Math.PI;

        // K(mu)
        final double ak0 = (3./7.) *(1.+2.*am0);
        final double ak1 = (3./7.) *(1.+2.*am);

        final double ax = 1.247;
        final double bx = 1.186;
        final double cx = 5.157;
        final double s1 = ax+bx*(am+am0)+cx*am*am0;
        final double s2 = 4.*(am+am0);
        final double px = 11.1*Math.exp(-0.087*theta)+1.1*Math.exp(-0.014*theta);
        final double r0 = (s1+px)/s2;

        // Calculations
        final double xx = r0/ak0/ak1;
        final double yy1 = Math.log(brr400/r0);
        final double yy2 = Math.log(brr1020/r0);
        final double y1 = xx*xx * yy1*yy1;
        final double y2 = xx*xx * yy2*yy2;
        final double v = y1/y2;

        // Impurity content (soot)
        final double soot =(13./6.)*(gi2*v-gi1)/(gs1-gs2*v);

        // Grain diameter
        final double grainDiam = y2/(13.*gi2+6.*gs2*soot);

        return new double[]{grainDiam/1000., soot*1.E6};  // in mm and PPM
    }

    static double computeR0PollutionThresh(double sza, double vza, double raa) {

        final double A = 1.247;
        final double B = 1.186;
        final double C = 5.157;

        final double mu_0 = Math.cos(sza * MathUtils.DTOR);
        final double mu = Math.cos(vza * MathUtils.DTOR);
        final double cosTheta = SnowUtils.calcScatteringCos(sza, vza, raa);
        final double theta = Math.acos(cosTheta) * MathUtils.RTOD;
        final double p = 11.1*Math.exp(-0.087*theta) + 1.1*Math.exp(-0.014*theta);

        return (A + B*(mu_0 + mu) + C*mu_0*mu + p)/(4.0*(mu_0 + mu));
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

}

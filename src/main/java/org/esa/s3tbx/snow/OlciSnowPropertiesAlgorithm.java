package org.esa.s3tbx.snow;

import org.esa.s3tbx.snow.math.Integrator;
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
     * Computes spectral spherical and planar albedos for CLEAN snow using AK latest algo from Oct 2018
     *
     * @param brr      - subset of BRR spectrum  (BRR_01, BRR_05, BRR_21, BRR_17)
     * @param deltaBrr - assumed uncertainty in Rayleigh corrected reflectances
     * @param sza      - sun zenith angle (deg)
     * @param vza      - view zenith angle (deg)
     * @return SpectralAlbedoResult
     */
    static SpectralAlbedoResult computeSpectralAlbedos(double[] brr, double deltaBrr, double sza, double vza) {
        return computeSpectralAlbedoFromTwoWavelengths(brr, deltaBrr, sza, vza, false);
    }

    /**
     * Computes spectral spherical and planar albedos in case of polluted snow.
     * Algorithm Tech Note: 'Manual_31_01_2018.docx', AK 20180131
     *
     * @param brr - subset of BRR spectrum  (BRR_01, BRR_06, BRR_21, BRR_17)
     * @param sza - sun zenith angle (deg)
     * @param vza - view zenith angle (deg)
     * @return double[][] spectralAlbedos : spherical and planar at OLCI wavelengths
     */
    static SpectralAlbedoResult computeSpectralAlbedosPolluted(double[] brr,
                                                               double sza,
                                                               double vza,
                                                               double deltaBrr) {
//        final int numWvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length;
//        return computeSpectralAlbedosPollutedFromFourWavelengths(brr, deltaBrr, sza, vza);
        // now using AK latest algo from Oct 2018:
        return computeSpectralAlbedoFromTwoWavelengths(brr, deltaBrr, sza, vza, true);
    }

    /**
     * Computes broadband albedos following AK algorithm given in
     * - 'TN_spectral_albedo_and_grain_size.docx' (20181015) for planar spectral albedo retrieval
     * - 'Technical note_BBA_DECEMBER_2017.doc' (20171204) for BB integration
     * <p>
     * See also KOKHANOVSKY et al. (2018) paper in The Cryosphere
     *
     * @param brr        - subset of BRR spectrum  (BRR_01, BRR_06, BRR_21, BRR_17)
     * @param deltaBrr   - assumed BRR uncertainty
     * @param isPolluted - indicator for clean or polluted snow
     * @param sza        - solar zenith angle
     * @param vza        - view zenith angle
     * @return double[] broadband albedos in 3 wavelength regions VIS, NIR, SW (300-700nm, 700-2400nm, 300-2400nm)
     */
    static SpectralAlbedoResult computeSpectralAlbedoFromTwoWavelengths(double[] brr, double deltaBrr,
                                                                        double sza, double vza,
                                                                        boolean isPolluted) {

        final int numWvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length;
        double[][] spectralAlbedos = new double[2][numWvl];

        final double[] refWvl = new double[]{400., 560., 865., 1020.};
        final double[] akappa = new double[]{2.365e-11, 2.839e-9, 2.3877e-7, 2.25e-6};
        double[] alpha = new double[4];
        for (int i = 0; i < alpha.length; i++) {
            alpha[i] = 4.0 * Math.PI * akappa[i] / refWvl[i];
        }

        final double consb = 0.3537;
        final double eps_1 = 1. / (1. - consb);
        final double eps_2 = 1. - eps_1;

        final double r0 = Math.pow(brr[2], eps_1) * Math.pow(brr[3], eps_2);

        final double p1 = Math.log(brr[0] / r0) * Math.log(brr[0] / r0);
        final double p2 = Math.log(brr[1] / r0) * Math.log(brr[1] / r0);
        final double m = Math.log(p1 / p2) / Math.log(refWvl[1] / refWvl[0]);

        final double amu1 = Math.cos(sza * MathUtils.DTOR);
        final double amu2 = Math.cos(vza * MathUtils.DTOR);

        final double u1 = SnowUtils.computeU(amu1);
        final double u2 = SnowUtils.computeU(amu2);

        final double x = u1 * u2 / r0;

        final double dlina = Math.log(brr[3] / r0) * Math.log(brr[3] / r0) / (x * x * alpha[3]);

        // dlina in mm:
        final double l = 1.E-6 * dlina;
        // diameter of grains (in mm):
        // final double aksi = 16.0 * 1.6 / 2.25;     // --> 11.38
        // final double grainDiam = l / aksi;         // not needed here

        // f (in 1/mm):
        final double SK = refWvl[0] / refWvl[3];
        final double f = p1 * Math.pow(SK, m) / (x * x * l);


        for (int i = 0; i < numWvl; i++) {
            final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];
            final double chi = OlciSnowPropertiesConstants.ICE_REFR_INDEX[i];
            final double wvlNm = 1000. * wvl;
            final double alka = 4.0 * Math.PI * chi / wvlNm;
            if (isPolluted) {
                // sd9=f*(wsk/ws(4))**(-m)
                // TT(j)=alka(j)*1.e+6+sd9
                // arr(j)=r0*exp(-x*sqrt(TT(j)*l))
                // plane1= (arr(j)/r0)**(r0/u2)
                // spher1=(arr(j)/r0)**(r0/u2/u1)
                final double sd9 = f * Math.pow(wvlNm / refWvl[3], -m);
                final double tt = alka * 1.E6 + sd9;
                final double arr = r0 * Math.exp(-x * Math.sqrt(tt * l));
                spectralAlbedos[0][i] = Math.pow(arr / r0, r0 / (u1 * u2));   // spectral spherical albedo
                spectralAlbedos[1][i] = Math.pow(arr / r0, r0 / u2);        // spectral planar albedo
            } else {
                // TT(j)=alka(j)*1.e+6
                // arr(j)=r0*exp(-x*sqrt(TT(j)*l))
                // plane=exp(-u1*sqrt(TT(j)*l))
                // spher=exp(-sqrt(TT(j)*l))
                final double tt = alka * 1.E6;
                spectralAlbedos[0][i] = Math.exp(-Math.sqrt(tt * l));        // spectral spherical albedo
                spectralAlbedos[1][i] = Math.exp(-u1 * Math.sqrt(tt * l));     // spectral planar albedo
            }
        }

        // relative error estimation for isPolluted snow following AK...
        if (isPolluted) {
            final double wvl_1 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[0];
            final double wvl_2 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[5];
            double[] z = new double[4];
            for (int i = 0; i < z.length; i++) {
                z[i] = 1.0 / Math.log(brr[i] / r0);
            }

            final double s = 0.0; // todo: unreadable, check with AK, set to 0.0 in the meantime

            final double nu_1 = 2.0 * eps_1 - 2.0 * eps_1 * z[3];
            final double nu_2 = 2.0 * eps_2 + 2.0 * eps_1 * 1.0 / Math.log(brr[3] / r0);

            double[] w = new double[4];
            w[0] = 2.0 * z[0] / (m * Math.log(wvl_2 / wvl_1));
            w[1] = -2.0 * z[1] / (m * Math.log(wvl_2 / wvl_1));
            w[2] = -(w[0] + w[1]) * eps_1;
            w[3] = -(w[0] + w[1]) * eps_2;

            double[] h = new double[4];
            h[0] = 2.0 * (1.0 + s) * z[0];
            h[1] = -2.0 * s * z[1];
            h[2] = -eps_1 * (h[0] + h[1] + 2.0 * z[3]);
            h[3] = 2.0 * eps_2 * (z[3] + s * z[1] - (1.0 + s * z[0]) * z[0]) - 2.0 * z[3];

            final double r0RelErr = computeR0RelErr(r0, brr, deltaBrr);
            final double lRelErr = computeLRelErr(l, brr, nu_1, nu_2, deltaBrr);
            final double mRelErr = computeMRelErr(m, brr, w, deltaBrr);
            final double fRelErr = computeFRelErr(f, brr, h, deltaBrr);

            return new SpectralAlbedoResult(spectralAlbedos, r0, f, l, m, r0RelErr, fRelErr, lRelErr, mRelErr);
        } else {
            return new SpectralAlbedoResult(spectralAlbedos, r0, f, l, m, 0.0, 0.0, 0.0, 0.0);
        }
    }

    /**
     * Computes broadband albedos following AK algorithm given in
     * - 'TN_spectral_albedo_and_grain_size.docx' (20181015) for planar spectral albedo retrieval
     * - 'Technical note_BBA_DECEMBER_2017.doc' (20171204) for BB integration
     * --> Applies now Simpson's rule for BB integration (request AK Oct 2018).
     * --> this code is used in version 2.0.11-SNAPSHOT, with table 'final_table_fluxes_angles.txt'
     *
     * @param mu_0                       - mu0
     * @param brr                        - subset of BRR spectrum  (BRR_01, BRR_06, BRR_21, BRR_17)
     * @param isPolluted                 - indicator for clean or polluted snow
     * @param refractiveIndexTable       - table with refractive indices
     * @param solarSpectrumExtendedTable - table with extended solar spectrum for sza = 0, ..., 88 deg
     * @param sza                        - solar zenith angle
     * @param vza                        - view zenith angle
     * @return double[] broadband albedos in 3 wavelength regions VIS, NIR, SW (300-700nm, 700-2400nm, 300-2400nm)
     */
//    static double[] computeBroadbandAlbedo_old(double mu_0,
//                                               double[] brr,
//                                               boolean isPolluted,
//                                               RefractiveIndexTable refractiveIndexTable,
//                                               SolarSpectrumExtendedTable solarSpectrumExtendedTable,
//                                               double sza, double vza) {
//
//        double bbaVis;
//        double bbaNir;
//        double bbaSw;
//        double numeratorVis;
//        double denominatorVis;
//        double numeratorNir;
//        double denominatorNir;
//        double numeratorSw;
//        double denominatorSw;
//
//        double[] wvlsFull = solarSpectrumExtendedTable.getWvl();
//        final double[] fLambda = SnowUtils.computeFLambda(solarSpectrumExtendedTable, sza);
//
//        double[] planarSpectralAlbedo = computeFullPlanarSpectralAlbedo(mu_0, brr,
//                                                                        refractiveIndexTable,
//                                                                        solarSpectrumExtendedTable, vza,
//                                                                        isPolluted);
//        double[] fLambdaTimesPlanarSpectralAlbedo = new double[wvlsFull.length];
//        for (int i = 0; i < fLambdaTimesPlanarSpectralAlbedo.length; i++) {
//            fLambdaTimesPlanarSpectralAlbedo[i] = planarSpectralAlbedo[i] * fLambda[i];
//        }
//
//        numeratorVis = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
//                                                   OlciSnowPropertiesConstants.BB_WVL_2, fLambdaTimesPlanarSpectralAlbedo, wvlsFull);
//        denominatorVis = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
//                                                     OlciSnowPropertiesConstants.BB_WVL_2, fLambda, wvlsFull);
//
//        numeratorNir = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_2,
//                                                   OlciSnowPropertiesConstants.BB_WVL_3, fLambdaTimesPlanarSpectralAlbedo, wvlsFull);
//        denominatorNir = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_2,
//                                                     OlciSnowPropertiesConstants.BB_WVL_3, fLambda, wvlsFull);
//
//        numeratorSw = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
//                                                  OlciSnowPropertiesConstants.BB_WVL_3, fLambdaTimesPlanarSpectralAlbedo, wvlsFull);
//        denominatorSw = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
//                                                    OlciSnowPropertiesConstants.BB_WVL_3, fLambda, wvlsFull);
//
//        bbaVis = numeratorVis / denominatorVis;
//        bbaNir = numeratorNir / denominatorNir;
//        bbaSw = numeratorSw / denominatorSw;
//
//        return new double[]{bbaVis, bbaNir, bbaSw};
//    }

    static double[] computeBroadbandAlbedo(double mu_0,
                                           double[] brr,
                                           boolean isPolluted,
                                           RefractiveIndexTable refractiveIndexTable,
                                           SolarSpectrumExtendedTable solarSpectrumExtendedTable,
                                           double sza, double vza) {

        double bbaVis;
        double bbaNir;
        double bbaSw;
        double numeratorVis;
        double denominatorVis;
        double numeratorNir;
        double denominatorNir;
        double numeratorSw;
        double denominatorSw;

        double[] wvlsFull = solarSpectrumExtendedTable.getWvl();
        final double[] fLambda = SnowUtils.computeFLambda(solarSpectrumExtendedTable, sza);

        double[] planarSpectralAlbedo = computeFullPlanarSpectralAlbedo(mu_0, brr,
                                                                        refractiveIndexTable,
                                                                        solarSpectrumExtendedTable, vza,
                                                                        isPolluted);
        double[] planarCleanSpectralAlbedo = null;
        if (isPolluted) {
            planarSpectralAlbedo = computeFullPlanarSpectralAlbedoPolluted(mu_0, brr,
                                                                                   refractiveIndexTable,
                                                                                   solarSpectrumExtendedTable, vza);
            planarCleanSpectralAlbedo = computeFullPlanarSpectralAlbedo(mu_0, brr,
                                                                        refractiveIndexTable,
                                                                        solarSpectrumExtendedTable, vza,
                                                                        false);
        }
        double[] fLambdaTimesPlanarSpectralAlbedo = new double[wvlsFull.length];
        final double wvlMax =
                OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[OlciSnowPropertiesConstants.OLCI_NUM_WVLS-1];
        for (int i = 0; i < fLambdaTimesPlanarSpectralAlbedo.length; i++) {
            if (isPolluted && wvlsFull[i] > wvlMax) {
                fLambdaTimesPlanarSpectralAlbedo[i] = planarCleanSpectralAlbedo[i] * fLambda[i];
            } else {
                fLambdaTimesPlanarSpectralAlbedo[i] = planarSpectralAlbedo[i] * fLambda[i];
            }
        }


        numeratorVis = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
                                                   OlciSnowPropertiesConstants.BB_WVL_2,
                                                   fLambdaTimesPlanarSpectralAlbedo,
                                                   wvlsFull);
        denominatorVis = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
                                                     OlciSnowPropertiesConstants.BB_WVL_2,
                                                     fLambda,
                                                     wvlsFull);
        numeratorSw = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
                                                  OlciSnowPropertiesConstants.BB_WVL_3,
                                                  fLambdaTimesPlanarSpectralAlbedo,
                                                  wvlsFull);
        denominatorSw = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
                                                    OlciSnowPropertiesConstants.BB_WVL_3,
                                                    fLambda,
                                                    wvlsFull);
        numeratorNir = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_2,
                                                   OlciSnowPropertiesConstants.BB_WVL_3,
                                                   fLambdaTimesPlanarSpectralAlbedo,
                                                   wvlsFull);
        denominatorNir = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_2,
                                                     OlciSnowPropertiesConstants.BB_WVL_3,
                                                     fLambda,
                                                     wvlsFull);

        bbaVis = numeratorVis / denominatorVis;
        bbaNir = numeratorNir / denominatorNir;
        bbaSw = numeratorSw / denominatorSw;

        return new double[]{bbaVis, bbaNir, bbaSw};
    }

    /**
     * Computes probability of photon absorption (PPA) at considered wavelengths.
     * Follows 'ppa_new.doc' (AK, 20181031)
     *
     * @param d - grain diameter in mm
     * @param f - f parameter for polluted snow
     * @param m - m parameter for polluted snow
     * @return double [] ppa
     */
    static double[] computeSpectralPPA_oct2018(double d, double f, double m) {
        double[] ppa = new double[OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length];

        final double lambda0 = 1.0;  // microns
        final double dMicrons = d * 1000.0;  // mm --> microns
        final double fPerMicrons = f / 1000.0;   // 1/mm --> 1/microns
        for (int i = 0; i < OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length; i++) {
            final double lambda = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];   // microns
            final double chi = OlciSnowPropertiesConstants.ICE_REFR_INDEX[i];
            final double alpha = 4.0 * Math.PI * chi / lambda;
            ppa[i] = 1.6 * dMicrons * (alpha + fPerMicrons * Math.pow(lambda / lambda0, -m));
        }

        return ppa;
    }

    /**
     * Computes threshold of R0 reflectance to identify polluted snow
     *
     * @param sza - sun zenith angle
     * @param vza - view zenith angle
     * @param raa - relative azimuth angle
     * @return threshold of R0 reflectance
     */
    static double computeR0ReflectancePollutionThresh(double sza, double vza, double raa) {

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

    //////////////////////////////// end of public //////////////////////////////////////////////////////////////

    static double[] computeFullPlanarSpectralAlbedo(double mu_0, double[] brr,
                                                    RefractiveIndexTable refractiveIndexTable,
                                                    SolarSpectrumExtendedTable solarSpectrumExtendedTable,
                                                    double vza,
                                                    boolean polluted) {

        // Same as computeSpectralAlbedoFromTwoWavelengths, but for computation of only planarSpectralAlbedo,
        // but here on fine 5nm grid rather than OLCI wavelengths only.
        final double[] refWvl = new double[]{400., 560., 865., 1020.};
        final double[] akappa = new double[]{2.365e-11, 2.839e-9, 2.3877e-7, 2.25e-6};
        double[] alpha = new double[4];
        for (int i = 0; i < alpha.length; i++) {
            alpha[i] = 4.0 * Math.PI * akappa[i] / refWvl[i];
        }

        final double consb = 0.3537;
        final double eps_1 = 1. / (1. - consb);
        final double eps_2 = 1. - eps_1;

        final double r0 = Math.pow(brr[2], eps_1) * Math.pow(brr[3], eps_2);

        final double p1 = Math.log(brr[0] / r0) * Math.log(brr[0] / r0);
        final double p2 = Math.log(brr[1] / r0) * Math.log(brr[1] / r0);
        final double m = Math.log(p1 / p2) / Math.log(refWvl[1] / refWvl[0]);

        final double mu_1 = Math.cos(vza * MathUtils.DTOR);

        final double u1 = SnowUtils.computeU(mu_0);
        final double u2 = SnowUtils.computeU(mu_1);

        final double x = u1 * u2  / r0;

        final double dlina = Math.log(brr[3] / r0) * Math.log(brr[3] / r0) / (x * x * alpha[3]);

        // dlina in mm:
        final double l = 1.E-6 * dlina;

        // f (in 1/mm):
        final double SK = refWvl[0] / refWvl[3];
        final double f = p1 * Math.pow(SK, m) / (x * x * l);

        double[] wvlsFull = solarSpectrumExtendedTable.getWvl();
        final int numWvl = wvlsFull.length;

        double[] planarSpectralAlbedo = new double[numWvl];

        for (int i = 0; i < numWvl; i++) {
            final double wvl = wvlsFull[i];
            final double chi = refractiveIndexTable.getRefractiveIndexImag(i);
            final double wvlNm = 1000. * wvl;
            final double alka = 4.0 * Math.PI * chi / wvlNm;
            if (polluted) {
                // FORTRAN:
                // sd9=f*(wsk/ws(4))**(-m)
                // TT(j)=alka(j)*1.e+6+sd9
                // arr(j)=r0*exp(-x*sqrt(TT(j)*l))
                // plane1= (arr(j)/r0)**(r0/u2)
                // spher1=(arr(j)/r0)**(r0/u2/u1)
                final double sd9 = f * Math.pow(wvlNm / refWvl[3], -m);
                final double tt = alka * 1.E6 + sd9;
                final double arr = r0 * Math.exp(-x * Math.sqrt(tt * l));
                planarSpectralAlbedo[i] = Math.pow(arr / r0, r0 / u2);
            } else {
                // FORTRAN:
                // TT(j)=alka(j)*1.e+6
                // arr(j)=r0*exp(-x*sqrt(TT(j)*l))
                // plane=exp(-u1*sqrt(TT(j)*l))
                // spher=exp(-sqrt(TT(j)*l))
                final double tt = alka * 1.E6;
                planarSpectralAlbedo[i] = Math.exp(-u1 * Math.sqrt(tt * l));
            }
        }
        return planarSpectralAlbedo;
    }

    static double[] computeFullPlanarSpectralAlbedoPolluted(double mu_0, double[] brr,
                                                            RefractiveIndexTable refractiveIndexTable,
                                                            SolarSpectrumExtendedTable solarSpectrumExtendedTable,
                                                            double vza) {

        // Same as computeSpectralAlbedoFromTwoWavelengths, but for computation of only planarSpectralAlbedo,
        // but here on fine 5nm grid rather than OLCI wavelengths only.
        final double[] refWvl = new double[]{400., 560., 865., 1020.};
        final double[] akappa = new double[]{2.365e-11, 2.839e-9, 2.3877e-7, 2.25e-6};
        double[] alpha = new double[4];
        for (int i = 0; i < alpha.length; i++) {
            alpha[i] = 4.0 * Math.PI * akappa[i] / refWvl[i];
        }

        final double consb = 0.3537;
        final double eps_1 = 1. / (1. - consb);
        final double eps_2 = 1. - eps_1;

        final double r0 = Math.pow(brr[2], eps_1) * Math.pow(brr[3], eps_2);

        final double p1 = Math.log(brr[0] / r0) * Math.log(brr[0] / r0);
        final double p2 = Math.log(brr[1] / r0) * Math.log(brr[1] / r0);
        final double m = Math.log(p1 / p2) / Math.log(refWvl[1] / refWvl[0]);

        final double mu_1 = Math.cos(vza * MathUtils.DTOR);

        final double u1 = SnowUtils.computeU(mu_0);
        final double u2 = SnowUtils.computeU(mu_1);

        final double x = u1 * u2 / r0;

        final double dlina = Math.log(brr[3] / r0) * Math.log(brr[3] / r0) / (x * x * alpha[3]);

        // dlina in mm:
        final double l = 1.E-6 * dlina;

        // f (in 1/mm):
        final double SK = refWvl[0] / refWvl[3];
        final double f = p1 * Math.pow(SK, m) / (x * x * l);

        double[] wvlsFull = solarSpectrumExtendedTable.getWvl();
        final int numWvl = wvlsFull.length;

        double[] planarSpectralAlbedo = new double[numWvl];

        for (int i = 10; i < numWvl; i++) {
            // this is for wvl >= 400nm
            final double wvl = wvlsFull[i];
            final double chi = refractiveIndexTable.getRefractiveIndexImag(i);
            final double wvlNm = 1000. * wvl;
            final double alka = 4.0 * Math.PI * chi / wvlNm;
            // FORTRAN:
            // sd9=f*(wsk/ws(4))**(-m)
            // TT(j)=alka(j)*1.e+6+sd9
            // arr(j)=r0*exp(-x*sqrt(TT(j)*l))
            // plane1= (arr(j)/r0)**(r0/u2)
            // spher1=(arr(j)/r0)**(r0/u2/u1)
            final double sd9 = f * Math.pow(wvlNm / refWvl[3], -m);
            final double tt = alka * 1.E6 + sd9;
            final double arr = r0 * Math.exp(-x * Math.sqrt(tt * l));
            planarSpectralAlbedo[i] = Math.pow(arr / r0, r0 / u2);
        }

        // NEW: this is for wvl < 400nm:
        final double psa1020 = planarSpectralAlbedo[72];
        final double psa560 = planarSpectralAlbedo[26];
        final double psa400 = planarSpectralAlbedo[10];

        final double c = ((psa560 - psa400) * (1.02 - 0.4) - (psa1020 - psa400) * (0.56 - 0.4)) /
                ((0.56 * 0.56 - 0.4 * 0.4) * (1.02 - 0.4) - (1.02 * 1.02 - 0.4 * 0.4) * (0.56 - 0.4));
        final double b = (psa1020 - psa400) / (1.02 - 0.4) - (0.4 + 1.02) * c;
        final double a = psa1020 - b * 1.02 - c * 1.02 * 1.02;
        for (int i = 0; i < 10; i++) {
            final double wvl = wvlsFull[i];
            planarSpectralAlbedo[i] = a + b * wvl + c * wvl * wvl;
        }
        return planarSpectralAlbedo;
    }

    private static double computeR0RelErr(double r0, double[] brr, double deltaBrr) {
        // AK: 'technical_note_JUNE_20_2018.docx', eq. (4)
        final double eps1 = 1.547269070091289;
        final double eps2 = -0.547269070091289;
        double r0RelErr = Math.sqrt(eps1 * eps1 * deltaBrr * deltaBrr / (brr[2] * brr[2]) +
                                            eps2 * eps2 * deltaBrr * deltaBrr / (brr[3] * brr[3]));
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
            sum += w[i] * w[i] * deltaBrr * deltaBrr / (brr[i] * brr[i]);
        }
        final double mRelErr = Math.sqrt(sum);
        // make sure error is positive and does not exceed value itself
        return Math.min(Math.abs(m), Math.abs(mRelErr));
    }

    private static double computeFRelErr(double f, double[] brr, double[] h, double deltaBrr) {
        // AK: 'technical_note_JUNE_20_2018.docx', eq. (4)
        double sum = 0.0;
        for (int i = 0; i < brr.length; i++) {
            sum += h[i] * h[i] * deltaBrr * deltaBrr / (brr[i] * brr[i]);
        }
        final double fRelErr = Math.sqrt(sum);
        return Math.min(Math.abs(f), Math.abs(fRelErr));
    }

}

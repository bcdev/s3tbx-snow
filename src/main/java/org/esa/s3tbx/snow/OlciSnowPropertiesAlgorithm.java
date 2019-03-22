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
     * @return SpectralAlbedoResult spectralAlbedos : spherical and planar at 21 OLCI wavelengths
     */
    static SpectralAlbedoResult computeSpectralAlbedosClean(double[] brr,
                                                            double deltaBrr,
                                                            double sza,
                                                            double vza) {
        return computeCoarseGridSpectralAlbedoWithErrorEstimates(brr, deltaBrr, sza, vza, 0.0, false);
    }

    /**
     * Computes spectral spherical and planar albedos for POLLUTED snow using AK latest algo from Oct 2018
     *
     * @param brr      - subset of BRR spectrum  (BRR_01, BRR_06, BRR_21, BRR_17)
     * @param deltaBrr - assumed BRR uncertainty
     * @param sza      - sun zenith angle (deg)
     * @param vza      - view zenith angle (deg)
     * @param r0Thresh - r0 must not go beyond r0Thresh for polluted snow
     * @return SpectralAlbedoResult spectralAlbedos : spherical and planar at 21 OLCI wavelengths
     */
    static SpectralAlbedoResult computeSpectralAlbedosPolluted(double[] brr,
                                                               double deltaBrr,
                                                               double sza,
                                                               double vza, double r0Thresh) {
        return computeCoarseGridSpectralAlbedoWithErrorEstimates(brr, deltaBrr, sza, vza, r0Thresh, true);
    }

    /**
     * Computes spectral spherical and planar albedos for POLLUTED snow using AK latest algo from Oct 2018
     *
     * @param brr        - subset of BRR spectrum  (BRR_01, BRR_06, BRR_21, BRR_17)
     * @param deltaBrr   - assumed BRR uncertainty
     * @param sza        - solar zenith angle
     * @param vza        - view zenith angle
     * @param r0Thresh   - r0 must not go beyond r0Thresh for polluted snow
     * @param isPolluted - indicator for clean or polluted snow
     * @return SpectralAlbedoResult: spherical/planar albedos at 21 OLCI wavelengths; r0, f, l, m and their errors
     */
    static SpectralAlbedoResult computeCoarseGridSpectralAlbedoWithErrorEstimates(double[] brr, double deltaBrr,
                                                                                  double sza, double vza,
                                                                                  double r0Thresh, boolean isPolluted) {

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
//                spectralAlbedos[0][i] = Math.pow(arr / r0, r0 / (u1 * u2));   // spectral spherical albedo
//                spectralAlbedos[1][i] = Math.pow(arr / r0, r0 / u2);
                spectralAlbedos[0][i] = Math.pow(arr / r0Thresh, r0Thresh / (u1 * u2));   // AK, 20190314
                spectralAlbedos[1][i] = Math.pow(arr / r0Thresh, r0Thresh / u2);       // AK, 20190314
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

        double[] z = new double[4];
        for (int i = 0; i < z.length; i++) {
            z[i] = 1.0 / Math.log(brr[i] / r0);
        }
        final double nu_1 = 2.0 * eps_1 - 2.0 * eps_1 * z[3];
        final double nu_2 = 2.0 * eps_2 + 2.0 * eps_1 * 1.0 / Math.log(brr[3] / r0);

        final double r0RelErr = computeR0RelErr(r0, brr, deltaBrr);
        final double lRelErr = computeLRelErr(l, brr, nu_1, nu_2, deltaBrr);

        // relative error estimation for isPolluted snow following AK...
        if (isPolluted) {
            final double wvl_1 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[0];
            final double wvl_2 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[5];

            final double s = 0.0; // todo: unreadable, check with AK, set to 0.0 in the meantime

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

            final double mRelErr = computeMRelErr(m, brr, w, deltaBrr);
            final double fRelErr = computeFRelErr(f, brr, h, deltaBrr);

            return new SpectralAlbedoResult(spectralAlbedos, r0, f, l, m, r0RelErr, fRelErr, lRelErr, mRelErr);
        } else {
            return new SpectralAlbedoResult(spectralAlbedos, r0, f, l, m, r0RelErr, Float.NaN, lRelErr, Float.NaN);
        }
    }

    /**
     * todo: document
     *
     * @param brr
     * @param refractiveIndexTable
     * @param solarSpectrumExtendedTable
     * @param sza
     * @param vza
     * @param isPolluted
     * @param r0Thresh
     * @return
     */
    static SpectralAlbedoResult computeFineGridSpectralAlbedo(double[] brr,
                                                              RefractiveIndexTable refractiveIndexTable,
                                                              SolarSpectrumExtendedTable solarSpectrumExtendedTable,
                                                              double sza, double vza,
                                                              boolean isPolluted, double r0Thresh) {

        // Same as computeCoarseGridSpectralAlbedo,
        // but here on fine 5nm grid rather than OLCI wavelengths only. No error estimates here.
        double[] wvlsFull = solarSpectrumExtendedTable.getWvl();
        final int numWvl = wvlsFull.length;
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

        // f (in 1/mm):
        final double SK = refWvl[0] / refWvl[3];
        final double f = p1 * Math.pow(SK, m) / (x * x * l);

        for (int i = 0; i < numWvl; i++) {
            final double wvl = wvlsFull[i];
            final double chi = refractiveIndexTable.getRefractiveIndexImag(i);
            final double wvlNm = 1000. * wvl;
            final double alka = 4.0 * Math.PI * chi / wvlNm;
            if (isPolluted) {
                // FORTRAN:
                // sd9=f*(wsk/ws(4))**(-m)
                // TT(j)=alka(j)*1.e+6+sd9
                // arr(j)=r0*exp(-x*sqrt(TT(j)*l))
                // plane1= (arr(j)/r0)**(r0/u2)
                // spher1=(arr(j)/r0)**(r0/u2/u1)
                final double sd9 = f * Math.pow(wvlNm / refWvl[3], -m);
                final double tt = alka * 1.E6 + sd9;
                final double arr = r0 * Math.exp(-x * Math.sqrt(tt * l));
                spectralAlbedos[0][i] = Math.pow(arr / r0Thresh, r0Thresh / (u1 * u2));   // AK, 20190314
                spectralAlbedos[1][i] = Math.pow(arr / r0Thresh, r0Thresh / u2);       // AK, 20190314
            } else {
                // FORTRAN:
                // TT(j)=alka(j)*1.e+6
                // arr(j)=r0*exp(-x*sqrt(TT(j)*l))
                // plane=exp(-u1*sqrt(TT(j)*l))
                // spher=exp(-sqrt(TT(j)*l))
                final double tt = alka * 1.E6;
                spectralAlbedos[0][i] = Math.exp(-Math.sqrt(tt * l));        // spectral spherical albedo
                spectralAlbedos[1][i] = Math.exp(-u1 * Math.sqrt(tt * l));     // spectral planar albedo
            }
        }
        return new SpectralAlbedoResult(spectralAlbedos, r0, f, l, m, Float.NaN, Float.NaN, Float.NaN, Float.NaN);
    }

    /**
     * Computes spherical and planar BB albedo in VIS, NIR, SW
     *
     * @param brr
     * @param isPolluted
     * @param refractiveIndexTable
     * @param solarSpectrumExtendedTable
     * @param sza
     * @param vza
     * @param r0thresh
     * @return
     */
    static double[][] computeBroadbandAlbedo(double[] brr,
                                             boolean isPolluted,
                                             RefractiveIndexTable refractiveIndexTable,
                                             SolarSpectrumExtendedTable solarSpectrumExtendedTable,
                                             double sza, double vza, double r0thresh) {

        double[] bbaVis = new double[2];
        double[] bbaNir = new double[2];
        double[] bbaSw = new double[2];

        double[] wvlsFull = solarSpectrumExtendedTable.getWvl();
        final double[] fLambda = SnowUtils.computeFLambda(solarSpectrumExtendedTable, sza);

        SpectralAlbedoResult fineGridSpectralAlbedo = computeFineGridSpectralAlbedo(brr,
                                                                                    refractiveIndexTable,
                                                                                    solarSpectrumExtendedTable,
                                                                                    sza, vza,
                                                                                    isPolluted, r0thresh);
        double[][] spectralAlbedo = fineGridSpectralAlbedo.getSpectralAlbedos(); // spherical: 0; planar: 1

        double[][] spectralAlbedoClean = null;
        if (isPolluted) {
            SpectralAlbedoResult fineGridSpectralAlbedoClean =
                    computeFineGridSpectralAlbedo(brr,
                                                  refractiveIndexTable,
                                                  solarSpectrumExtendedTable, sza, vza,
                                                  false, r0thresh);
            spectralAlbedoClean = fineGridSpectralAlbedoClean.getSpectralAlbedos();
        }
        double[][] fLambdaTimesSpectralAlbedo = new double[2][wvlsFull.length];
        final double wvlMax =
                OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[OlciSnowPropertiesConstants.OLCI_NUM_WVLS - 1];
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < fLambdaTimesSpectralAlbedo[0].length; i++) {
                if (isPolluted && wvlsFull[i] > wvlMax) {
                    fLambdaTimesSpectralAlbedo[j][i] = spectralAlbedoClean[j][i] * fLambda[i];
                } else {
                    fLambdaTimesSpectralAlbedo[j][i] = spectralAlbedo[j][i] * fLambda[i];
                }
            }
        }

        for (int i = 0; i < 2; i++) {
            final double numeratorVis = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
                                                                    OlciSnowPropertiesConstants.BB_WVL_2,
                                                                    fLambdaTimesSpectralAlbedo[i],
                                                                    wvlsFull);
            final double denominatorVis = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
                                                                      OlciSnowPropertiesConstants.BB_WVL_2,
                                                                      fLambda,
                                                                      wvlsFull);
            final double numeratorNir = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_2,
                                                                    OlciSnowPropertiesConstants.BB_WVL_3,
                                                                    fLambdaTimesSpectralAlbedo[i],
                                                                    wvlsFull);
            final double denominatorNir = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_2,
                                                                      OlciSnowPropertiesConstants.BB_WVL_3,
                                                                      fLambda,
                                                                      wvlsFull);
            final double numeratorSw = numeratorVis + numeratorNir;
            final double denominatorSw = denominatorVis + denominatorNir;

            bbaVis[i] = numeratorVis / denominatorVis;
            bbaNir[i] = numeratorNir / denominatorNir;
            bbaSw[i] = numeratorSw / denominatorSw;
        }

        final double[] sphericalBBA = {bbaVis[0], bbaNir[0], bbaSw[0]};
        final double[] planarBBA = {bbaVis[1], bbaNir[1], bbaSw[1]};
        
        return new double[][]{sphericalBBA, planarBBA};
    }

    /**
     * Computes probability of photon absorption (PPA) at considered wavelengths.
     * Follows 'ppa_new2_TECHNICAL_NOTE_28_11_2018.doc' (AK, 20181128)
     *
     * @param spectralAlbedoResult - length from spectral snow albedo result
     * @param isPollutedSnow
     * @return double [] ppa
     */
    static double[] computeSpectralPPA_nov2018(SpectralAlbedoResult spectralAlbedoResult, boolean isPollutedSnow) {
        double[] ppa = new double[OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length];

        final double gamma = 3.0 / 64.0;
        final double l = spectralAlbedoResult.getL();
        for (int i = 0; i < OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length; i++) {
            final double lambda = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];   // microns
            if (isPollutedSnow) {
                final double sphericalSpectralAlbedo = spectralAlbedoResult.getSpectralAlbedos()[0][i];
                ppa[i] = gamma * Math.log(sphericalSpectralAlbedo) * Math.log(sphericalSpectralAlbedo);
            } else {
                final double chi = OlciSnowPropertiesConstants.ICE_REFR_INDEX[i];
                final double alpha = 4.0 * Math.PI * chi / lambda;
                ppa[i] = gamma * alpha * l;
            }
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

    //////////////////////////////// end of package local //////////////////////////////////////////////////////////////

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

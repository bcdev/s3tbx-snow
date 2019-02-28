package org.esa.s3tbx.snow;

import org.esa.s3tbx.snow.math.Integrator;
import org.esa.snap.core.util.math.MathUtils;

/**
 * @author olafd
 */
class OlciSicePropertiesAlgorithm {

    static int computeSnowFlags() {
        // todo
        return 0;
    }

    /**
     * todo
     * requires brr[0], brr[5], brr[9], brr[10], brr[20] (400, 560, 681, 709, 1020)
     *
     * @param brr
     * @param r0
     * @param xx
     * @return SiceSnowPropertiesResult
     */
    static SiceSnowPropertiesResult computeGeneralSnowProperties(double[] brr, double r0, double xx) {

        // snow grain size:
        final double alam4 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[20];  // 1.02;
        final double akap4 = 2.25E-6;
        final double alpha4 = 4. * Math.PI * akap4 / alam4;
        final double brr1020 = brr[20];

        final double effAbsLength = Math.log(brr1020 / r0) * Math.log(brr1020 / r0) / (xx * xx * alpha4);   // in microns
        final double effAbsLengthMillimeter = effAbsLength / 1000.0;

        final double psi = 0.06;       // one number in breadboard 'psi.dat'...
        final double grainDiam = effAbsLengthMillimeter * psi;    // in mm

        // snow specific area:
        final double grainDiamMetres = grainDiam / 1000.0;    // in metres
        final double snowSpecificArea = computeSnowSpecificArea(grainDiamMetres);

        // snow pollution
        // requires brr[0], brr[5], brr[9], brr[10] (400, 560, 681, 709)
        double relImpurityLoad = computeRelativeImpurityLoad(brr, r0, xx, effAbsLengthMillimeter);

        // fill result with numbers we have up to now
        return new SiceSnowPropertiesResult(effAbsLength, grainDiam, snowSpecificArea, relImpurityLoad, 0.0, null);
    }

    /**
     * Provides spectral spherical and planar albedos.
     * We need all 21 rtoa here, but only brr400!
     *
     * @param snowProperties - result array which should already contain at least effective absorption length
     * @param rtoa           - rtoa
     * @param brr            - brr
     * @param sza            - sza
     * @param vza            - vza
     * @param raa            - raa
     * @return void
     */
    static void computeSpectralAlbedos(SiceSnowPropertiesResult snowProperties,
                                       double[] rtoa, double[] brr,
                                       double sza, double vza, double raa) {
        final double camu1 = Math.cos(sza * MathUtils.DTOR);
        final double camu2 = Math.cos(vza * MathUtils.DTOR);
        final double samu1 = Math.sin(sza * MathUtils.DTOR);
        final double samu2 = Math.sin(vza * MathUtils.DTOR);
        final double co = -camu1 * camu2 + samu1 * samu2 * Math.cos(raa * MathUtils.DTOR);
        final double u1 = SnowUtils.computeU(camu1);
        final double u2 = SnowUtils.computeU(camu2);

        double r0a1 = falex1(camu1, camu2, co);

        double[][] spectralAlbedos = new double[2][OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length];
        // todo: ask Alex what the two ways of spectral albedo retrieval mean
        for (int i = 0; i < OlciSnowPropertiesConstants.OLCI_NUM_WVLS; i++) {
            if (rtoa[i] > r0a1) {
                r0a1 = rtoa[i];
            }
            // todo: is the following line correct?? Note that r0a1 can change while going through the loop of 21 rtoa!
            // todo: why rtoa here but brr below??
            final double rs = Math.pow(rtoa[i] / r0a1, r0a1 / (u1 * u2));
            final double rp = Math.pow(rs, u1);
            spectralAlbedos[0][i] = rs;
            spectralAlbedos[1][i] = rp;
        }

        final double thresh = r0a1 - 0.1;
        snowProperties.setR0a1Thresh(thresh);
        if (brr[0] >= thresh) {
            // in this case, override the spectral albedos computed above
            final double effAbsLength = snowProperties.getEffAbsLength();
            for (int i = 0; i < OlciSnowPropertiesConstants.OLCI_NUM_WVLS; i++) {
                final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[i];
                final double iceRefrIndex = OlciSnowPropertiesConstants.ICE_REFR_INDEX_SICE_OLCI[i];
                final double dega = effAbsLength * 4.0 * Math.PI * iceRefrIndex / wvl;
                final double sqrtDega = Math.sqrt(dega);
                final double rsalex = sqrtDega > 1.E-6 ? Math.exp(-sqrtDega) : 1.0;
                final double rpalex = Math.pow(rsalex, u1);
                spectralAlbedos[0][i] = rsalex;    // todo: what is physically new here ?
                spectralAlbedos[1][i] = rpalex;
            }
        }
        snowProperties.setSpectralAlbedos(spectralAlbedos);
    }

    static double[] computePlanarBroadbandAlbedo() {
        // todo
//        !         Step 4.1 planar BBA
//
//        !     planar BBA
//                NSOLO = 0
//        x1 = 0.4425
//        x2 = 0.70875
//        x3 = 1.020
//
//        y1 = rp(3)
//        y2 = rp(11)
//        y3 = rp(21)
//
//        d1 = (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1)
//        d2 = (x3 * x3 - x1 * x1) * (x2 - x1) - (x2 * x2 - x1 * x1) * (x3 - x1)
//
//        !           second order polynomial coefficients  for planar albedo:
//        as = d1 / d2
//        bs = (y3 - y2 - as * (x3 * x3 - x2 * x2)) / (x3 - x2)
//        cs = y3 - as * x3 * x3 - bs * x3
//
//        !           limits of integration
//        at = 0.3
//        bt = 2.4
//        aat = 0.7
//
//        call qsimp(fun1, at, bt, ss1)
//        call qsimp(fun2, at, bt, ss2)
//
//        answer1 = ss1 / ss2
//
//        call qsimp(fun1, at, aat, ss1)
//        call qsimp(fun2, at, aat, ss2)
//
//        answer2 = ss1 / ss2
//
//        call qsimp(fun1, aat, bt, ss1)
//        call qsimp(fun2, aat, bt, ss2)
//
//        answer3 = ss1 / ss2

//        numeratorVis = Integrator.integrateSimpson(OlciSnowPropertiesConstants.BB_WVL_1,
//                                                   OlciSnowPropertiesConstants.BB_WVL_2,
//                                                   fLambdaTimesPlanarSpectralAlbedo,
//                                                   wvlsFull);


        return null;
    }

    static double[] computeSphericalBroadbandAlbedo() {
        // todo
        return null;
    }

    static double computeSolarLightSpectrum() {
        // todo
        return 0;
    }

    /**
     * some magic maths by Alex... // todo: clarify what this means
     *
     * @param am1 - ?
     * @param am2 - ?
     * @param co  - cosine of scattering angle
     * @return - ?
     */
    private static double falex1(double am1, double am2, double co) {
        final double a = 1.247;
        final double b = 1.186;
        final double c = 5.157;
        final double a1 = 0.087;
        final double a2 = 0.014;
        final double scat = Math.acos(co) * MathUtils.RTOD;
        final double p = 11.1 * Math.exp(-a1 * scat) + 1.1 * Math.exp(-a2 * scat);

        return (a + b * (am1 + am2) + c * am1 * am2 + p) / 4. / (am1 + am2);
    }


    static double computeR0(double[] brr, double sza, double vza) {
        // todo
        final double alam3 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[16];  // 0.865
        final double alam4 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[20];  // 1.02;

        final double akap3 = 2.4E-7;
        final double akap4 = 2.25E-6;

        final double alpha3 = 4. * Math.PI * akap3 / alam3;
        final double alpha4 = 4. * Math.PI * akap4 / alam4;

        final double eps = Math.sqrt(alpha3 / alpha4);
        final double ax1 = 1. / (1. - eps);
        final double ax2 = 1. / (1. - 1. / eps);

        final double brr865 = brr[16];
        final double brr1020 = brr[20];

        return Math.pow(brr865, ax1) * Math.pow(brr1020, ax2);
    }

    static double computeXX(double r0, double sza, double vza) {
        final double amu1 = Math.cos(sza * MathUtils.DTOR);
        final double amu2 = Math.cos(vza * MathUtils.DTOR);

        final double u1 = SnowUtils.computeU(amu1);
        final double u2 = SnowUtils.computeU(amu2);

        return u1 * u2 / r0;
    }

//    static double fun1(double x, double brr400, double r0, double xx, double effAbsLength, double r0a1Thresh, double sza,
//                       double as, double bs, double cs, boolean planar) {
//
//        double[] a = new double[6];
//        final double[][] bbb = OlciSnowPropertiesConstants.BBB_COEFFS_SICE;
//        final double[] bbbBounds = OlciSnowPropertiesConstants.BBB_COEFFS_SICE_BOUNDS;
//        for (int i = 0; i < 9; i++) {
//            for (int j = 0; j < 6; j++) {
//                if (x >= bbbBounds[i] && x < bbbBounds[i + 1]) {
//                    a[j] = bbb[i][j];
//                }
//            }
//        }
//
//        double astra;
//        if (x < 0.4) {
//            astra = 2.0E-11;
//        } else {
//            astra = a[0] + a[1] * x + a[2] * x * x + a[3] * x * x * x + a[4] * x * x * x * x + a[5] * x * x * x * x * x;
//        }
//
//        final double dega = effAbsLength * 4.0 * Math.PI * astra / x;
//        final double sqrtDega = Math.sqrt(dega);
//        final double rsd = sqrtDega > 1.E-6 ? Math.exp(-sqrtDega) : 1.0;
//        final double camu1 = Math.cos(sza * MathUtils.DTOR);
//        final double um1 = SnowUtils.computeU(camu1);
//
//        double f1;
//        if (brr400 <= r0a1Thresh && x <= 1.02) {
//            f1 = as * x * x + bs * x + cs;
//        } else {
//            f1 = planar ? Math.pow(rsd, um1) : rsd;
//        }
//
//        final double p0 = 32.38;
//        final double p1 = -160140.33;
//        final double p2 = 7959.53;
//        final double t1 = 85.34 * 1.e-3;
//        final double t2 = 401.79 * 1.e-3;
//
//        double funcs;
//        if (x <= 0.4) {
//            funcs = p0 + p1 * Math.exp(-0.4 / t1) + p2 * Math.exp(-0.4 / t2);
//        } else {
//            funcs = p0 + p1 * Math.exp(-x / t1) + p2 * Math.exp(-x / t2);
//        }
//
//        return f1 * funcs;
//    }

    private static double computeRelativeImpurityLoad(double[] brr, double r0, double x, double effAbsLengthMillimeter) {
        final double wvl400 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[0];
        final double wvl560 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[5];

        final double brr400 = brr[0];
        final double brr560 = brr[5];

        final double ap1 = Math.log(brr400 / r0) * Math.log(brr400 / r0);
        final double ap2 = Math.log(brr560 / r0) * Math.log(brr560 / r0);

        final double ang1 = Math.log(ap1 / ap2);
        final double ang2 = Math.log(wvl560 / wvl400);
        final double ang = ang1 / ang2;

        final double af = ap1 * Math.pow(0.4, ang) / (effAbsLengthMillimeter * x * x);

        // todo: make the following nicer, taken 1:1 from breadboard
        int ipol = 1;
        int ntype = 4;

        if (ang < 0.5 || ang > 10.0) {
            ipol = 0;
        }
        if (ang >= 0.85 && ang < 1.15) {
            ntype = 1;
        }
        if (ang >= 1.15) {
            ntype = 2;
        }

        final double r709 = brr[10];
        final double r681 = brr[9];
        final double ratio = r709 / r681;

        if (ratio > 1.0) {
            ntype = 3;
        }

        final double bbsoot = 0.9;
        final double sootka = 0.47;
        final double shape = 1.6;
        final double absoot = 1000.0 * bbsoot * 4.0 * Math.PI * sootka;
        final double abdust = 600.0;

        double conc = 0.0;
        if (ntype == 1) {
            conc = af * shape / absoot;
        }
        if (ntype == 2) {
            conc = af * shape / abdust;
        }
        if (ipol == 0) {
            conc = 0.0;
        }
        final double a680 = 1.E7;
        final double b680 = 600.0;
        if (ntype == 3) {
            conc = a680 * a680 * Math.log(ratio) * Math.log(ratio) - b680;
        }
        return conc;
    }

    private static double computeSnowSpecificArea(double grainDiamMetres) {
        return 6.0 / (OlciSnowPropertiesConstants.RHO_ICE * grainDiamMetres);
    }

}

package org.esa.s3tbx.snow;

import org.esa.snap.core.util.math.MathUtils;

/**
 * @author olafd
 */
class OlciSicePropertiesAlgorithm {

    static int computeSnowFlags() {
        // todo
        return 0;
    }

    static SiceSnowPropertiesResult computeGeneralSnowProperties(double[] brr, double sza, double vza) {

        // snow grain size:
        final double r0 = computeR0(brr, sza, vza);

        final double alam4 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[20];  // 1.02;
        final double akap4 = 2.25E-6;
        final double alpha4 = 4. * Math.PI * akap4 / alam4;
        final double brr1020 = brr[20];

        final double amu1 = Math.cos(sza * MathUtils.DTOR);
        final double amu2 = Math.cos(vza * MathUtils.DTOR);

        final double u1 = SnowUtils.computeU(amu1);
        final double u2 = SnowUtils.computeU(amu2);

        final double x = u1 * u2 / r0;

        final double effAbsLength = Math.log(brr1020 / r0) * Math.log(brr1020 / r0) / (x * x * alpha4);   // in microns
        final double effAbsLengthMillimeter = effAbsLength / 1000.0;

        final double psi = 0.06;       // one number in breadboard 'psi.dat'...
        final double grainDiam = effAbsLengthMillimeter * psi;    // in mm

        // snow specific area:
        final double grainDiamMetres = grainDiam / 1000.0;    // in metres
        final double snowSpecificArea = computeSnowSpecificArea(grainDiamMetres);

        // snow pollution
        double relImpurityLoad = computeRelativeImpurityLoad(brr, r0, x, effAbsLengthMillimeter);

        return new SiceSnowPropertiesResult(null, grainDiam, snowSpecificArea, relImpurityLoad);
    }

    static SpectralAlbedoResult computeSpectralAlbedos() {
        // todo
        return null;
    }

    static double[] computePlanarBroadbandAlbedo() {
        // todo
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

    static double falex1() {
        // todo
        return 0;
    }


    private static double computeR0(double[] brr, double sza, double vza) {
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

package org.esa.s3tbx.snow;

import org.esa.s3tbx.snow.math.Integrator;
import org.esa.s3tbx.snow.math.SiceFun1Function;
import org.esa.s3tbx.snow.math.SiceFun1InterpolInsideFunction;
import org.esa.s3tbx.snow.math.SiceFun2Function;
import org.esa.snap.core.gpf.Tile;
import org.esa.snap.core.util.math.MathUtils;

/**
 * @author olafd
 */
class OlciSiceSnowPropertiesAlgorithm {

    static int computeSnowFlags() {
        // todo : NDSI, NDBI
//        !         NDSI
//        andsi=(r(17)-r(21))/(r(17)+r(21))
//        if (aNDSI.gt.0.03.and.r(1).gt.0.5) indexs=1
//
//        !         NDBI
//        andbi=(r(1)-r(21))/(r(1)+r(21))
//        if (aNDBI.gt.0.33) indexi=1
//        if (aNDBI.gt.0.66) indexd=1


        return 0;
    }

    /**
     * Provides general snow properties:
     * - effective absorption length
     * - grain diameter
     * - snow specific area
     * - relative impurity load
     *
     * @param brr400  -
     * @param brr560  -
     * @param brr681  -
     * @param brr709  -
     * @param brr1020 -
     * @param r0      -
     * @param xx      -
     * @return SiceSnowPropertiesResult
     */
    static SiceSnowPropertiesResult computeGeneralSnowProperties(double brr400,
                                                                 double brr560,
                                                                 double brr681,
                                                                 double brr709,
                                                                 double brr1020,
                                                                 double r0,
                                                                 double xx) {

        // snow grain size:
        final double alam4 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[20];  // 1.02;
        final double akap4 = 2.25E-6;
        final double alpha4 = 4. * Math.PI * akap4 / alam4;

        final double effAbsLength = Math.log(brr1020 / r0) * Math.log(brr1020 / r0) / (xx * xx * alpha4);   // in microns
        final double effAbsLengthMillimeter = effAbsLength / 1000.0;

        final double psi = 0.06;       // one number in breadboard 'psi.dat'...
        final double grainDiam = effAbsLengthMillimeter * psi;    // in mm

        // snow specific area:
        final double grainDiamMetres = grainDiam / 1000.0;    // in metres
        final double snowSpecificArea = computeSnowSpecificArea(grainDiamMetres);

        // snow impurity
        // requires brr400, brr560, brr681, brr709
        SiceSnowImpurity snowImpurity =
                computeSnowImpurity(brr400, brr560, brr681, brr709, r0, xx, effAbsLengthMillimeter);

        // fill result with numbers we have up to now
        return new SiceSnowPropertiesResult(effAbsLength, grainDiam, snowSpecificArea, snowImpurity, 0.0, 0.0, null, null);
    }

    /**
     * Provides spectral spherical and planar albedos.
     * We need all 21 rtoa here, but only brr400!
     *
     * @param snowProperties - result array which should already contain at least effective absorption length
     * @param rtoa           - rtoa
     * @param brr400         - brr400
     * @param sza            - sza
     * @param vza            - vza
     * @param raa            - raa
     */
    static void computeSpectralAlbedos(SiceSnowPropertiesResult snowProperties,
                                       double[] rtoa, double brr400,
                                       double sza, double vza, double raa) {
        final double camu1 = Math.cos(sza * MathUtils.DTOR);
        final double camu2 = Math.cos(vza * MathUtils.DTOR);
        final double samu1 = Math.sin(sza * MathUtils.DTOR);
        final double samu2 = Math.sin(vza * MathUtils.DTOR);
        final double co = -camu1 * camu2 + samu1 * samu2 * Math.cos(raa * MathUtils.DTOR);
        final double u1 = SnowUtils.computeU(camu1);
        final double u2 = SnowUtils.computeU(camu2);

        final double scat = Math.acos(co) * MathUtils.RTOD;
        snowProperties.setScatteringAngle(scat);
        double r0a1 = falex1(camu1, camu2, scat);

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
        if (brr400 >= thresh) {
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
        snowProperties.setSphericalSpectralAlbedos(spectralAlbedos[0]);
        snowProperties.setPlanarSpectralAlbedos(spectralAlbedos[1]);
    }

    /**
     * Provides broadband albedos
     * @param snowProperties -
     * @param brr400         -
     * @param sza            -
     * @param refractiveIndexTable
     * @param wvlFullGrid    -  0.3 + i*0.005 in [0.3, 2.4], todo: set up array in initialize method!
     * @param astraFullGrid
     */
    static void computeBroadbandAlbedos(SiceSnowPropertiesResult snowProperties,
                                        double brr400,
                                        double sza,
                                        RefractiveIndexTable refractiveIndexTable,
                                        double[] wvlFullGrid,
                                        double[] astraFullGrid) {
        computePlanarBroadbandAlbedo(snowProperties, brr400, sza, refractiveIndexTable, wvlFullGrid, astraFullGrid);
        computeSphericalBroadbandAlbedo(snowProperties, brr400, sza, refractiveIndexTable, wvlFullGrid, astraFullGrid);
    }

    static double computeR0(double brr865, double brr1020) {
        final double alam3 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[16];  // 0.865
        final double alam4 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[20];  // 1.02;

        final double akap3 = 2.4E-7;
        final double akap4 = 2.25E-6;

        final double alpha3 = 4. * Math.PI * akap3 / alam3;
        final double alpha4 = 4. * Math.PI * akap4 / alam4;

        final double eps = Math.sqrt(alpha3 / alpha4);
        final double ax1 = 1. / (1. - eps);
        final double ax2 = 1. / (1. - 1. / eps);

        return Math.pow(brr865, ax1) * Math.pow(brr1020, ax2);
    }

    static double computeXX(double r0, double sza, double vza) {
        final double amu1 = Math.cos(sza * MathUtils.DTOR);
        final double amu2 = Math.cos(vza * MathUtils.DTOR);

        final double u1 = SnowUtils.computeU(amu1);
        final double u2 = SnowUtils.computeU(amu2);

        return u1 * u2 / r0;
    }

    static void setPollutionTypeFlag(int x, int y,
                                     Tile sicePollutionFlagTile,
                                     SiceSnowPropertiesResult siceSnowProperties,
                                     double ndbi) {
        if (ndbi <= 0.33) {
            // otherwise we have clean ice,i.e. no pollution at all
            final int pollutionType = siceSnowProperties.getSnowImpurity().getPollutionType();
            sicePollutionFlagTile.setSample(x, y, pollutionType, true);
        }
    }

    static void setGroundTypeFlag(int x, int y,
                                  Tile siceGroundFlagTile,
                                  SiceSnowPropertiesResult siceSnowProperties,
                                  double rtoa400, double rtoa1020,
                                  double ndsi, double ndbi) {

        if (ndsi > 0.03 && rtoa400 > 0.5) {
            siceGroundFlagTile.setSample(x, y, OlciSnowPropertiesConstants.SICE_SNOW, true);
        } else if (ndbi > 0.33) {
            siceGroundFlagTile.setSample(x, y, OlciSnowPropertiesConstants.SICE_BARE_ICE_CLEAN, true);
        } else if (ndbi > 0.66) {
            siceGroundFlagTile.setSample(x, y, OlciSnowPropertiesConstants.SICE_BARE_ICE_POLLUTED, true);
        } else if (siceSnowProperties.getSnowGrainSize() <= 0.01 || (rtoa400 - rtoa1020 <= 0.13)) {
            siceGroundFlagTile.setSample(x, y, OlciSnowPropertiesConstants.SICE_UNCERTAIN, true);
        }
    }

    static double getAstra(double x, double[] xa, double[] ya) {
        double mindiff = Double.MAX_VALUE;
        int kss = 0;
        for (int k = 0; k < xa.length; k++) {
            double dif = Math.abs(x - xa[k]);
            if (dif < mindiff) {
                mindiff = dif;
                kss = k;
                if (dif <= 1.E-6) {
                    return ya[k];
                }
                break;
            }
            if (dif <= 1.E-6) {
                return ya[k];
            }
        }
        double delta = x - xa[kss];
        int kss2 = delta <= 0 ? kss - 1 : kss;
        double x0, x1, y0, y1;
        x0 = xa[kss2];
        x1 = xa[kss2 + 1];
        y0 = ya[kss2];
        y1 = ya[kss + 1];

        return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    }

    static double getAstra2(double x, double[] xa, double[] ya) {
        double mindiff = Double.MAX_VALUE;
        int kss = 0;
        for (int k = 0; k < xa.length; k++) {
            double dif = Math.abs(x - xa[k]);
            if (dif < mindiff) {
                mindiff = dif;
                kss = k;
                if (dif <= 1.E-6) {
                    return ya[k];
                }
                break;
            }
            if (dif <= 1.E-6) {
                return ya[k];
            }
        }
        double delta = x - xa[kss];
        int kss2 = delta <= 0 ? kss - 1 : kss;
        double x0, x1, y0, y1;
        x0 = xa[kss2];
        x1 = xa[kss2 + 1];
        y0 = ya[kss2];
        y1 = ya[kss + 1];

        return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    }


    private static void computePlanarBroadbandAlbedo(SiceSnowPropertiesResult snowProperties,
                                                     double brr400,
                                                     double sza,
                                                     RefractiveIndexTable refractiveIndexTable,
                                                     double[] wvlFullGrid,
                                                     double[] astraFullGrid) {
        final double x1 = 0.4425;
        final double x2 = 0.70875;
        final double x3 = 1.020;

        final double y1 = snowProperties.getPlanarSpectralAlbedos()[2];
        final double y2 = snowProperties.getPlanarSpectralAlbedos()[10];
        final double y3 = snowProperties.getPlanarSpectralAlbedos()[20];

        final double d1 = (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1);
        final double d2 = (x3 * x3 - x1 * x1) * (x2 - x1) - (x2 * x2 - x1 * x1) * (x3 - x1);

        // second order polynomial coefficients for planar albedo:
        final double as = d1 / d2;
        final double bs = (y3 - y2 - as * (x3 * x3 - x2 * x2)) / (x3 - x2);
        final double cs = y3 - as * x3 * x3 - bs * x3;

        // limits of integration
        final double at = OlciSnowPropertiesConstants.BB_WVL_1;
        final double aat = OlciSnowPropertiesConstants.BB_WVL_2;
        final double bt = OlciSnowPropertiesConstants.BB_WVL_3;

//        final SiceFun1Function fun1 = new SiceFun1Function();
        final double[] xa = refractiveIndexTable.getWvl();
        final double[] ya = refractiveIndexTable.getRefractiveIndexImag();
        final SiceFun1InterpolInsideFunction fun1 = new SiceFun1InterpolInsideFunction(xa, ya);
        final SiceFun2Function fun2 = new SiceFun2Function();



        // fun1 params are:
        // double brr400, double effAbsLength, double r0a1Thresh, double cosSza,
        // double as, double bs, double cs, double planar
        final double effAbsLength = snowProperties.getEffAbsLength();
        final double r0a1Thresh = snowProperties.getR0a1Thresh();
        final double camu1 = Math.cos(sza * MathUtils.DTOR);
        double[] paramsFun1Planar = new double[]{brr400, effAbsLength, r0a1Thresh, camu1, as, bs, cs, 1.0};
        double[] paramsFun2 = new double[]{}; // no parameters needed

        final double numeratorVisPlanar = Integrator.integrateSimpsonSiceAlex(at, aat, fun1, paramsFun1Planar, wvlFullGrid);
        final double denominatorVisPlanar = Integrator.integrateSimpsonSiceAlex(at, aat, fun2, paramsFun2, wvlFullGrid);
//        paramsFun1Planar = new double[]{brr400, effAbsLength, r0a1Thresh, camu1, as, bs, cs, 1.0, 0.0};
//        paramsFun2 = new double[]{0.0};
//        final double numeratorVisPlanar = Integrator.integrateSimpsonSiceOlaf(at, aat, fun1, paramsFun1Planar, wvlFullGrid, astraFullGrid);
//        final double denominatorVisPlanar = Integrator.integrateSimpsonSiceOlaf(at, aat, fun2, paramsFun2, wvlFullGrid, astraFullGrid);
        final double bbVisPlanar = numeratorVisPlanar / denominatorVisPlanar;

        final double numeratorNirPlanar = Integrator.integrateSimpsonSiceAlex(aat, bt, fun1, paramsFun1Planar, wvlFullGrid);
        final double denominatorNirPlanar = Integrator.integrateSimpsonSiceAlex(aat, bt, fun2, paramsFun2, wvlFullGrid);
//        final double numeratorNirPlanar = Integrator.integrateSimpsonSiceOlaf(aat, bt, fun1, paramsFun1Planar, wvlFullGrid, astraFullGrid);
//        final double denominatorNirPlanar = Integrator.integrateSimpsonSiceOlaf(aat, bt, fun2, paramsFun2, wvlFullGrid, astraFullGrid);
        final double bbNirPlanar = numeratorNirPlanar / denominatorNirPlanar;

        final double numeratorSwPlanar = numeratorVisPlanar + numeratorNirPlanar;
        final double denominatorSwPlanar = denominatorVisPlanar + denominatorNirPlanar;
        final double bbSwPlanar = numeratorSwPlanar / denominatorSwPlanar;

        final double[] planarBBAlbedo = new double[]{bbVisPlanar, bbNirPlanar, bbSwPlanar};

        snowProperties.setPlanarBroadbandAlbedos(planarBBAlbedo);
    }

    private static void computeSphericalBroadbandAlbedo(SiceSnowPropertiesResult snowProperties,
                                                        double brr400,
                                                        double sza,
                                                        RefractiveIndexTable refractiveIndexTable, double[] wvlFullGrid,
                                                        double[] astraFullGrid) {
        final double x1 = 0.4425;
        final double x2 = 0.70875;
        final double x3 = 1.020;

        final double y1 = snowProperties.getSphericalSpectralAlbedos()[2];
        final double y2 = snowProperties.getSphericalSpectralAlbedos()[10];
        final double y3 = snowProperties.getSphericalSpectralAlbedos()[20];

        final double d1 = (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1);
        final double d2 = (x3 * x3 - x1 * x1) * (x2 - x1) - (x2 * x2 - x1 * x1) * (x3 - x1);

        // second order polynomial coefficients for planar albedo:
        final double as = d1 / d2;
        final double bs = (y3 - y2 - as * (x3 * x3 - x2 * x2)) / (x3 - x2);
        final double cs = y3 - as * x3 * x3 - bs * x3;

        // limits of integration
        final double at = OlciSnowPropertiesConstants.BB_WVL_1;
        final double aat = OlciSnowPropertiesConstants.BB_WVL_2;
        final double bt = OlciSnowPropertiesConstants.BB_WVL_3;

        final SiceFun1Function fun1 = new SiceFun1Function();
        final SiceFun2Function fun2 = new SiceFun2Function();

        // fun1 params are:
        // double brr400, double effAbsLength, double r0a1Thresh, double cosSza,
        // double as, double bs, double cs, double planar
        final double effAbsLength = snowProperties.getEffAbsLength();
        final double r0a1Thresh = snowProperties.getR0a1Thresh();
        final double camu1 = Math.cos(sza * MathUtils.DTOR);
        double[] paramsFun1Spherical = new double[]{brr400, effAbsLength, r0a1Thresh, camu1, as, bs, cs, 0.0};
        double[] paramsFun2 = new double[]{}; // no parameters needed

        final double numeratorVisSpherical = Integrator.integrateSimpsonSiceAlex(at, aat, fun1, paramsFun1Spherical, wvlFullGrid);
        final double denominatorVisSpherical = Integrator.integrateSimpsonSiceAlex(at, aat, fun2, paramsFun2, wvlFullGrid);
//        paramsFun1Spherical = new double[]{brr400, effAbsLength, r0a1Thresh, camu1, as, bs, cs, 0.0, 0.0};
//        paramsFun2 = new double[]{0.0}; // no parameters needed
//        final double numeratorVisSpherical = Integrator.integrateSimpsonSiceOlaf(at, aat, fun1, paramsFun1Spherical, wvlFullGrid, astraFullGrid);
//        final double denominatorVisSpherical = Integrator.integrateSimpsonSiceOlaf(at, aat, fun2, paramsFun2, wvlFullGrid, astraFullGrid);
        final double bbVisSpherical = numeratorVisSpherical / denominatorVisSpherical;

        final double numeratorNirSpherical = Integrator.integrateSimpsonSiceAlex(aat, bt, fun1, paramsFun1Spherical, wvlFullGrid);
        final double denominatorNirSpherical = Integrator.integrateSimpsonSiceAlex(aat, bt, fun2, paramsFun2, wvlFullGrid);
//        final double numeratorNirSpherical = Integrator.integrateSimpsonSiceOlaf(aat, bt, fun1, paramsFun1Spherical, wvlFullGrid, astraFullGrid);
//        final double denominatorNirSpherical = Integrator.integrateSimpsonSiceOlaf(aat, bt, fun2, paramsFun2, wvlFullGrid, astraFullGrid);
        final double bbNirSpherical = numeratorNirSpherical / denominatorNirSpherical;

        final double numeratorSwSpherical = numeratorVisSpherical + numeratorNirSpherical;
        final double denominatorSwSpherical = denominatorVisSpherical + denominatorNirSpherical;
        final double bbSwSpherical = numeratorSwSpherical / denominatorSwSpherical;

        final double[] sphericalBBAlbedo = new double[]{bbVisSpherical, bbNirSpherical, bbSwSpherical};

        snowProperties.setSphericalBroadbandAlbedos(sphericalBBAlbedo);
    }

    /**
     * some magic maths by Alex... // todo: clarify what this means
     */
    private static double falex1(double am1, double am2, double scat) {
        final double a = 1.247;
        final double b = 1.186;
        final double c = 5.157;
        final double a1 = 0.087;
        final double a2 = 0.014;
        final double p = 11.1 * Math.exp(-a1 * scat) + 1.1 * Math.exp(-a2 * scat);

        return (a + b * (am1 + am2) + c * am1 * am2 + p) / 4. / (am1 + am2);
    }

    private static SiceSnowImpurity computeSnowImpurity(double brr400, double brr560, double brr681, double brr709,
                                                        double r0, double x, double effAbsLengthMillimeter) {
        final double wvl400 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[0];
        final double wvl560 = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[5];

        final double ap1 = Math.log(brr400 / r0) * Math.log(brr400 / r0);
        final double ap2 = Math.log(brr560 / r0) * Math.log(brr560 / r0);

        final double ang1 = Math.log(ap1 / ap2);
        final double ang2 = Math.log(wvl560 / wvl400);
        final double ang = ang1 / ang2;

        final double af = ap1 * Math.pow(0.4, ang) / (effAbsLengthMillimeter * x * x);

        // todo: make the following nicer, taken 1:1 from breadboard
        int ipol = 1;    // this means polluted snow
        int ntype = 0;   // meaning of ntype: type of pollutants --> 1-soot, 2-dust, 3-algae, 0-uncertain

        if (ang < 0.5 || ang > 10.0) {
            ipol = 0;
        }
        if (ang >= 0.85 && ang < 1.15) {
            ntype = 1;
        }
        if (ang >= 1.15) {
            ntype = 2;
        }

        final double ratio = brr709 / brr681;

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
        if (ipol == 0) {    // this means clean snow
            conc = 0.0;
        }
        final double a680 = 1.E7;
        final double b680 = 600.0;
        if (ntype == 3) {
            conc = a680 * a680 * Math.log(ratio) * Math.log(ratio) - b680;
        }

        // absorption Angstroem exponent, normalized absorption coefficient at 1000nm, concentration of pollutants, flag
        return new SiceSnowImpurity(conc, ang, af, ntype);
    }

    private static double computeSnowSpecificArea(double grainDiamMetres) {
        return 6.0 / (OlciSnowPropertiesConstants.RHO_ICE * grainDiamMetres);
    }

}

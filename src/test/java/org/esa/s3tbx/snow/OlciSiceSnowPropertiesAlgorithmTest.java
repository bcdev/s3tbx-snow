package org.esa.s3tbx.snow;

import org.esa.s3tbx.snow.math.Integrator;
import org.esa.s3tbx.snow.math.SiceFun1Function;
import org.esa.s3tbx.snow.math.SiceFun2Function;
import org.junit.Before;
import org.junit.Test;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;


public class OlciSiceSnowPropertiesAlgorithmTest {

    private double sza_1;
    private double vza_1;
    private double saa_1;
    private double vaa_1;
    private double[] rtoa_1;
    private double[] brr_1;

    private double sza_5;
    private double vza_5;
    private double saa_5;
    private double vaa_5;
    private double[] rtoa_5;
    private double[] brr_5;

    private double sza_30;
    private double vza_30;
    private double saa_30;
    private double vaa_30;
    private double[] rtoa_30;
    private double[] brr_30;

    private double[] wvlFullGrid;
    private int numWvl;

    @Before
    public void setUp() throws Exception {

//        0.3 + i*0.005 in [0.3, 2.4], todo: set up array in initialize method!
        numWvl = (int) ((OlciSnowPropertiesConstants.BB_WVL_3 - OlciSnowPropertiesConstants.BB_WVL_1) / 0.005 + 1);
        wvlFullGrid = new double[numWvl];
        for (int i = 0; i < numWvl; i++) {
            wvlFullGrid[i] = OlciSnowPropertiesConstants.BB_WVL_1 + i * 0.005;
        }

        // first test input: first lines of input1_TOA.dat (angles), input2_BOA.dat (BRR)
        sza_1 = 50.48;
        vza_1 = 52.06;
        saa_1 = 141.9;
        vaa_1 = 95.9;
        rtoa_1 = new double[]{0.9220E+00, 0.9240E+00, 0.9270E+00, 0.9180E+00, 0.8930E+00, 0.8200E+00, 0.8000E+00,
                0.8530E+00, 0.8620E+00, 0.8630E+00, 0.8470E+00, 0.8350E+00, 0.1950E+00, 0.3980E+00, 0.7430E+00,
                0.7990E+00, 0.7490E+00, 0.7030E+00, 0.5540E+00, 0.2860E+00, 0.4500E+00};
        brr_1 = new double[]{0.9830E+00, 0.9820E+00, 0.9760E+00, 0.9730E+00, 0.9650E+00, 0.9430E+00, 0.9150E+00,
                0.9070E+00, 0.9050E+00, 0.9000E+00, 0.8500E+00, 0.8400E+00, 0.1840E+00, 0.3940E+00, 0.7440E+00,
                0.8000E+00, 0.7440E+00, 0.6980E+00, 0.5490E+00, 0.2810E+00, 0.4460E+00};

        // a second test input: fifth lines of input1_TOA.dat (angles), input2_BOA.dat (BRR)
        sza_5 = 49.63;
        vza_5 = 53.91;
        saa_5 = 140.4;
        vaa_5 = 95.1;
        rtoa_5 = new double[]{0.8550E+00, 0.8610E+00, 0.8690E+00, 0.8670E+00, 0.8480E+00, 0.7870E+00, 0.7700E+00,
                0.8180E+00, 0.8250E+00, 0.8270E+00, 0.8060E+00, 0.7960E+00, 0.2550E+00, 0.4330E+00, 0.7150E+00,
                0.7660E+00, 0.7250E+00, 0.6890E+00, 0.5460E+00, 0.3210E+00, 0.4780E+00};
        brr_5 = new double[]{0.8940E+00, 0.9010E+00, 0.9060E+00, 0.9140E+00, 0.9110E+00, 0.8990E+00, 0.8750E+00,
                0.8670E+00, 0.8640E+00, 0.8600E+00, 0.8090E+00, 0.8000E+00, 0.2440E+00, 0.4290E+00, 0.7150E+00,
                0.7680E+00, 0.7200E+00, 0.6840E+00, 0.5410E+00, 0.3160E+00, 0.4730E+00};

        // third test input for 'retrieval 1' of spectral albedos...
        sza_30 = 42.31;
        vza_30 = 41.0;
        saa_30 = 144.2;
        vaa_30 = 99.92;
        rtoa_30 = new double[]{0.7540E+00, 0.7570E+00, 0.7670E+00, 0.7700E+00, 0.7560E+00, 0.7120E+0, 0.6950E+00,
                0.7240E+00, 0.7290E+00, 0.7290E+00, 0.7030E+00, 0.6800E+00, 0.1580E+00, 0.3030E+00, 0.5820E+00,
                0.6410E+00, 0.5860E+00, 0.5420E+00, 0.4010E+00, 0.1710E+00, 0.3230E+00};

        brr_30 = new double[]{0.7800E+00, 0.7830E+00, 0.7920E+00, 0.8030E+00, 0.8000E+00, 0.7930E+00, 0.7710E+00,
                0.7580E+00, 0.7560E+00, 0.7520E+00, 0.7000E+00, 0.6820E+00, 0.1490E+00, 0.2980E+00, 0.5810E+00,
                0.6410E+00, 0.5820E+00, 0.5380E+00, 0.3980E+00, 0.1670E+00, 0.3200E+00};
    }

    @Test
    public void testComputeSnowFlags() {
        // todo
    }

    @Test
    public void testComputeR0() {
        assertEquals(0.9856, OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1, sza_1, vza_1), 1.E-3);
    }

    @Test
    public void testComputeXX() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1, sza_1, vza_1);
        assertEquals(0.9443, OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1), 1.E-3);
    }

    @Test
    public void testComputeSnowGrainSize() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1, sza_1, vza_1);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1);
        SiceSnowPropertiesResult generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_1, r0, xx);
        // result: first line, third column of output_flags.dat
        assertEquals(1.526, generalSnowProperties.getSnowGrainSize(), 1.E-3);

        r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_5, sza_5, vza_5);
        xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_5, vza_5);
        generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_5, r0, xx);
        // result: fifth line, third column of output_flags.dat
        assertEquals(0.895, generalSnowProperties.getSnowGrainSize(), 1.E-3);
    }

    @Test
    public void testComputeSnowSpecificArea() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1, sza_1, vza_1);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1);
        SiceSnowPropertiesResult generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_1, r0, xx);
        // result: first line, fourth column of output_flags.dat
        assertEquals(4.287, generalSnowProperties.getSnowSpecificArea(), 1.E-3);

        r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_5, sza_5, vza_5);
        xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_5, vza_5);
        generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_5, r0, xx);
        // result: fifth line, fourth column of output_flags.dat
        assertEquals(7.311, generalSnowProperties.getSnowSpecificArea(), 1.E-3);
    }

    @Test
    public void testComputeRelativeImpurityLoad() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1, sza_1, vza_1);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1);
        SiceSnowPropertiesResult generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_1, r0, xx);
        // result: first line, fourth column of output_impurity.dat  (clean case)
        assertEquals(0.0, generalSnowProperties.getRelativeImpurityLoad(), 1.E-3);

        r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_5, sza_5, vza_5);
        xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_5, vza_5);
        generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_5, r0, xx);
        // result: fifth line, fourth column of output_impurity.dat
        assertEquals(0.255E-8, generalSnowProperties.getRelativeImpurityLoad(), 1.E-3);
    }

    @Test
    public void testComputeSpectralAlbedos() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1, sza_1, vza_1);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1);
        SiceSnowPropertiesResult siceSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_1, r0, xx);
        double raa = SnowUtils.getRelAziSice(saa_1, vaa_1);
        OlciSiceSnowPropertiesAlgorithm.computeSpectralAlbedos(siceSnowProperties, rtoa_1, brr_1, sza_1, vza_1, raa);
        double[][] spectralAlbedos = siceSnowProperties.getSpectralAlbedos();
        assertNotNull(spectralAlbedos);
        assertEquals(2, spectralAlbedos.length);
        assertEquals(OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length, spectralAlbedos[0].length);
        // spherical, 'retrieval 2' (brr[0] >= thresh):
        checkSphericalSpectralAlbedos_retrieval2(spectralAlbedos[0]);
        // planar, 'retrieval 2' (brr[0] >= thresh):
        checkPlanarSpectralAlbedos_retrieval2(spectralAlbedos[1]);

        // 'retrieval 1':
        r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_30, sza_30, vza_30);
        xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_30, vza_30);
        siceSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_30, r0, xx);
        raa = SnowUtils.getRelAziSice(saa_30, vaa_30);
        OlciSiceSnowPropertiesAlgorithm.computeSpectralAlbedos(siceSnowProperties, rtoa_30, brr_30, sza_30, vza_30, raa);
        spectralAlbedos = siceSnowProperties.getSpectralAlbedos();
        assertNotNull(spectralAlbedos);
        assertEquals(2, spectralAlbedos.length);
        assertEquals(OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length, spectralAlbedos[0].length);
        // spherical, 'retrieval 1' (brr[0] < thresh):
        checkSphericalSpectralAlbedos_retrieval1(spectralAlbedos[0]);
        // planar, 'retrieval 1' (brr[0] < thresh):
        checkPlanarSpectralAlbedos_retrieval1(spectralAlbedos[1]);

    }

    private void checkPlanarSpectralAlbedos_retrieval1(double[] spectralAlbedo) {
        assertEquals(0.7609, spectralAlbedo[0], 1.E-3);
        assertEquals(0.7637, spectralAlbedo[1], 1.E-3);
        assertEquals(0.7732, spectralAlbedo[2], 1.E-3);
        assertEquals(0.7760, spectralAlbedo[3], 1.E-3);
        assertEquals(0.7628, spectralAlbedo[4], 1.E-3);
        assertEquals(0.7210, spectralAlbedo[5], 1.E-3);
        assertEquals(0.7049, spectralAlbedo[6], 1.E-3);
        assertEquals(0.7325, spectralAlbedo[7], 1.E-3);
        assertEquals(0.7372, spectralAlbedo[8], 1.E-3);
        assertEquals(0.7372, spectralAlbedo[9], 1.E-3);
        assertEquals(0.7125, spectralAlbedo[10], 1.E-3);
        assertEquals(0.6906, spectralAlbedo[11], 1.E-3);
        assertEquals(0.1756, spectralAlbedo[12], 1.E-3);
        assertEquals(0.3235, spectralAlbedo[13], 1.E-3);
        assertEquals(0.5968, spectralAlbedo[14], 1.E-3);
        assertEquals(0.6534, spectralAlbedo[15], 1.E-3);
        assertEquals(0.6007, spectralAlbedo[16], 1.E-3);
        assertEquals(0.5582, spectralAlbedo[17], 1.E-3);
        assertEquals(0.4208, spectralAlbedo[18], 1.E-3);
        assertEquals(0.1891, spectralAlbedo[19], 1.E-3);
        assertEquals(0.3435, spectralAlbedo[20], 1.E-3);
    }

    private void checkSphericalSpectralAlbedos_retrieval1(double[] spectralAlbedo) {
        assertEquals(0.7732, spectralAlbedo[0], 1.E-3);
        assertEquals(0.7759, spectralAlbedo[1], 1.E-3);
        assertEquals(0.7850, spectralAlbedo[2], 1.E-3);
        assertEquals(0.7877, spectralAlbedo[3], 1.E-3);
        assertEquals(0.7750, spectralAlbedo[4], 1.E-3);
        assertEquals(0.7350, spectralAlbedo[5], 1.E-3);
        assertEquals(0.7196, spectralAlbedo[6], 1.E-3);
        assertEquals(0.7460, spectralAlbedo[7], 1.E-3);
        assertEquals(0.7505, spectralAlbedo[8], 1.E-3);
        assertEquals(0.7505, spectralAlbedo[9], 1.E-3);
        assertEquals(0.7268, spectralAlbedo[10], 1.E-3);
        assertEquals(0.7058, spectralAlbedo[11], 1.E-3);
        assertEquals(0.1945, spectralAlbedo[12], 1.E-3);
        assertEquals(0.3457, spectralAlbedo[13], 1.E-3);
        assertEquals(0.6152, spectralAlbedo[14], 1.E-3);
        assertEquals(0.6699, spectralAlbedo[15], 1.E-3);
        assertEquals(0.6189, spectralAlbedo[16], 1.E-3);
        assertEquals(0.5777, spectralAlbedo[17], 1.E-3);
        assertEquals(0.4427, spectralAlbedo[18], 1.E-3);
        assertEquals(0.2086, spectralAlbedo[19], 1.E-3);
        assertEquals(0.3658, spectralAlbedo[20], 1.E-3);
    }

    private void checkPlanarSpectralAlbedos_retrieval2(double[] spectralAlbedo) {
        assertEquals(0.9958, spectralAlbedo[0], 1.E-3);
        assertEquals(0.9956, spectralAlbedo[1], 1.E-3);
        assertEquals(0.9931, spectralAlbedo[2], 1.E-3);
        assertEquals(0.9841, spectralAlbedo[3], 1.E-3);
        assertEquals(0.9784, spectralAlbedo[4], 1.E-3);
        assertEquals(0.9615, spectralAlbedo[5], 1.E-3);
        assertEquals(0.9373, spectralAlbedo[6], 1.E-3);
        assertEquals(0.9138, spectralAlbedo[7], 1.E-3);
        assertEquals(0.9105, spectralAlbedo[8], 1.E-3);
        assertEquals(0.9078, spectralAlbedo[9], 1.E-3);
        assertEquals(0.8879, spectralAlbedo[10], 1.E-3);
        assertEquals(0.8535, spectralAlbedo[11], 1.E-3);
        assertEquals(0.8452, spectralAlbedo[12], 1.E-3);
        assertEquals(0.8398, spectralAlbedo[13], 1.E-3);
        assertEquals(0.8359, spectralAlbedo[14], 1.E-3);
        assertEquals(0.8218, spectralAlbedo[15], 1.E-3);
        assertEquals(0.7482, spectralAlbedo[16], 1.E-3);
        assertEquals(0.7025, spectralAlbedo[17], 1.E-3);
        assertEquals(0.6865, spectralAlbedo[18], 1.E-3);
        assertEquals(0.6555, spectralAlbedo[19], 1.E-3);
        assertEquals(0.4414, spectralAlbedo[20], 1.E-3);
    }

    private void checkSphericalSpectralAlbedos_retrieval2(double[] spectralAlbedo) {
        assertEquals(0.9957, spectralAlbedo[0], 1.E-3);
        assertEquals(0.9954, spectralAlbedo[1], 1.E-3);
        assertEquals(0.9929, spectralAlbedo[2], 1.E-3);
        assertEquals(0.9836, spectralAlbedo[3], 1.E-3);
        assertEquals(0.9778, spectralAlbedo[4], 1.E-3);
        assertEquals(0.9605, spectralAlbedo[5], 1.E-3);
        assertEquals(0.9357, spectralAlbedo[6], 1.E-3);
        assertEquals(0.9117, spectralAlbedo[7], 1.E-3);
        assertEquals(0.9082, spectralAlbedo[8], 1.E-3);
        assertEquals(0.9055, spectralAlbedo[9], 1.E-3);
        assertEquals(0.8851, spectralAlbedo[10], 1.E-3);
        assertEquals(0.8499, spectralAlbedo[11], 1.E-3);
        assertEquals(0.8414, spectralAlbedo[12], 1.E-3);
        assertEquals(0.8359, spectralAlbedo[13], 1.E-3);
        assertEquals(0.8319, spectralAlbedo[14], 1.E-3);
        assertEquals(0.8175, spectralAlbedo[15], 1.E-3);
        assertEquals(0.7424, spectralAlbedo[16], 1.E-3);
        assertEquals(0.6959, spectralAlbedo[17], 1.E-3);
        assertEquals(0.6796, spectralAlbedo[18], 1.E-3);
        assertEquals(0.6481, spectralAlbedo[19], 1.E-3);
        assertEquals(0.4318, spectralAlbedo[20], 1.E-3);
    }

    @Test
    public void testComputeBroadbandAlbedo() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_30, sza_30, vza_30);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_30, vza_30);
        SiceSnowPropertiesResult siceSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_30, r0, xx);
        double raa = SnowUtils.getRelAziSice(saa_1, vaa_1);
        OlciSiceSnowPropertiesAlgorithm.computeSpectralAlbedos(siceSnowProperties, rtoa_30, brr_30, sza_30, vza_30, raa);
        OlciSiceSnowPropertiesAlgorithm.computePlanarBroadbandAlbedo(siceSnowProperties, brr_30[0], sza_30, wvlFullGrid);
        OlciSiceSnowPropertiesAlgorithm.computeSphericalBroadbandAlbedo(siceSnowProperties, brr_30[0], sza_30, wvlFullGrid);
        assertNotNull(siceSnowProperties.getPlanarBroadbandAlbedos());
        assertEquals(3, siceSnowProperties.getPlanarBroadbandAlbedos().length);
        assertEquals(0.57, siceSnowProperties.getPlanarBroadbandAlbedos()[0], 1.E-2);
        assertEquals(0.76, siceSnowProperties.getPlanarBroadbandAlbedos()[1], 1.E-2);
        assertEquals(0.37, siceSnowProperties.getPlanarBroadbandAlbedos()[2], 1.E-2);
        assertNotNull(siceSnowProperties.getSphericalBroadbandAlbedos());
        assertEquals(3, siceSnowProperties.getSphericalBroadbandAlbedos().length);
        assertEquals(0.58, siceSnowProperties.getSphericalBroadbandAlbedos()[0], 1.E-2);
        assertEquals(0.77, siceSnowProperties.getSphericalBroadbandAlbedos()[1], 1.E-2);
        assertEquals(0.37, siceSnowProperties.getSphericalBroadbandAlbedos()[2], 1.E-2);

        r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1, sza_1, vza_1);
        xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1);
        siceSnowProperties = OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties(brr_1, r0, xx);
        raa = SnowUtils.getRelAziSice(saa_1, vaa_1);
        OlciSiceSnowPropertiesAlgorithm.computeSpectralAlbedos(siceSnowProperties, rtoa_1, brr_1, sza_1, vza_1, raa);
        OlciSiceSnowPropertiesAlgorithm.computePlanarBroadbandAlbedo(siceSnowProperties, brr_1[0], sza_1, wvlFullGrid);
        OlciSiceSnowPropertiesAlgorithm.computeSphericalBroadbandAlbedo(siceSnowProperties, brr_1[0], sza_1, wvlFullGrid);
        assertNotNull(siceSnowProperties.getPlanarBroadbandAlbedos());
        assertEquals(3, siceSnowProperties.getPlanarBroadbandAlbedos().length);
        assertEquals(0.68, siceSnowProperties.getPlanarBroadbandAlbedos()[0], 4.E-2);     // 0.642
        assertEquals(0.97, siceSnowProperties.getPlanarBroadbandAlbedos()[1], 2.E-2);     // 0.969
        assertEquals(0.35, siceSnowProperties.getPlanarBroadbandAlbedos()[2], 1.E-2);     // 0.298
        assertNotNull(siceSnowProperties.getSphericalBroadbandAlbedos());
        assertEquals(3, siceSnowProperties.getSphericalBroadbandAlbedos().length);
        assertEquals(0.67, siceSnowProperties.getSphericalBroadbandAlbedos()[0], 1.E-2);    // 0.646
        assertEquals(0.98, siceSnowProperties.getSphericalBroadbandAlbedos()[1], 1.E-2);    // 0.970
        assertEquals(0.35, siceSnowProperties.getSphericalBroadbandAlbedos()[2], 1.E-2);    // 0.304
        System.out.println();
    }

    @Test
    public void testFun1() {
        double x = 2.38671637;  // a very large wavelength, 2368nm
        // params:
        // double brr400, double effAbsLength, double r0a1Thresh, double cosSza,
        // double as, double bs, double cs, double planar
        double brr400 = 0.819;
        double effAbsLength = 130.823105;
        double r0a1Thresh = 0.787318468;
        double cosSza = 0.749880135;
        double as = -0.110871777;
        double bs = 0.085336022;
        double cs = 0.778457224;
        double planar = 1.0;
        double[] params = new double[]{brr400, effAbsLength, r0a1Thresh, cosSza, as, bs, cs, planar};

        SiceFun1Function fun1 = new SiceFun1Function();
        double value = fun1.value(x, params);
        // at this wvl we obviously lose some precision in the breadboard with the 'real' numbers, so we need a large delta...
        assertEquals(26.8077, value, 1.0);  // Java: 26.2951

        x = 0.900050223;         // wavelength 900nm
        brr400 = 0.363;
        effAbsLength = 2952.22119;
        r0a1Thresh = 0.887;
        cosSza = 0.575433493;
        as = -0.743040442;
        bs = 0.625367403;
        cs = 0.319069743;
        planar = 0.0;
        fun1 = new SiceFun1Function();
        params = new double[]{brr400, effAbsLength, r0a1Thresh, cosSza, as, bs, cs, planar};
        value = fun1.value(x, params);
        assertEquals(245.1264, value, 1.E-3);      // this is better: 245.1263
    }

    @Test
    public void testFun2() {
        double x = 2.08124995;
        final double[] dummy = new double[]{0};
        SiceFun2Function fun2 = new SiceFun2Function();
        assertEquals(77.1787, fun2.value(x, dummy), 1.E-3);

        x = 0.753125012;
        fun2 = new SiceFun2Function();
        assertEquals(1230.1737, fun2.value(x, dummy), 1.E-3);
    }

    @Test
    public void testIntegrateSimpsonSice() {
        // limits of integration
        final double at = OlciSnowPropertiesConstants.BB_WVL_1;
        final double aat = OlciSnowPropertiesConstants.BB_WVL_2;
        final double bt = OlciSnowPropertiesConstants.BB_WVL_3;

        double x = 0.900050223;         // wavelength 900nm
        double brr400 = 0.363;
        double effAbsLength = 2952.22119;
        double r0a1Thresh = 0.887;
        double cosSza = 0.575433493;
        double as = -0.743040442;
        double bs = 0.625367403;
        double cs = 0.319069743;
        double planar = 0.0;
        SiceFun1Function fun1 = new SiceFun1Function();
        double[] params = new double[]{brr400, effAbsLength, r0a1Thresh, cosSza, as, bs, cs, planar};
        final double simpsonSice = Integrator.integrateSimpsonSice(at, bt, fun1, params, wvlFullGrid);
//        final double simpsonSiceAlex = Integrator.integrateSimpsonSiceAlex(at, bt, fun1, params, wvlFullGrid);
        System.out.println("simpsonSice = " + simpsonSice);
//        System.out.println("simpsonSiceAlex = " + simpsonSiceAlex);
    }


}

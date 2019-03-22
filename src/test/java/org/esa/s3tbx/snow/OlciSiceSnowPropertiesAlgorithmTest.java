package org.esa.s3tbx.snow;

import org.esa.s3tbx.snow.math.*;
import org.esa.snap.core.gpf.Tile;
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

    private double sza_xx;
    private double vza_xx;
    private double saa_xx;
    private double vaa_xx;
    private double[] rtoa_xx;
    private double[] brr_xx;

    private double[] wvlFullGrid;
    private double[] astraFullGrid;

    private RefractiveIndexTable refractiveIndexTable;

    @Before
    public void setUp() {

//        double wvlStep = 0.005;
        double wvlStep = 0.03;

//        0.3 + i*0.005 in [0.3, 2.4], todo: set up array in initialize method!
        int numWvl = (int) ((OlciSnowPropertiesConstants.BB_WVL_3 - OlciSnowPropertiesConstants.BB_WVL_1) / wvlStep + 1);
        wvlFullGrid = new double[numWvl];
        for (int i = 0; i < numWvl; i++) {
            wvlFullGrid[i] = OlciSnowPropertiesConstants.BB_WVL_1 + i * wvlStep;
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

        // new single line test input, 20190321
        sza_xx = 57.9;
        vza_xx = 1.38;
        saa_xx = 164.0;
        vaa_xx = 110.0;
        rtoa_xx = new double[]{0.331, 0.3277, 0.3199, 0.3167, 0.3166,
                0.3096, 0.3113, 0.3185, 0.3184, 0.3177, 0.3071, 0.3060, 0.0688, 0.1412, 0.2634, 0.2838, 0.2477, 0.2206,
                0.1713, 0.08725, 0.105965};

        brr_xx = new double[]{0.331, 0.3277, 0.3199, 0.3167, 0.3166, 0.3096, 0.3113, 0.3185, 0.3184, 0.3177, 0.3071,
                0.3060, 0.0688, 0.1412, 0.2634, 0.2838, 0.2477, 0.2206, 0.1713, 0.08725, 0.105965
        };


        refractiveIndexTable = new RefractiveIndexTable();
        astraFullGrid = SnowUtils.linearInterpolate(wvlFullGrid, refractiveIndexTable.getWvl(),
                                                    refractiveIndexTable.getRefractiveIndexImag());

    }

    @Test
    public void testComputeSnowFlags() {
        // todo
    }

    @Test
    public void testComputeR0() {
        assertEquals(0.9856, OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1[16], brr_1[20]), 1.E-3);
    }

    @Test
    public void testComputeXX() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1[16], brr_1[20]);
        assertEquals(0.9443, OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1), 1.E-3);
    }

    @Test
    public void testComputeSnowGrainSize() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1[16], brr_1[20]);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1);
        SiceSnowPropertiesResult generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_1[0], brr_1[5], brr_1[9], brr_1[10], brr_1[20], r0, xx);
        // result: first line, third column of output_flags.dat
        assertEquals(1.526, generalSnowProperties.getSnowGrainSize(), 1.E-3);

        r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_5[16], brr_5[20]);
        xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_5, vza_5);
        generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_5[0], brr_5[5], brr_5[9], brr_5[10], brr_5[20], r0, xx);
        // result: fifth line, third column of output_flags.dat
        assertEquals(0.895, generalSnowProperties.getSnowGrainSize(), 1.E-3);
    }

    @Test
    public void testComputeSnowSpecificArea() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1[16], brr_1[20]);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1);
        SiceSnowPropertiesResult generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_1[0], brr_1[5], brr_1[9], brr_1[10], brr_1[20], r0, xx);
        // result: first line, fourth column of output_flags.dat
        assertEquals(4.287, generalSnowProperties.getSnowSpecificArea(), 1.E-3);

        r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_5[16], brr_5[20]);
        xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_5, vza_5);
        generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_5[0], brr_5[5], brr_5[9], brr_5[10], brr_5[20], r0, xx);
        // result: fifth line, fourth column of output_flags.dat
        assertEquals(7.311, generalSnowProperties.getSnowSpecificArea(), 1.E-3);
    }

    @Test
    public void testComputeRelativeImpurityLoad() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1[16], brr_1[20]);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1);
        SiceSnowPropertiesResult generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_1[0], brr_1[5], brr_1[9], brr_1[10], brr_1[20], r0, xx);
        // result: first line, fourth column of output_impurity.dat  (clean case)
        assertEquals(0.0, generalSnowProperties.getSnowImpurity().getConcentrationOfPollutants(), 1.E-3);
        assertEquals(0, generalSnowProperties.getSnowImpurity().getPollutionType());

        r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_5[16], brr_5[20]);
        xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_5, vza_5);
        generalSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_5[0], brr_5[5], brr_5[9], brr_5[10], brr_5[20], r0, xx);
        // result: fifth line, fourth column of output_impurity.dat
        assertEquals(0.255E-8, generalSnowProperties.getSnowImpurity().getConcentrationOfPollutants(), 1.E-3);
        assertEquals(2, generalSnowProperties.getSnowImpurity().getPollutionType());
    }

    @Test
    public void testComputeSpectralAlbedos() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_1[16], brr_1[20]);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_1, vza_1);
        SiceSnowPropertiesResult siceSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_1[0], brr_1[5], brr_1[9], brr_1[10], brr_1[20], r0, xx);
        double raa = SnowUtils.getRelAziSice(saa_1, vaa_1);
        OlciSiceSnowPropertiesAlgorithm.computeSpectralAlbedos(siceSnowProperties, rtoa_1, brr_1[0], sza_1, vza_1, raa);
        // spherical, 'retrieval 2' (brr[0] >= thresh):
        double[] sphericalSpectralAlbedos = siceSnowProperties.getSphericalSpectralAlbedos();
        assertNotNull(sphericalSpectralAlbedos);
        assertEquals(OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length, sphericalSpectralAlbedos.length);
        checkSphericalSpectralAlbedos_retrieval2(sphericalSpectralAlbedos);
        // planar, 'retrieval 2' (brr[0] >= thresh):
        double[] planarSpectralAlbedos = siceSnowProperties.getPlanarSpectralAlbedos();
        assertNotNull(planarSpectralAlbedos);
        assertEquals(OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length, planarSpectralAlbedos.length);
        checkPlanarSpectralAlbedos_retrieval2(planarSpectralAlbedos);

        // 'retrieval 1':
        r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_30[16], brr_30[20]);
        xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_30, vza_30);
        siceSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_30[0], brr_30[5], brr_30[9], brr_30[10], brr_30[20], r0, xx);
        raa = SnowUtils.getRelAziSice(saa_30, vaa_30);
        OlciSiceSnowPropertiesAlgorithm.computeSpectralAlbedos(siceSnowProperties, rtoa_30, brr_30[0], sza_30, vza_30, raa);
        // spherical, 'retrieval 1' (brr[0] < thresh):
        sphericalSpectralAlbedos = siceSnowProperties.getSphericalSpectralAlbedos();
        assertNotNull(sphericalSpectralAlbedos);
        assertEquals(OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length, sphericalSpectralAlbedos.length);
        checkSphericalSpectralAlbedos_retrieval1(sphericalSpectralAlbedos);
        // planar, 'retrieval 1' (brr[0] < thresh):
        planarSpectralAlbedos = siceSnowProperties.getPlanarSpectralAlbedos();
        assertNotNull(planarSpectralAlbedos);
        assertEquals(OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length, planarSpectralAlbedos.length);
        checkPlanarSpectralAlbedos_retrieval1(planarSpectralAlbedos);

    }

    @Test
    public void testComputeSpectralAlbedos_xx() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_xx[16], brr_xx[20]);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_xx, vza_xx);
        SiceSnowPropertiesResult siceSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_xx[0], brr_xx[5], brr_xx[9], brr_xx[10], brr_xx[20], r0, xx);
        double raa = SnowUtils.getRelAziSice(saa_xx, vaa_xx);
        OlciSiceSnowPropertiesAlgorithm.computeSpectralAlbedos(siceSnowProperties, rtoa_xx, brr_xx[0], sza_xx, vza_xx, raa);
        // spherical, 'retrieval 2' (brr[0] >= thresh):
        double[] sphericalSpectralAlbedos = siceSnowProperties.getSphericalSpectralAlbedos();
        assertNotNull(sphericalSpectralAlbedos);
        assertEquals(OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length, sphericalSpectralAlbedos.length);
        checkSphericalSpectralAlbedos_retrieval_xx(sphericalSpectralAlbedos);
        // planar, 'retrieval 2' (brr[0] >= thresh):
        double[] planarSpectralAlbedos = siceSnowProperties.getPlanarSpectralAlbedos();
        assertNotNull(planarSpectralAlbedos);
        assertEquals(OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length, planarSpectralAlbedos.length);
        checkPlanarSpectralAlbedos_retrieval_xx(planarSpectralAlbedos);
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

    private void checkPlanarSpectralAlbedos_retrieval_xx(double[] spectralAlbedo) {
        assertEquals(0.4374, spectralAlbedo[0], 1.E-3);
        assertEquals(0.4341, spectralAlbedo[1], 1.E-3);
        assertEquals(0.4262, spectralAlbedo[2], 1.E-3);
        assertEquals(0.4229, spectralAlbedo[3], 1.E-3);
        assertEquals(0.4228, spectralAlbedo[4], 1.E-3);
        assertEquals(0.4157, spectralAlbedo[5], 1.E-3);
        assertEquals(0.4174, spectralAlbedo[6], 1.E-3);
        assertEquals(0.4248, spectralAlbedo[7], 1.E-3);
        assertEquals(0.4247, spectralAlbedo[8], 1.E-3);
        assertEquals(0.4240, spectralAlbedo[9], 1.E-3);
        assertEquals(0.4131, spectralAlbedo[10], 1.E-3);
        assertEquals(0.4120, spectralAlbedo[11], 1.E-3);
        assertEquals(0.1321, spectralAlbedo[12], 1.E-3);
        assertEquals(0.2285, spectralAlbedo[13], 1.E-3);
        assertEquals(0.3675, spectralAlbedo[14], 1.E-3);
        assertEquals(0.3890, spectralAlbedo[15], 1.E-3);
        assertEquals(0.3507, spectralAlbedo[16], 1.E-3);
        assertEquals(0.3211, spectralAlbedo[17], 1.E-3);
        assertEquals(0.2648, spectralAlbedo[18], 1.E-3);
        assertEquals(0.1583, spectralAlbedo[19], 1.E-3);
        assertEquals(0.1836, spectralAlbedo[20], 1.E-3);
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

    private void checkSphericalSpectralAlbedos_retrieval_xx(double[] spectralAlbedo) {
        assertEquals(0.3925, spectralAlbedo[0], 1.E-3);
        assertEquals(0.3891, spectralAlbedo[1], 1.E-3);
        assertEquals(0.3811, spectralAlbedo[2], 1.E-3);
        assertEquals(0.3778, spectralAlbedo[3], 1.E-3);
        assertEquals(0.3777, spectralAlbedo[4], 1.E-3);
        assertEquals(0.3705, spectralAlbedo[5], 1.E-3);
        assertEquals(0.3722, spectralAlbedo[6], 1.E-3);
        assertEquals(0.3796, spectralAlbedo[7], 1.E-3);
        assertEquals(0.3795, spectralAlbedo[8], 1.E-3);
        assertEquals(0.3788, spectralAlbedo[9], 1.E-3);
        assertEquals(0.3679, spectralAlbedo[10], 1.E-3);
        assertEquals(0.3668, spectralAlbedo[11], 1.E-3);
        assertEquals(0.1013, spectralAlbedo[12], 1.E-3);
        assertEquals(0.1883, spectralAlbedo[13], 1.E-3);
        assertEquals(0.3223, spectralAlbedo[14], 1.E-3);
        assertEquals(0.3437, spectralAlbedo[15], 1.E-3);
        assertEquals(0.3057, spectralAlbedo[16], 1.E-3);
        assertEquals(0.2766, spectralAlbedo[17], 1.E-3);
        assertEquals(0.2224, spectralAlbedo[18], 1.E-3);
        assertEquals(0.1243, spectralAlbedo[19], 1.E-3);
        assertEquals(0.1470, spectralAlbedo[20], 1.E-3);
    }

    @Test
    public void testComputeBroadbandAlbedo_xx() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_xx[16], brr_xx[20]);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_xx, vza_xx);
        SiceSnowPropertiesResult siceSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_xx[0], brr_xx[5], brr_xx[9], brr_xx[10], brr_xx[20], r0, xx);
        double raa = SnowUtils.getRelAziSice(saa_xx, vaa_xx);
        OlciSiceSnowPropertiesAlgorithm.computeSpectralAlbedos(siceSnowProperties, rtoa_xx, brr_xx[0], sza_xx, vza_xx, raa);
        OlciSiceSnowPropertiesAlgorithm.computeBroadbandAlbedos(siceSnowProperties, brr_xx[0], sza_xx,
                                                                refractiveIndexTable);
        assertNotNull(siceSnowProperties.getPlanarBroadbandAlbedos());
        assertEquals(3, siceSnowProperties.getPlanarBroadbandAlbedos().length);
        assertEquals(0.4231, siceSnowProperties.getPlanarBroadbandAlbedos()[0], 1.E-2);
        assertEquals(0.3409, siceSnowProperties.getPlanarBroadbandAlbedos()[1], 1.E-2);
        assertEquals(0.3835, siceSnowProperties.getPlanarBroadbandAlbedos()[2], 1.E-2);
        assertNotNull(siceSnowProperties.getSphericalBroadbandAlbedos());
        assertEquals(3, siceSnowProperties.getSphericalBroadbandAlbedos().length);
        assertEquals(0.3779, siceSnowProperties.getSphericalBroadbandAlbedos()[0], 1.E-2);
        assertEquals(0.3055, siceSnowProperties.getSphericalBroadbandAlbedos()[1], 1.E-2);
        assertEquals(0.343, siceSnowProperties.getSphericalBroadbandAlbedos()[2], 1.E-2);
    }

    @Test
    public void testComputeFlags_xx() {
        double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr_xx[16], brr_xx[20]);
        double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza_xx, vza_xx);
        SiceSnowPropertiesResult siceSnowProperties =
                OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                        (brr_xx[0], brr_xx[5], brr_xx[9], brr_xx[10], brr_xx[20], r0, xx);

        double ndsi = (rtoa_xx[16] - rtoa_xx[20]) / (rtoa_xx[16] + rtoa_xx[20]);
        double ndbi = (rtoa_xx[0] - rtoa_xx[20]) / (rtoa_xx[0] + rtoa_xx[20]);
        final int groundTypeFlag =
                OlciSiceSnowPropertiesAlgorithm.computeGroundTypeFlag(siceSnowProperties,
                                                                      rtoa_xx[0], rtoa_xx[20], ndsi, ndbi);
        assertEquals(2, groundTypeFlag);

        final int pollutionTypeFlag = OlciSiceSnowPropertiesAlgorithm.computePollutionTypeFlag(siceSnowProperties, ndbi);
        assertEquals(0, pollutionTypeFlag);
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
    public void testFun1Performance() {
        // here we show that FUN1 via internal astra interpolation is ~25% slower than using a static polynom!
        double x;
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

        SiceFun1InterpolInsideFunction fun1Interpol = new SiceFun1InterpolInsideFunction(refractiveIndexTable.getWvl(),
                                                                                         refractiveIndexTable.getRefractiveIndexImag());
        long t1 = System.currentTimeMillis();
        for (int i = 0; i < refractiveIndexTable.getWvl().length - 1; i++) {
            x = refractiveIndexTable.getWvl(i) + 0.002;
            for (int j = 0; j < 1.E5; j++) {
                fun1Interpol.value(x, params);
            }
        }
        long t2 = System.currentTimeMillis();
        System.out.println("t2-t1 for fun1 INTERPOL = " + (t2 - t1));

        SiceFun1PolynomFunction fun1Poly = new SiceFun1PolynomFunction(null, null);
        long t3 = System.currentTimeMillis();
        for (int i = 0; i < refractiveIndexTable.getWvl().length - 1; i++) {
            x = refractiveIndexTable.getWvl(i) + 0.002;
            for (int j = 0; j < 1.E5; j++) {
                fun1Poly.value(x, params);
            }
        }
        long t4 = System.currentTimeMillis();
        System.out.println("t4-t3 for fun1 POLY = " + (t4 - t3));
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
    public void test_fun1_interpolate() {
        double x = 0.4;
        double cosSza = 0.5;
        double as = 0.188057452;
        double bs = -0.967674673;
        double cs = 1.19137311;
        double brr400 = 0.95;
        double r0a1Thresh = 0.97;
        double effAbsLength = 14761.8789;
        double[] xa = refractiveIndexTable.getWvl();
        double[] ya = refractiveIndexTable.getRefractiveIndexImag();
        SiceFun1InterpolInsideFunction fun1 = new SiceFun1InterpolInsideFunction(xa, ya);
        double[] params = new double[]{brr400, effAbsLength, r0a1Thresh, cosSza, as, bs, cs, 1.0};
        double result = fun1.value(x, params);
        assertEquals(1250.103, result, 1.E-3);
    }


    @Test
    public void testIntegrateSimpsonSice_fun2() {
        // limits of integration
        final double at = OlciSnowPropertiesConstants.BB_WVL_1;
        final double bt = OlciSnowPropertiesConstants.BB_WVL_3;

        SiceFun2Function fun2 = new SiceFun2Function();
        double simpsonSiceAlex = 0.0;
        double[] params2 = new double[]{0};
        final long t3 = System.currentTimeMillis();
        for (int i = 0; i < 10000; i++) {
            simpsonSiceAlex = Integrator.integrateSimpsonSiceAlex(at, bt, fun2, params2);
        }
        final long t4 = System.currentTimeMillis();
        System.out.println("t4-t3 for fun2 Simpson SICE Alex = " + (t4 - t3));
        // for function 'fun2', the 'Simpson Alex' is almost the same as Fortran BB
        System.out.println("fun2 simpsonSiceAlex = " + simpsonSiceAlex);          // J: 1262.2802; F: 1262.2793
    }

    @Test
    public void testLinearInterpolation() {
        // preparation
        // x^2
        double[] x = {0, 10, 20, 30, 40, 50};
        double[] y = {0, 100, 400, 900, 1600, 2500};
        double[] xi = {16};
        double[] yi = SnowUtils.linearInterpolate(xi, x, y);
        assertEquals(1, yi.length);
        assertEquals(280.0, yi[0], 1.E-3);     // !!!
    }

    @Test
    public void testGetAstra() {
        double astra1 = SnowUtils.linearInterpolate(0.795, refractiveIndexTable.getWvl(),
                                                    refractiveIndexTable.getRefractiveIndexImag());

        double[] astra2 = SnowUtils.linearInterpolate(new double[]{0.795}, refractiveIndexTable.getWvl(),
                                                      refractiveIndexTable.getRefractiveIndexImag());

        System.out.println("astra1 = " + astra1);
        System.out.println("astra2 = " + astra2[0]);
        assertEquals(1.26E-7, astra1, 1.E-9);
        assertEquals(1.26E-7, astra2[0], 1.E-9);
    }

    @Test
    public void testGetAstraFullGrid() {
        double[] x = wvlFullGrid;
        double[] astra = SnowUtils.linearInterpolate(x, refractiveIndexTable.getWvl(),
                                                     refractiveIndexTable.getRefractiveIndexImag());
        assertNotNull(astra);
//        assertEquals(421, astra.length);
//        assertEquals(2.0E-011, astra[0], 1.E-12);
//        assertEquals(1.26E-7, astra[99], 1.E-9);
//        assertEquals(5.7805E-004, astra[420], 1.E-5);

        assertEquals(71, astra.length);
        assertEquals(2.0E-011, astra[0], 1.E-12);
        assertEquals(3.84E-6, astra[29], 1.E-9);
        assertEquals(5.7805E-004, astra[70], 1.E-5);
    }
}

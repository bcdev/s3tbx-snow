package org.esa.s3tbx.snow;

import org.junit.Before;
import org.junit.Test;

import static junit.framework.TestCase.assertEquals;


public class OlciSiceSnowPropertiesAlgorithmTest {

    private double sza_1;
    private double vza_1;
    private double saa_1;
    private double vaa_1;
    private double[] rtoa_1;
    private double[] brr_1;

    private double sza_2;
    private double vza_2;
    private double saa_2;
    private double vaa_2;
    private double[] rtoa_2;
    private double[] brr_2;

    @Before
    public void setUp() throws Exception {
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
        sza_2 = 49.63;
        vza_2 = 53.91;
        saa_2 = 140.4;
        vaa_2 = 95.1;
        rtoa_2 = new double[]{0.8550E+00, 0.8610E+00, 0.8690E+00, 0.8670E+00, 0.8480E+00, 0.7870E+00, 0.7700E+00,
                0.8180E+00, 0.8250E+00, 0.8270E+00, 0.8060E+00, 0.7960E+00, 0.2550E+00, 0.4330E+00, 0.7150E+00,
                0.7660E+00, 0.7250E+00, 0.6890E+00, 0.5460E+00, 0.3210E+00, 0.4780E+00};
        brr_2 = new double[]{0.8940E+00, 0.9010E+00, 0.9060E+00, 0.9140E+00, 0.9110E+00, 0.8990E+00, 0.8750E+00,
                0.8670E+00, 0.8640E+00, 0.8600E+00, 0.8090E+00, 0.8000E+00, 0.2440E+00, 0.4290E+00, 0.7150E+00,
                0.7680E+00, 0.7200E+00, 0.6840E+00, 0.5410E+00, 0.3160E+00, 0.4730E+00};
    }

    @Test
    public void testComputeSnowFlags() {
        // todo
    }

    @Test
    public void testComputeSnowGrainSize() {
        SiceSnowPropertiesResult generalSnowProperties =
                OlciSicePropertiesAlgorithm.computeGeneralSnowProperties(brr_1, sza_1, vza_1);
        // result: first line, third column of output_flags.dat
        assertEquals(1.526, generalSnowProperties.getSnowGrainSize(), 1.E-3);

        generalSnowProperties =
                OlciSicePropertiesAlgorithm.computeGeneralSnowProperties(brr_2, sza_2, vza_2);
        // result: fifth line, third column of output_flags.dat
        assertEquals(0.895, generalSnowProperties.getSnowGrainSize(), 1.E-3);
    }

    @Test
    public void testComputeSnowSpecificArea() {
        SiceSnowPropertiesResult generalSnowProperties =
                OlciSicePropertiesAlgorithm.computeGeneralSnowProperties(brr_1, sza_1, vza_1);
        // result: first line, fourth column of output_flags.dat
        assertEquals(4.287, generalSnowProperties.getSnowSpecificArea(), 1.E-3);

        generalSnowProperties =
                OlciSicePropertiesAlgorithm.computeGeneralSnowProperties(brr_2, sza_2, vza_2);
        // result: fifth line, fourth column of output_flags.dat
        assertEquals(7.311, generalSnowProperties.getSnowSpecificArea(), 1.E-3);
    }

    @Test
    public void testComputeRelativeImpurityLoad() {
        SiceSnowPropertiesResult generalSnowProperties =
                OlciSicePropertiesAlgorithm.computeGeneralSnowProperties(brr_1, sza_1, vza_1);
        // result: first line, fourth column of output_impurity.dat
        assertEquals(0.0, generalSnowProperties.getRelativeImpurityLoad(), 1.E-3);

        generalSnowProperties =
                OlciSicePropertiesAlgorithm.computeGeneralSnowProperties(brr_2, sza_2, vza_2);
        // result: fifth line, fourth column of output_impurity.dat
        assertEquals(0.255E-8, generalSnowProperties.getRelativeImpurityLoad(), 1.E-3);
    }

    @Test
    public void testComputeSpectralAlbedos() {
        // todo
    }

    @Test
    public void testComputePlanarBroadbandAlbedo() {
        // todo
    }

    @Test
    public void testComputeSphericalBroadbandAlbedo() {
        // todo
    }

    @Test
    public void testComputeSolarLightSpectrum() {
        // todo
    }

    @Test
    public void testFalex1() {
        // todo
    }


}

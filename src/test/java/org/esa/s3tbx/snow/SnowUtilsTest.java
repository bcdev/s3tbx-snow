package org.esa.s3tbx.snow;

import org.junit.Ignore;
import org.junit.Test;

import static junit.framework.Assert.assertNotNull;
import static junit.framework.TestCase.assertEquals;


public class SnowUtilsTest {

    @Test
    public void testInterpolateSpectralAlbedos() {
        //preparation
        double[] x = {0, 50, 100};
        double[] y = {0, 50, 200};
        double[] xi = {10, 25, 75, 80, 100};

        double[] yi = SnowUtils.splineInterpolate(x, y, xi);

        //assertion
        assertEquals(5, yi.length, 1.E-2);
        assertEquals(10.0, yi[0]);
        assertEquals(25.0, yi[1]);
        assertEquals(125.0, yi[2]);
        assertEquals(140.0, yi[3]);
        assertEquals(200.0, yi[4]);
    }

    @Test
    public void testSetupRcSourceBandNames() {
        String[] namesAlbedo = new String[]{"Oa01_radiance", "Oa17_radiance"};
        String[] namesPpa = new String[]{"Oa01_radiance", "Oa07_radiance", "Oa13_radiance", "Oa21_radiance"};
        final String[] rcSourceBands = SnowUtils.setupRcSourceBands(namesAlbedo, namesPpa);
        assertNotNull(rcSourceBands);
        assertEquals(5, rcSourceBands.length);
        assertEquals("Oa01_radiance", rcSourceBands[0]);
        assertEquals("Oa17_radiance", rcSourceBands[1]);
        assertEquals("Oa07_radiance", rcSourceBands[2]);
        assertEquals("Oa13_radiance", rcSourceBands[3]);
        assertEquals("Oa21_radiance", rcSourceBands[4]);
    }

    @Test
    public void testRoundTo4DecimalPlaces() {
        assertEquals(1.2346, SnowUtils.cutTo4DecimalPlaces(1.2345678), 1.E-8);
        assertEquals(7.6543, SnowUtils.cutTo4DecimalPlaces(7.654321), 1.E-8);
    }

    @Test
    @Ignore
    public void testGetInterpolFlux() {
        SolarSpectrumExtendedTable solarSpectrumExtendedTable = new SolarSpectrumExtendedTable();
        double[][] solarSpectrum = solarSpectrumExtendedTable.getSolarSpectrum();

        int wvlIndex = 0;   // first row in 'final_table_fluxes_angles.txt' file
        double sza = 60.0;
        int lowerIndex = (int) sza/15;
        int upperIndex = lowerIndex + 1;

        double interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(530.89, interpolFlux, 1.E-3);

        sza = 65.0;
        interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(366.816, interpolFlux, 1.E-3);

        sza = 67.5;
        interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(282.1417, interpolFlux, 1.E-3);

        sza = 70.0;
        interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(195.922, interpolFlux, 1.E-3);

        sza = 74.999999;
        interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(19.51, interpolFlux, 1.E-3);

        wvlIndex = 73;   // some row in 'final_table_fluxes_angles.txt' file
        sza = 15.0;
        lowerIndex = (int) sza/15;
        upperIndex = lowerIndex + 1;

        interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(1316.1, interpolFlux, 1.E-3);

        sza = 20.0;
        interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(1278.549, interpolFlux, 1.E-3);

        sza = 22.5;
        interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(1255.913, interpolFlux, 1.E-3);

        sza = 25.;
        interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(1230.761, interpolFlux, 1.E-3);

        sza = 30.0;
        interpolFlux = SnowUtils.getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, wvlIndex);
        assertEquals(1173.1, interpolFlux, 1.E-3);
    }

    @Test
    @Ignore
    public void testGetExtrapolFlux() {
        SolarSpectrumExtendedTable solarSpectrumExtendedTable = new SolarSpectrumExtendedTable();
        double[][] solarSpectrum = solarSpectrumExtendedTable.getSolarSpectrum();

        double sza = 75.0;
        int lowerIndex = (int) sza/15;
        double extrapolFlux = SnowUtils.getExtrapolFlux(solarSpectrum[lowerIndex][0], sza);
        assertEquals(19.51, extrapolFlux, 1.E-3);

        sza = 80.0;
        extrapolFlux = SnowUtils.getExtrapolFlux(solarSpectrum[lowerIndex][0], sza);
        assertEquals(13.089, extrapolFlux, 1.E-3);

        sza = 85.0;
        extrapolFlux = SnowUtils.getExtrapolFlux(solarSpectrum[lowerIndex][0], sza);
        assertEquals(6.569, extrapolFlux, 1.E-3);

        sza = 90.0;
        extrapolFlux = SnowUtils.getExtrapolFlux(solarSpectrum[lowerIndex][0], sza);
        assertEquals(0.0, extrapolFlux, 1.E-3);
    }
}

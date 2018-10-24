package org.esa.s3tbx.snow;

import org.junit.Ignore;
import org.junit.Test;

import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertNotNull;


public class SnowAuxdataTest {

    @Test
    public void testLoadRefractiveIndexData() {
        final RefractiveIndexTable refractiveIndexTable = new RefractiveIndexTable();
        assertNotNull(refractiveIndexTable);
        assertEquals(162, refractiveIndexTable.getRefractiveIndexImag().length);
        assertEquals(2.0E-11, refractiveIndexTable.getRefractiveIndexImag(0), 1.E-12);
        assertEquals(0.25, refractiveIndexTable.getWvl(0), 1.E-12);
        assertEquals(1.18E-7, refractiveIndexTable.getRefractiveIndexImag(43), 1.E-12);
        assertEquals(0.79, refractiveIndexTable.getWvl(43), 1.E-12);
        assertEquals(3.04E-6, refractiveIndexTable.getRefractiveIndexImag(80), 1.E-12);
        assertEquals(1.16, refractiveIndexTable.getWvl(80), 1.E-12);
        assertEquals(5.956E-4, refractiveIndexTable.getRefractiveIndexImag(161), 1.E-12);
        assertEquals(2.41, refractiveIndexTable.getWvl(161), 1.E-12);

    }


    @Test
    public void testLoadSolarSpectrumData() {
        // solar_spectrum_2:
        final SolarSpectrumTable solarSpectrumTable = new SolarSpectrumTable();
        assertNotNull(solarSpectrumTable);
        assertEquals(411, solarSpectrumTable.getSolarSpectrum().length);
        assertEquals(530.9, solarSpectrumTable.getSolarSpectrum(0), 1.E-12);
        assertEquals(0.35, solarSpectrumTable.getWvl(0), 1.E-12);
        assertEquals(584.8, solarSpectrumTable.getSolarSpectrum(84), 1.E-12);
        assertEquals(0.77, solarSpectrumTable.getWvl(84), 1.E-12);
        assertEquals(21.07, solarSpectrumTable.getSolarSpectrum(410), 1.E-12);
        assertEquals(2.4, solarSpectrumTable.getWvl(410), 1.E-12);
    }

    @Test
    @Ignore
    public void testLoadExtendedSolarSpectrumData() {
        // final_table_fluxes_angles:
        final SolarSpectrumExtendedTable solarSpectrumExtendedTable = new SolarSpectrumExtendedTable();
        assertNotNull(solarSpectrumExtendedTable);
        assertEquals(6, solarSpectrumExtendedTable.getSolarSpectrum().length);
        assertEquals(411, solarSpectrumExtendedTable.getSolarSpectrum()[0].length);

        assertEquals(1323.6, solarSpectrumExtendedTable.getSolarSpectrum(0, 0), 1.E-12);
        assertEquals(0.35, solarSpectrumExtendedTable.getWvl(0), 1.E-12);
        assertEquals(1214.1, solarSpectrumExtendedTable.getSolarSpectrum(0, 84), 1.E-12);
        assertEquals(0.77, solarSpectrumExtendedTable.getWvl(84), 1.E-12);
        assertEquals(45.899, solarSpectrumExtendedTable.getSolarSpectrum(0, 410), 1.E-12);
        assertEquals(2.4, solarSpectrumExtendedTable.getWvl(410), 1.E-12);

        assertEquals(1267.3, solarSpectrumExtendedTable.getSolarSpectrum(1, 0), 1.E-12);
        assertEquals(1171.3, solarSpectrumExtendedTable.getSolarSpectrum(1, 84), 1.E-12);
        assertEquals(44.184, solarSpectrumExtendedTable.getSolarSpectrum(1, 410), 1.E-12);

        assertEquals(1077.9, solarSpectrumExtendedTable.getSolarSpectrum(2, 3), 1.E-12);
        assertEquals(1084.7, solarSpectrumExtendedTable.getSolarSpectrum(2, 76), 1.E-12);
        assertEquals(94.081, solarSpectrumExtendedTable.getSolarSpectrum(2, 213), 1.E-12);

        assertEquals(1451.5, solarSpectrumExtendedTable.getSolarSpectrum(3, 21), 1.E-12);
        assertEquals(1227.5, solarSpectrumExtendedTable.getSolarSpectrum(3, 46), 1.E-12);
        assertEquals(593.6, solarSpectrumExtendedTable.getSolarSpectrum(3, 111), 1.E-12);

        assertEquals(903.23, solarSpectrumExtendedTable.getSolarSpectrum(4, 33), 1.E-12);
        assertEquals(356.27, solarSpectrumExtendedTable.getSolarSpectrum(4, 133), 1.E-12);
        assertEquals(131.47, solarSpectrumExtendedTable.getSolarSpectrum(4, 233), 1.E-12);

        assertEquals(202.25, solarSpectrumExtendedTable.getSolarSpectrum(5, 44), 1.E-12);
        assertEquals(53.157, solarSpectrumExtendedTable.getSolarSpectrum(5, 244), 1.E-12);
        assertEquals(13.134, solarSpectrumExtendedTable.getSolarSpectrum(5, 344), 1.E-12);
    }

    @Test
    @Ignore
    public void testLoadExtendedSolarSpectrumData_oct2018() {
        // final_table_fluxes_angles_oct2018:
        final SolarSpectrumExtendedTable solarSpectrumExtendedTable = new SolarSpectrumExtendedTable();
        assertNotNull(solarSpectrumExtendedTable);
        assertEquals(89, solarSpectrumExtendedTable.getSolarSpectrum().length);
        assertEquals(411, solarSpectrumExtendedTable.getSolarSpectrum()[0].length);

        assertEquals(1323.9, solarSpectrumExtendedTable.getSolarSpectrum(0, 0), 1.E-12);
        assertEquals(0.35, solarSpectrumExtendedTable.getWvl(0), 1.E-12);
        assertEquals(1214.0, solarSpectrumExtendedTable.getSolarSpectrum(0, 84), 1.E-12);
        assertEquals(0.77, solarSpectrumExtendedTable.getWvl(84), 1.E-12);
        assertEquals(45.868, solarSpectrumExtendedTable.getSolarSpectrum(0, 410), 1.E-12);
        assertEquals(2.4, solarSpectrumExtendedTable.getWvl(410), 1.E-12);

        assertEquals(1633.1, solarSpectrumExtendedTable.getSolarSpectrum(33, 33), 1.E-12);
        assertEquals(528.96, solarSpectrumExtendedTable.getSolarSpectrum(43, 133), 1.E-12);
        assertEquals(160.13, solarSpectrumExtendedTable.getSolarSpectrum(53, 233), 1.E-12);

        assertEquals(11.904, solarSpectrumExtendedTable.getSolarSpectrum(88, 44), 1.E-12);
        assertEquals(4.415, solarSpectrumExtendedTable.getSolarSpectrum(88, 244), 1.E-12);
        assertEquals(0.5293, solarSpectrumExtendedTable.getSolarSpectrum(88, 344), 1.E-12);
    }
}

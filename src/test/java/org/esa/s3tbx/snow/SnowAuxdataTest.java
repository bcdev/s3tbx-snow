package org.esa.s3tbx.snow;

import org.junit.Ignore;
import org.junit.Test;

import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertNotNull;


public class SnowAuxdataTest {

    @Test
    public void testLoadRefractiveIndexData() throws Exception {
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

//    @Test
//    public void testLoadSolarSpectrumData() throws Exception {
//        final SolarSpectrumTable solarSpectrumTable = new SolarSpectrumTable();
//        assertNotNull(solarSpectrumTable);
//        assertEquals(7923, solarSpectrumTable.getSolarSpectrum().length);
//        assertEquals(474.57, solarSpectrumTable.getSolarSpectrum(0), 1.E-12);
//        assertEquals(299.8E-3, solarSpectrumTable.getWvl(0), 1.E-12);
//        assertEquals(917.91, solarSpectrumTable.getSolarSpectrum(884), 1.E-12);
//        assertEquals(892.14E-3, solarSpectrumTable.getWvl(884), 1.E-12);
//        assertEquals(728.41, solarSpectrumTable.getSolarSpectrum(2097), 1.E-12);
//        assertEquals(1000.4E-3, solarSpectrumTable.getWvl(2097), 1.E-12);
//        assertEquals(60.45, solarSpectrumTable.getSolarSpectrum(7922), 1.E-12);
//        assertEquals(2397.51E-3, solarSpectrumTable.getWvl(7922), 1.E-12);
//
//    }

    @Test
    @Ignore
    public void testLoadSolarSpectrumData() throws Exception {
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
}

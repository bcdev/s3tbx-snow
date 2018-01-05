package org.esa.s3tbx.snow;

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
    public void testSetupRcSourceBandNames() throws Exception {
        String[] namesAlbedo = new String[]{"Oa01_radiance", "Oa17_radiance"};;
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
}

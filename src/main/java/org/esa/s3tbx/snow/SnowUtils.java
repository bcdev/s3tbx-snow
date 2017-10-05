package org.esa.s3tbx.snow;

import static org.esa.s3tbx.snow.OlciSnowAlbedoConstants.*;

/**
 * S3 Snow utility class
 *
 * @author olafd
 */
public class SnowUtils {

    public static String[] getReflectanceTypeBandNames(Sensor sensor, int reflType) {
        return sensor.getRequiredBrrBandNames();
//        if (reflType == SensorConstants.REFL_TYPE_BRR) {
//            return sensor.getRequiredBrrBandNames();
//        } else if (reflType == SensorConstants.REFL_TYPE_TOA) {
//            return sensor.getRequiredReflBandNames();
//        } else {
//            throw new IllegalArgumentException("reflType " + reflType + " not supported.");
//        }
    }

    public static double[] computeKappa2() {
        double[] kappa2 = new double[WAVELENGTH_GRID_OLCI.length];

        // first interval, bands 1-5
        for (int i = 0; i < 5; i++) {
            final double wvl = WAVELENGTH_GRID_OLCI[i];
            kappa2[i] = computeKappa2(wvl, 0);
        }
        // second interval, bands 6-9
        for (int i = 6; i < 9; i++) {
            final double wvl = WAVELENGTH_GRID_OLCI[i];
            kappa2[i] = computeKappa2(wvl, 1);
        }
        // third interval, bands 10-21
        for (int i = 10; i < WAVELENGTH_GRID_OLCI.length; i++) {
            final double wvl = WAVELENGTH_GRID_OLCI[i];
            kappa2[i] = computeKappa2(wvl, 2);
        }

        return kappa2;
    }

    public static double computeKappa2(double wvl, int rangeIndex) {
        if (rangeIndex < 0 || rangeIndex > 3) {
            throw new IllegalArgumentException("rangeIndex must be in [0,3]");
        }
        return c0[rangeIndex] + c1[rangeIndex] * (wvl - LAMBDA_0[rangeIndex]) / H[rangeIndex] +
                c2[rangeIndex] * Math.pow((wvl - LAMBDA_0[rangeIndex]) / H[rangeIndex], 2.0) +
                c3[rangeIndex] * Math.pow((wvl - LAMBDA_0[rangeIndex]) / H[rangeIndex], 3.0);
    }

}

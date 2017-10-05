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
            kappa2[i] = c0[0] + c1[0] * (wvl - LAMBDA_0[0]) / H[0] +
                    c2[0] * Math.pow((wvl - LAMBDA_0[0]) / H[0], 2.0) +
                    c3[0] * Math.pow((wvl - LAMBDA_0[0]) / H[0], 3.0);
        }
        // second interval, bands 6-9
        for (int i = 6; i < 9; i++) {
            final double wvl = WAVELENGTH_GRID_OLCI[i];
            kappa2[i] = c0[1] + c1[1] * (wvl - LAMBDA_0[1]) / H[1] +
                    c2[1] * Math.pow((wvl - LAMBDA_0[1]) / H[1], 2.0) +
                    c3[1] * Math.pow((wvl - LAMBDA_0[1]) / H[1], 3.0);
        }
        // third interval, bands 10-21
        for (int i = 10; i < WAVELENGTH_GRID_OLCI.length; i++) {
            final double wvl = WAVELENGTH_GRID_OLCI[i];
            kappa2[i] = c0[2] + c1[2] * (wvl - LAMBDA_0[2]) / H[2] +
                    c2[2] * Math.pow((wvl - LAMBDA_0[2]) / H[2], 2.0) +
                    c3[2] * Math.pow((wvl - LAMBDA_0[2]) / H[2], 3.0);
        }

        return kappa2;
    }

}

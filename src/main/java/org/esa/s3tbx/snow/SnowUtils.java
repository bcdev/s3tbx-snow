package org.esa.s3tbx.snow;

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

}

package org.esa.s3tbx.snow;

/**
 * Constants for supported sensors for snow albedo retrieval (just OLCI so far).
 *
 * @author olafd
 */
public class SensorConstants {

    public static final int REFL_TYPE_TOA = 1;
    public static final int REFL_TYPE_BRR = 2;

    public static final String OLCI_SZA_NAME = "SZA";
    public static final String OLCI_SAA_NAME = "SAA";
    public static final String OLCI_VZA_NAME = "OZA";
    public static final String OLCI_VAA_NAME = "OAA";
    public static final String OLCI_L1B_FLAGS_NAME = "quality_flags";
    public static final int OLCI_INVALID_BIT = 25;
    public static final int OLCI_LAND_BIT = 31;

    public final static String OLCI_BRR_BAND_PREFIX = "rBRR";

    public final static String[] OLCI_TARGET_TPGS = {
            "TP_latitude", "TP_longitude",
            "SZA", "SAA", "OZA", "OAA",
            "horizontal_wind_vector_1", "horizontal_wind_vector_2"
    };

}

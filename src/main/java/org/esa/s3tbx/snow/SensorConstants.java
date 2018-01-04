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
    public static final String OLCI_VZA_NAME = "OZA";
    public static final String OLCI_SAA_NAME = "SAA";
    public static final String OLCI_VAA_NAME = "OAA";
    public static final String OLCI_L1B_FLAGS_NAME = "quality_flags";
    public static final int OLCI_INVALID_BIT = 25;
    public static final int OLCI_LAND_BIT = 31;

    public final static String OLCI_REFL_BAND_SUFFIX = "reflectance";
    public final static String OLCI_BRR_BAND_PREFIX = "rBRR";

    // we need OLCI BRR bands 1, 2, 3, 4, 5, 12, 17, 21;
//    public final static String[] OLCI_REQUIRED_RADIANCE_BAND_NAMES = new String[]{
//            "Oa01_radiance", "Oa02_radiance", "Oa03_radiance", "Oa04_radiance", "Oa05_radiance",
//            "Oa12_radiance", "Oa17_radiance", "Oa21_radiance"
//    };
//
//    public final static String[] OLCI_REQUIRED_REFL_BAND_NAMES = new String[]{
//            "Oa01_reflectance", "Oa02_reflectance", "Oa03_reflectance", "Oa04_reflectance", "Oa05_reflectance",
//            "Oa12_reflectance", "Oa17_reflectance", "Oa21_reflectance"
//    };
//
//    public final static String[] OLCI_REQUIRED_BRR_BAND_NAMES = new String[]{
//            "rBRR_01", "rBRR_02", "rBRR_03", "rBRR_04", "rBRR_05",
//            "rBRR_12", "rBRR_17", "rBRR_21"
//    };

    // for SIMPLE_APPROXIMATION we only need OLCI BRR bands 1 + 18 or 21   (for 865 or 1020 choice)
    public final static String[] OLCI_DEFAULT_REQUIRED_RADIANCE_BAND_NAMES = new String[]{
            "Oa01_radiance", "Oa21_radiance"
    };

    public final static String[] OLCI_DEFAULT_REQUIRED_REFL_BAND_NAMES = new String[]{
            "Oa01_reflectance", "Oa21_reflectance"
    };

    public final static String[] OLCI_DEFAULT_REQUIRED_BRR_BAND_NAMES = new String[]{
            "rBRR_01", "rBRR_21"
    };

    public final static String[] OLCI_TARGET_TPGS = {
            "TP_latitude", "TP_longitude",
            "SZA", "SAA", "OZA", "OAA",
            "horizontal_wind_vector_1", "horizontal_wind_vector_2"
    };

}

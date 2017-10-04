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

    public final static String[] OLCI_REQUIRED_RADIANCE_BAND_NAMES = new String[]{
            "Oa01_radiance", "Oa02_radiance", "Oa03_radiance", "Oa04_radiance", "Oa05_radiance",
            "Oa06_radiance", "Oa07_radiance", "Oa08_radiance", "Oa09_radiance", "Oa10_radiance",
            "Oa11_radiance", "Oa12_radiance", "Oa13_radiance", "Oa14_radiance", "Oa15_radiance",
            "Oa16_radiance", "Oa17_radiance", "Oa18_radiance", "Oa19_radiance", "Oa20_radiance",
            "Oa21_radiance"
    };

    public final static String OLCI_REFL_BAND_SUFFIX = "reflectance";

    public final static String[] OLCI_REQUIRED_REFL_BAND_NAMES = new String[]{
            "Oa01_reflectance", "Oa02_reflectance", "Oa03_reflectance", "Oa04_reflectance", "Oa05_reflectance",
            "Oa06_reflectance", "Oa07_reflectance", "Oa08_reflectance", "Oa09_reflectance", "Oa10_reflectance",
            "Oa11_reflectance", "Oa12_reflectance", "Oa13_reflectance", "Oa14_reflectance", "Oa15_reflectance",
            "Oa16_reflectance", "Oa17_reflectance", "Oa18_reflectance", "Oa19_reflectance", "Oa20_reflectance",
            "Oa21_reflectance"
    };

    public final static String[] OLCI_REQUIRED_BRR_BAND_NAMES = new String[]{
            "rBRR_01", "rBRR_02", "rBRR_03", "rBRR_04", "rBRR_05",
            "rBRR_06", "rBRR_07", "rBRR_08", "rBRR_09", "rBRR_10",
            "rBRR_11", "rBRR_12", "rBRR_13", "rBRR_14", "rBRR_15",
            "rBRR_16", "rBRR_17", "rBRR_18", "rBRR_19", "rBRR_20",
            "rBRR_21"
    };
}

package org.esa.s3tbx.snow;

/**
 * Constants for OLCI snow albedo retrieval
 *
 * @author olafd
 */
public class OlciSnowPropertiesConstants {


    public static final double[] WAVELENGTH_GRID_OLCI = {
            0.400,         // OA01
            0.4125,        // OA02
            0.4425,        // OA03
            0.490,         // OA04
            0.510,         // OA05
            0.560,         // OA06
            0.620,         // OA07
            0.665,         // OA08
            0.67375,       // OA09
            0.68125,       // OA10
            0.70875,       // OA11
            0.75375,       // OA12
            0.76125,       // OA13
            0.764375,      // OA14
            0.7675,        // OA15
            0.77875,       // OA16
            0.865,         // OA17
            0.885,         // OA18
            0.900,         // OA19
            0.940,         // OA20
            1.02           // OA21
    };

    public static final int OLCI_NUM_WVLS = WAVELENGTH_GRID_OLCI.length;

    // Table 1 from spectral_albedo_october_5_2017.doc:
    // no longer used
//    static final double[] LAMBDA_0 = {0.35, 0.525, 0.7, 1.98};
//    static final double[] H = {0.175, 0.175, 1.28, 0.4675};
//    static final double[] c0 = {3.85E-9, 2.2E-9, 1.30689, 1.27525};
//    static final double[] c1 = {-5.287E-9, 3.384E-9, -0.026, -0.026};
//    static final double[] c2 = {2.2594E-9, 6.41E-9, 0.017956, -0.00149};
//    static final double[] c3 = {1.38411E-9, 1.63665E-8, -0.024, -0.0128};

    public static final double[] ICE_REFR_INDEX = {
            2.365E-011,2.669E-011,6.268E-011,4.172E-010,8.036E-010,
            2.839E-009,8.580E-009,1.770E-008,1.890E-008,2.090E-008,
            3.440E-008,5.870E-008,7.080E-008,7.500E-008,8.080E-008,
            1.020E-007,2.400E-007,3.600E-007,4.200E-007,5.530E-007,
            2.250E-006
    };

    // from file 'kap_olci.dat'. They are partly different from ICE_REFR_INDEX. todo: ask Alex why
    public static final double[] ICE_REFR_INDEX_SICE_OLCI = {
            2.365E-011,2.669E-011,7.0E-011,4.17E-010,8.04E-010,
            2.84E-009,8.580E-009,1.78E-008,1.95E-008,2.1E-008,
            3.3E-008,6.23E-008,7.1E-008,7.68E-008,8.13E-008,
            9.88E-008,2.400E-007,3.64E-007,4.200E-007,5.530E-007,
            2.250E-006
    };

    // from file 'kap_ice.dat'.
    public static final double[][] BBB_COEFFS_SICE = {
            {-4.37756e-6, 4.11459e-5, -1.53323e-4, 2.83421e-4, -2.60375e-4, 9.54029e-5},
            {0.01918, -0.10443, 0.2267, -0.2452, 0.13211, -0.02835},
            {0.36039, -1.62512, 2.92662, -2.6307, 1.18018, -0.21136},
            {-4.16499, 15.79931, -23.97152, 18.18461, -6.89715, 1.04638},
            {-94.86896, 263.99299, -275.39228, 127.638, -22.17592, 0.},
            {-26.53327, 80.49568, -97.54006, 59.01673, -17.83125, 2.15238},
            {-115.76303, 245.91454, -195.75858, 69.20755, -9.16815, 0.},
            {-43.87256, 82.08948, -57.53144, 17.9008, -2.08659, 0.},
            {-78.84221, 174.34897, -153.78386, 67.64246, -14.83937, 1.29914}
    };

    public static final double[] BBB_COEFFS_SICE_BOUNDS =
            new double[]{0.4, 0.8, 1.02, 1.25, 1.4, 1.5, 1.8, 2.0, 2.22, 2.5};


    public static final double[][] KAPPA_2_SPLINE_COEFFS = {
            // from fit to https://atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat
            {
                    // 4th order spline fit, 400-525nm:
                    4.913273573521765E-8 ,
                    -5.446912137768008E-7,
                    2.2340617964072435E-6,
                    -4.029955265312061E-6,
                    2.70443404633965E-6  ,
            },
            {
                    // 4th order spline fit, 525-700nm:
                     7.871975805631149E-7  ,
                     -5.5075218680064285E-6,
                     1.4735203245000926E-5 ,
                     -1.7996480996799473E-5,
                     8.533638647997086E-6  ,
            },
            {
                    // 4th order spline fit, 700-1020nm:
                   3.20662332955914E-4   ,
                    -0.0016037321152243488,
                    0.0030016823137293783,
                    -0.0024933767575023914,
                    7.763358059437843E-4,
            }
    };

    // constants for planar BB albedo retrieval (AK, 20171201):
//    public static final double BB_WVL_0 = 0.3;
    public static final double BB_WVL_1 = 0.3;
    public static final double BB_WVL_2 = 0.7;
    public static final double BB_WVL_3 = 2.4;

    public static final double RHO_ICE = 917.0;  // in kg*m^-3

//    public static final double BARE_ICE_THRESH = 1./3.;   // AK, 20181001

    // define here, we do not want a dependency to Idepix just for this:
    public static final String IDEPIX_CLASSIF_BAND_NAME = "pixel_classif_flags";

    public static final int S3_SNOW_SZA_HIGH = 0;
    public static final int S3_SNOW_GLINT = 1;
    public static final int S3_SNOW_BACKSCATTERING = 2;

    public static final String S3_SNOW_SZA_HIGH_DESCR_TEXT =
            "SZA high (> 75 deg), retrievals may be doubtful";
    public static final String S3_SNOW_GLINT_DESCR_TEXT =
            "Pixel is in glint region (relative azimuth close to 180 deg)";
    public static final String S3_SNOW_BACKSCATTERING_DESCR_TEXT =
            "Pixel is in backscattering region (relative azimuth close to 0 deg)";
}

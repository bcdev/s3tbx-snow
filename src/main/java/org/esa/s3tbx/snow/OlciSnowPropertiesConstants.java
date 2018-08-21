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
    public static final double BB_WVL_1 = 0.3;
    public static final double BB_WVL_2 = 0.7;
    public static final double BB_WVL_3 = 2.4;

    public static final double RHO_ICE = 917.0;  // in kg*m^-3

    // define here, we do not want a dependency to Idepix just for this:
    public static final String IDEPIX_CLASSIF_BAND_NAME = "pixel_classif_flags";
}

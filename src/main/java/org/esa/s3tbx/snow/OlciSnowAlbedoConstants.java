package org.esa.s3tbx.snow;

/**
 * Constants for OLCI snow albedo retrieval
 *
 * @author olafd
 */
public class OlciSnowAlbedoConstants {

//    public static final double[] WAVELENGTH_GRID_OLCI_PARTIAL = {
//            0.299,         // 'outside' wavelength needed for extrapolation
//            0.400,         // OA01
//            0.4125,        // OA02
//            0.4425,        // OA03
//            0.490,         // OA04
//            0.510,         // OA05
//            0.560,         // OA06
//            0.620,         // OA07
//            0.665,         // OA08
//            0.68125,       // OA10
//            0.70875,       // OA11
//            0.77875,       // OA16
//            0.865,         // OA17
//            0.885,         // OA18
//            1.02           // OA21
//    };

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

    public static final double[] WAVELENGTH_GRID_OLCI_EXTENDED = {
            0.299, 0.308, 0.318, 0.328, 0.339, 0.351, 0.364, 0.374, 0.381, 0.388,
            0.396, 0.404, 0.412, 0.421, 0.43, 0.44, 0.449, 0.46, 0.471, 0.482,
            0.494, 0.506, 0.52, 0.533, 0.548, 0.563, 0.58, 0.597, 0.616, 0.635,
            0.652, 0.665, 0.678,  0.692, 0.707, 0.722, 0.738, 0.755, 0.772, 0.791,
            0.81, 0.83, 0.851, 0.874, 0.897, 0.915, 0.927, 0.939, 0.952, 0.964,
            0.985, 1.011, 1.02
    };

    public static final double[] F_LAMBDA_EXTENDED = {
            3.43E-6, 2.905E-4, 0.00243, 0.00607, 0.00795, 0.00906, 0.00878, 0.00697, 0.00552, 0.00585,
            0.00759, 0.01035, 0.01172, 0.01255, 0.01235, 0.01445, 0.01694, 0.01933, 0.019, 0.02038,
            0.02076, 0.022, 0.0215, 0.02336, 0.02455, 0.02499, 0.02612, 0.02695, 0.02783, 0.02592,
            0.02094, 0.01858, 0.01898, 0.01772, 0.01936, 0.0165, 0.01957, 0.02027, 0.0207, 0.02024,
            0.0196, 0.01799, 0.02066, 0.02134, 0.01408, 0.00884, 0.00752, 0.00324, 0.00383, 0.00691,
            0.01564, 0.01451, 0.00884
    };

    public static final int[] SPECTRAL_ALBEDO_OUTPUT_WAVELENGTHS = {
            400, 560, 665, 865, 1020
    };
    public static final int[] SPECTRAL_ALBEDO_OUTPUT_WAVELENGTH_INDICES = {
            0, 5, 7, 16, 20
    };

    // Table 1 from spectral_albedo_october_5_2017.doc:
    public static final double[] LAMBDA_0 = {0.35, 0.525, 0.7, 1.98};
    public static final double[] LAMBDA_MAX = {0.525, 0.7, 1.98, 2.4475};
    public static final double[] H = {0.175, 0.175, 1.28, 0.4675};
    public static final double[] c0 = {3.85E-9, 2.2E-9, 1.30689, 1.27525};
    public static final double[] c1 = {-5.287E-9, 3.384E-9, -0.026, -0.026};
    public static final double[] c2 = {2.2594E-9, 6.41E-9, 0.017956, -0.00149};
    public static final double[] c3 = {1.38411E-9, 1.63665E-8, -0.024, -0.0128};

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
}

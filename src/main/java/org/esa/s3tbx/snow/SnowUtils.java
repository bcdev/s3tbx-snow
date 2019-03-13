package org.esa.s3tbx.snow;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.esa.snap.core.datamodel.FlagCoding;
import org.esa.snap.core.datamodel.Mask;
import org.esa.snap.core.datamodel.Product;
import org.esa.snap.core.util.BitSetter;
import org.esa.snap.core.util.math.MathUtils;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 * S3 Snow utility class
 *
 * @author olafd
 */
public class SnowUtils {

    public static double computeU(double mu) {
        return 3.0 * (1.0 + 2.0 * mu) / 7.0;
    }

    public static RefractiveIndexTable getRefractiveIndexInterpolated(RefractiveIndexTable refractiveIndexTable,
                                                                      SolarSpectrumExtendedTable solarSpectrumExtendedTable) {
        double[] wvlsFull = solarSpectrumExtendedTable.getWvl();

        RefractiveIndexTable refractiveIndexTableInterpolated = new RefractiveIndexTable(wvlsFull.length);
        refractiveIndexTableInterpolated.setWvl(wvlsFull);
        refractiveIndexTableInterpolated.setRefractiveIndexImag(wvlsFull);

        final double[] wvlsRefrIndexAux = refractiveIndexTable.getWvl();
        refractiveIndexTableInterpolated.setRefractiveIndexImag(linearInterpolateWithSplineFunction(wvlsRefrIndexAux,
                                                                                  refractiveIndexTable.getRefractiveIndexImag(),
                                                                                  wvlsFull));
        return refractiveIndexTableInterpolated;
    }

    public static double[] computeFLambda(SolarSpectrumExtendedTable solarSpectrumTable, double sza) {

        final double[][] solarSpectrum = solarSpectrumTable.getSolarSpectrum();

//        final double[] solarSpectrumInterpolated = getSolarSpectrumInterpolated(solarSpectrum, sza);
        final int szaIndex = (int) Math.min(Math.round(sza), 88);
        final double[] solarSpectrumInterpolated = solarSpectrum[szaIndex];

        double[] fLambda = new double[solarSpectrumInterpolated.length];

        for (int i = 0; i < fLambda.length; i++) {
            final double trans = 1.0;   // for the new solar spectrum
            fLambda[i] = trans * solarSpectrumInterpolated[i] / 1000;
        }

        return fLambda;
    }

    public static double[] linearInterpolateWithSplineFunction(double[] x, double[] y, double[] xi) {
        final LinearInterpolator linearInterpolator = new LinearInterpolator();
        PolynomialSplineFunction psf = linearInterpolator.interpolate(x, y);

        double[] yi = new double[xi.length];
        for (int i = 0; i < xi.length; i++) {
            yi[i] = psf.value(xi[i]);
        }
        return yi;
    }

    public static double[] linearInterpolate(double[] x, double[] xa, double[] ya) {
        double[] result = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            result[i] = linearInterpolate(x[i], xa, ya);
        }
        return result;
    }

    public static double linearInterpolate(double x, double[] xa, double[] ya) {
        int kk = -1;
        for (int k = 0; k < xa.length - 1; k++) {
            if (x >= xa[k] && x <= xa[k + 1]) {
                kk = k;
                break;
            }
        }

        if (kk < 0) {
            throw new OutOfRangeException(x, xa[0], xa[xa.length - 1]);
        }

        final double x0 = xa[kk];
        final double x1 = xa[kk + 1];
        final double y0 = ya[kk];
        final double y1 = ya[kk + 1];

        return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    }

    public static String[] setupRcSourceBands(String[] requiredRadianceBandNamesAlbedo, String[] requiredRadianceBandNamesPpa) {
        ArrayList<String> rcSourceBands = new ArrayList<>();
        Collections.addAll(rcSourceBands, requiredRadianceBandNamesAlbedo);
        if (requiredRadianceBandNamesPpa != null) {
            for (String bandName : requiredRadianceBandNamesPpa) {
                if (!rcSourceBands.contains(bandName)) {
                    rcSourceBands.add(bandName);
                }
            }
        }
        return rcSourceBands.toArray(new String[rcSourceBands.size()]);
    }

    public static double cutTo4DecimalPlaces(double value) {
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            return Double.NaN;
        }
        final double x1 = value * 10000.0;
        final double x2 = Math.round(x1);
        return x2 / 10000.0;
    }

    public static double cutTo7DecimalPlaces(double value) {
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            return Double.NaN;
        }
        final double x1 = value * 10000000.0;
        final double x2 = Math.round(x1);
        return x2 / 10000000.0;
    }

    public static double getRelAzi(double saa, double vaa) {
        final double saaRad = Math.toRadians(saa);
        final double vaaRad = Math.toRadians(vaa);
        return Math.toDegrees(Math.acos(Math.cos(saaRad) * Math.cos(vaaRad) + Math.sin(saaRad) * Math.sin(vaaRad)));
    }

    public static double getRelAziSice(double saa, double vaa) {
        // this is the definition as in Fortran breadboard for SICE
        return Math.abs(180.0 - (vaa - saa));
    }

    public static double calcScatteringCos(double sza, double vza, double raa) {
        final double sins = (float) Math.sin(sza * MathUtils.DTOR);
        final double sinv = (float) Math.sin(vza * MathUtils.DTOR);
        final double coss = (float) Math.cos(sza * MathUtils.DTOR);
        final double cosv = (float) Math.cos(vza * MathUtils.DTOR);

        // Compute the geometric conditions
        final double cosphi = Math.cos((180. - raa) * MathUtils.DTOR);  // AK, manual_31_01_2018.docx: use 180 - raa !

        // cos of scattering angle
        return -coss * cosv - sins * sinv * cosphi;
    }

    static double getExtrapolFlux(double lowerFlux, double sza) {
        final double cosSza75 = Math.cos(75. * MathUtils.DTOR);
        final double cosSza = Math.cos(sza * MathUtils.DTOR);
        final double frac = (cosSza75 - cosSza) / cosSza75;  // sza = 90deg --> cos(sza) = 0
        final double upperFlux = 0.0; // sza = 90deg
        return lowerFlux + frac * (upperFlux - lowerFlux);
    }

    static double getInterpolFlux(double[][] solarSpectrum, double sza, int lowerIndex, int upperIndex, int i) {
        final double cosLowerSza = Math.cos(15. * lowerIndex * MathUtils.DTOR);
        final double cosUpperSza = Math.cos(15. * upperIndex * MathUtils.DTOR);
        final double cosSza = Math.cos(sza * MathUtils.DTOR);
        final double frac = (cosSza - cosLowerSza) / (cosUpperSza - cosLowerSza);
        final double lowerFlux = solarSpectrum[lowerIndex][i];
        final double upperFlux = solarSpectrum[upperIndex][i];
        return lowerFlux + frac * (upperFlux - lowerFlux);
    }

    // currently not used
//    private static double[] getSolarSpectrumInterpolated(double[][] solarSpectrum, double sza) {
//        int lowerIndex = (int) sza / 15;   // 0, 1, 2, 3, 4, 5
//        int upperIndex = lowerIndex + 1;
//
//        double[] solarSpectrumInterpolated = new double[solarSpectrum[0].length];
//        for (int i = 0; i < solarSpectrumInterpolated.length; i++) {
//            if (upperIndex <= 5) {
//                final double interpolFlux = getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, i);
//                solarSpectrumInterpolated[i] = interpolFlux;
//            } else {
//                final double extrapolFlux = getExtrapolFlux(solarSpectrum[5][i], sza);
//                solarSpectrumInterpolated[i] = extrapolFlux;
//            }
//        }
//
//        return solarSpectrumInterpolated;
//    }

    /**
     * Provides S3 Snow flag coding
     *
     * @param flagId - the flag ID
     * @return - the flag coding
     */
    public static FlagCoding createS3SnowFlagCoding(String flagId) {
        FlagCoding flagCoding = new FlagCoding(flagId);

        flagCoding.addFlag("S3_SNOW_SZA_HIGH", BitSetter.setFlag(0, OlciSnowPropertiesConstants.S3_SNOW_SZA_HIGH),
                           OlciSnowPropertiesConstants.S3_SNOW_SZA_HIGH_DESCR_TEXT);
        flagCoding.addFlag("S3_SNOW_GLINT", BitSetter.setFlag(0, OlciSnowPropertiesConstants.S3_SNOW_GLINT),
                           OlciSnowPropertiesConstants.S3_SNOW_GLINT_DESCR_TEXT);
        flagCoding.addFlag("S3_SNOW_SCATTERING_ANGLE", BitSetter.setFlag(0, OlciSnowPropertiesConstants.S3_SNOW_BACKSCATTERING),
                           OlciSnowPropertiesConstants.S3_SNOW_BACKSCATTERING_DESCR_TEXT);

        return flagCoding;
    }

    /**
     * Provides OLCI pixel classification flag bitmask
     *
     * @param s3snowProduct - the S3 Snow product
     */
    public static void setupS3SnowBitmask(Product s3snowProduct) {

        int index = 0;
        int w = s3snowProduct.getSceneRasterWidth();
        int h = s3snowProduct.getSceneRasterHeight();
        Mask mask;

        mask = Mask.BandMathsType.create("S3_SNOW_SZA_HIGH",
                                         OlciSnowPropertiesConstants.S3_SNOW_SZA_HIGH_DESCR_TEXT,
                                         w, h,
                                         "s3snow_flags.S3_SNOW_SZA_HIGH",
                                         Color.red, 0.5f);
        s3snowProduct.getMaskGroup().add(index++, mask);

        mask = Mask.BandMathsType.create("S3_SNOW_GLINT",
                                         OlciSnowPropertiesConstants.S3_SNOW_GLINT_DESCR_TEXT,
                                         w, h,
                                         "s3snow_flags.S3_SNOW_GLINT",
                                         Color.yellow, 0.5f);
        s3snowProduct.getMaskGroup().add(index++, mask);

        mask = Mask.BandMathsType.create("S3_SNOW_SCATTERING_ANGLE",
                                         OlciSnowPropertiesConstants.S3_SNOW_BACKSCATTERING_DESCR_TEXT,
                                         w, h,
                                         "s3snow_flags.S3_SNOW_SCATTERING_ANGLE",
                                         Color.blue, 0.5f);
        s3snowProduct.getMaskGroup().add(index, mask);
    }

    /**
     * Provides SICE pollution type flag coding
     *
     * @param flagId - the flag ID
     * @return - the flag coding
     */
    public static FlagCoding createSicePollutionTypeFlagCoding(String flagId) {
        FlagCoding flagCoding = new FlagCoding(flagId);

        flagCoding.addFlag("SICE_POLLUTION_UNCERTAIN", BitSetter.setFlag(0, OlciSnowPropertiesConstants.SICE_POLLUTION_UNCERTAIN),
                           OlciSnowPropertiesConstants.SICE_POLLUTION_UNCERTAIN_DESCR_TEXT);
        flagCoding.addFlag("SICE_POLLUTION_SOOT", BitSetter.setFlag(0, OlciSnowPropertiesConstants.SICE_POLLUTION_SOOT),
                           OlciSnowPropertiesConstants.SICE_POLLUTION_SOOT_DESCR_TEXT);
        flagCoding.addFlag("SICE_POLLUTION_DUST", BitSetter.setFlag(0, OlciSnowPropertiesConstants.SICE_POLLUTION_DUST),
                           OlciSnowPropertiesConstants.SICE_POLLUTION_DUST_DESCR_TEXT);
        flagCoding.addFlag("SICE_POLLUTION_ALGAE", BitSetter.setFlag(0, OlciSnowPropertiesConstants.SICE_POLLUTION_ALGAE),
                           OlciSnowPropertiesConstants.SICE_POLLUTION_ALGAE_DESCR_TEXT);

        return flagCoding;
    }

    /**
     * Provides SICE pollution type flag bitmask
     *
     * @param siceProduct - the SICE product
     */
    public static void setupSicePollutionTypeBitmask(Product siceProduct) {

        int index = 0;
        int w = siceProduct.getSceneRasterWidth();
        int h = siceProduct.getSceneRasterHeight();
        Mask mask;

        mask = Mask.BandMathsType.create("SICE_POLLUTION_UNCERTAIN",
                                         OlciSnowPropertiesConstants.SICE_POLLUTION_UNCERTAIN_DESCR_TEXT,
                                         w, h,
                                         "sice_pollution_type_flags.SICE_POLLUTION_UNCERTAIN",
                                         Color.red, 0.5f);
        siceProduct.getMaskGroup().add(index++, mask);

        mask = Mask.BandMathsType.create("SICE_POLLUTION_SOOT",
                                         OlciSnowPropertiesConstants.SICE_POLLUTION_SOOT_DESCR_TEXT,
                                         w, h,
                                         "sice_pollution_type_flags.SICE_POLLUTION_SOOT",
                                         Color.yellow, 0.5f);
        siceProduct.getMaskGroup().add(index++, mask);

        mask = Mask.BandMathsType.create("SICE_POLLUTION_DUST",
                                         OlciSnowPropertiesConstants.SICE_POLLUTION_DUST_DESCR_TEXT,
                                         w, h,
                                         "sice_pollution_type_flags.SICE_POLLUTION_DUST",
                                         Color.yellow, 0.5f);
        siceProduct.getMaskGroup().add(index++, mask);

        mask = Mask.BandMathsType.create("SICE_POLLUTION_ALGAE",
                                         OlciSnowPropertiesConstants.SICE_POLLUTION_ALGAE_DESCR_TEXT,
                                         w, h,
                                         "sice_pollution_type_flags.SICE_POLLUTION_ALGAE",
                                         Color.blue, 0.5f);
        siceProduct.getMaskGroup().add(index, mask);
    }


    /**
     * Provides SICE ground type flag coding
     *
     * @param flagId - the flag ID
     * @return - the flag coding
     */
    public static FlagCoding createSiceGroundTypeFlagCoding(String flagId) {
        FlagCoding flagCoding = new FlagCoding(flagId);

        flagCoding.addFlag("SICE_UNCERTAIN", BitSetter.setFlag(0, OlciSnowPropertiesConstants.SICE_UNCERTAIN),
                           OlciSnowPropertiesConstants.SICE_UNCERTAIN_DESCR_TEXT);
        flagCoding.addFlag("SICE_SNOW", BitSetter.setFlag(0, OlciSnowPropertiesConstants.SICE_SNOW),
                           OlciSnowPropertiesConstants.SICE_SNOW_DESCR_TEXT);
        flagCoding.addFlag("SICE_BARE_ICE_CLEAN", BitSetter.setFlag(0, OlciSnowPropertiesConstants.SICE_BARE_ICE_CLEAN),
                           OlciSnowPropertiesConstants.SICE_BARE_ICE_CLEAN_DESCR_TEXT);
        flagCoding.addFlag("SICE_BARE_ICE_POLLUTED", BitSetter.setFlag(0, OlciSnowPropertiesConstants.SICE_BARE_ICE_POLLUTED),
                           OlciSnowPropertiesConstants.SICE_BARE_ICE_POLLUTED_DESCR_TEXT);

        return flagCoding;
    }

    /**
     * Provides SICE ground type flag bitmask
     *
     * @param siceProduct - the SICE product
     */
    public static void setupSiceGroundTypeBitmask(Product siceProduct) {

        int index = 0;
        int w = siceProduct.getSceneRasterWidth();
        int h = siceProduct.getSceneRasterHeight();
        Mask mask;

        mask = Mask.BandMathsType.create("SICE_UNCERTAIN",
                                         OlciSnowPropertiesConstants.SICE_UNCERTAIN_DESCR_TEXT,
                                         w, h,
                                         "sice_ground_type_flags.SICE_UNCERTAIN",
                                         Color.red, 0.5f);
        siceProduct.getMaskGroup().add(index++, mask);

        mask = Mask.BandMathsType.create("SICE_SNOW",
                                         OlciSnowPropertiesConstants.SICE_SNOW_DESCR_TEXT,
                                         w, h,
                                         "sice_ground_type_flags.SICE_SNOW",
                                         Color.yellow, 0.5f);
        siceProduct.getMaskGroup().add(index++, mask);

        mask = Mask.BandMathsType.create("SICE_BARE_ICE_CLEAN",
                                         OlciSnowPropertiesConstants.SICE_BARE_ICE_CLEAN_DESCR_TEXT,
                                         w, h,
                                         "sice_ground_type_flags.SICE_BARE_ICE_CLEAN",
                                         Color.yellow, 0.5f);
        siceProduct.getMaskGroup().add(index++, mask);

        mask = Mask.BandMathsType.create("SICE_BARE_ICE_POLLUTED",
                                         OlciSnowPropertiesConstants.SICE_BARE_ICE_POLLUTED_DESCR_TEXT,
                                         w, h,
                                         "sice_ground_type_flags.SICE_BARE_ICE_POLLUTED",
                                         Color.blue, 0.5f);
        siceProduct.getMaskGroup().add(index, mask);
    }


}

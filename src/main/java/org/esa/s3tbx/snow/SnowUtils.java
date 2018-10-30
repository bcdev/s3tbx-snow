package org.esa.s3tbx.snow;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.esa.snap.core.util.math.MathUtils;

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
        refractiveIndexTableInterpolated.setRefractiveIndexImag(splineInterpolate(wvlsRefrIndexAux,
                                                                                  refractiveIndexTable.getRefractiveIndexImag(),
                                                                                  wvlsFull));
        return refractiveIndexTableInterpolated;
    }

    public static double[] computeFLambda(SolarSpectrumExtendedTable solarSpectrumTable, double sza) {

        final double[][] solarSpectrum = solarSpectrumTable.getSolarSpectrum();

        final double[] solarSpectrumInterpolated = getSolarSpectrumInterpolated(solarSpectrum, sza);
//        final int szaIndex = (int) Math.min(Math.round(sza), 88);
//        final double[] solarSpectrumInterpolated = solarSpectrum[szaIndex];

        double[] fLambda = new double[solarSpectrumInterpolated.length];

        for (int i = 0; i < fLambda.length; i++) {
            final double trans = 1.0;   // for the new solar spectrum
            fLambda[i] = trans * solarSpectrumInterpolated[i] / 1000;
        }

        return fLambda;
    }

    public static double[] splineInterpolate(double[] x, double[] y, double[] xi) {
        final LinearInterpolator linearInterpolator = new LinearInterpolator();
        PolynomialSplineFunction psf = linearInterpolator.interpolate(x, y);

        double[] yi = new double[xi.length];
        for (int i = 0; i < xi.length; i++) {
            yi[i] = psf.value(xi[i]);
        }
        return yi;
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
        return x2/10000.0;
    }

    public static double cutTo7DecimalPlaces(double value) {
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            return Double.NaN;
        }
        final double x1 = value * 10000000.0;
        final double x2 = Math.round(x1);
        return x2/10000000.0;
    }

    public static double getRelAzi(double saa, double vaa) {
        final double saaRad = Math.toRadians(saa);
        final double vaaRad = Math.toRadians(vaa);
        return Math.toDegrees(Math.acos(Math.cos(saaRad) * Math.cos(vaaRad) + Math.sin(saaRad) * Math.sin(vaaRad)));
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

    private static double[] getSolarSpectrumInterpolated(double[][] solarSpectrum, double sza) {
        int lowerIndex = (int) sza/15;   // 0, 1, 2, 3, 4, 5
        int upperIndex = lowerIndex + 1;

        double[] solarSpectrumInterpolated = new double[solarSpectrum[0].length];
        for (int i = 0; i < solarSpectrumInterpolated.length; i++) {
            if (upperIndex <= 5) {
                final double interpolFlux = getInterpolFlux(solarSpectrum, sza, lowerIndex, upperIndex, i);
                solarSpectrumInterpolated[i] = interpolFlux;
            } else {
                final double extrapolFlux = getExtrapolFlux(solarSpectrum[5][i], sza);
                solarSpectrumInterpolated[i] = extrapolFlux;
            }
        }

        return solarSpectrumInterpolated;
    }

}

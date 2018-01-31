package org.esa.s3tbx.snow;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.ArrayList;
import java.util.Collections;

import static org.esa.s3tbx.snow.OlciSnowAlbedoConstants.*;

/**
 * S3 Snow utility class
 *
 * @author olafd
 */
public class SnowUtils {

    public static double computeKappa2(double wvl, int rangeIndex) {
        if (rangeIndex < 0 || rangeIndex > 3) {
            throw new IllegalArgumentException("rangeIndex must be in [0,3]");
        }
        return c0[rangeIndex] + c1[rangeIndex] * (wvl - LAMBDA_0[rangeIndex]) / H[rangeIndex] +
                c2[rangeIndex] * Math.pow((wvl - LAMBDA_0[rangeIndex]) / H[rangeIndex], 2.0) +
                c3[rangeIndex] * Math.pow((wvl - LAMBDA_0[rangeIndex]) / H[rangeIndex], 3.0);
    }

    public static double computeU(double mu) {
        return 3.0 * (1.0 + 2.0 * mu) / 7.0;
    }

    public static RefractiveIndexTable getRefractiveIndexInterpolated(RefractiveIndexTable refractiveIndexTable,
                                                                      SolarSpectrumTable solarSpectrumTable) {
        double[] wvlsFull = solarSpectrumTable.getWvl();

        RefractiveIndexTable refractiveIndexTableInterpolated = new RefractiveIndexTable(wvlsFull.length);
        refractiveIndexTableInterpolated.setWvl(wvlsFull);
        refractiveIndexTableInterpolated.setRefractiveIndexImag(wvlsFull);

        final double[] wvlsRefrIndexAux = refractiveIndexTable.getWvl();
        refractiveIndexTableInterpolated.setRefractiveIndexImag(splineInterpolate(wvlsRefrIndexAux,
                                                                                  refractiveIndexTable.getRefractiveIndexImag(),
                                                                                  wvlsFull));
        return refractiveIndexTableInterpolated;
    }

    public static double getRefractiveIndex(RefractiveIndexTable refractiveIndexTable,
                                            double wvl) {
        final double[] wvlsTable = refractiveIndexTable.getWvl();
        for (int i = 0; i < wvlsTable.length - 1; i++) {
            if (wvl >= wvlsTable[i] && wvl < wvlsTable[i + 1]) {
                return refractiveIndexTable.getRefractiveIndexImag(i);
            }
        }
        return -1;
    }

    public static double[] getFLambda(SolarSpectrumTable solarSpectrumTable) {

        final double[] solarSpectrum = solarSpectrumTable.getSolarSpectrum();

        double[] fLambda = new double[solarSpectrum.length];

        for (int i = 0; i < fLambda.length; i++) {
            final double trans = 1.0;   // for the new solar spectrum
            fLambda[i] = trans * solarSpectrum[i] / 1000;
        }

        return fLambda;
    }

    public static double getFLambda(SolarSpectrumTable solarSpectrumTable, double mu0, double wvl) {
        final double[] fLambda = getFLambda(solarSpectrumTable);
        final double[] wvlsTable = solarSpectrumTable.getWvl();
        for (int i = 0; i < wvlsTable.length - 1; i++) {
            if (wvl >= wvlsTable[i] && wvl < wvlsTable[i + 1]) {
                return fLambda[i];
            }
        }
        return -1;
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

}

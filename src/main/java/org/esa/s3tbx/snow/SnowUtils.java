package org.esa.s3tbx.snow;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import static org.esa.s3tbx.snow.OlciSnowAlbedoConstants.*;

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

    public static double[] computeKappa2() {
        double[] kappa2 = new double[WAVELENGTH_GRID_OLCI.length];

        // first interval, bands 1-5
        for (int i = 0; i < 5; i++) {
            final double wvl = WAVELENGTH_GRID_OLCI[i];
            kappa2[i] = computeKappa2(wvl, 0);
        }
        // second interval, bands 6-9
        for (int i = 6; i < 9; i++) {
            final double wvl = WAVELENGTH_GRID_OLCI[i];
            kappa2[i] = computeKappa2(wvl, 1);
        }
        // third interval, bands 10-21
        for (int i = 10; i < WAVELENGTH_GRID_OLCI.length; i++) {
            final double wvl = WAVELENGTH_GRID_OLCI[i];
            kappa2[i] = computeKappa2(wvl, 2);
        }

        return kappa2;
    }

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

//    public static double[] getFLambda(SolarSpectrumTable solarSpectrumTable, double mu0) {
//
//        final double[] wvl = solarSpectrumTable.getWvl();
//        final double[] solarSpectrum = solarSpectrumTable.getSolarSpectrum();
//
//        double[] fLambda = new double[wvl.length];
//
//        for (int i = 0; i < fLambda.length; i++) {
//            final double p1 = OlciSnowAlbedoConstants.TAU_MOL_P1;
//            final double p2 = OlciSnowAlbedoConstants.TAU_MOL_P2;
//            final double p3 = OlciSnowAlbedoConstants.TAU_MOL_P3;
//            final double beta = OlciSnowAlbedoConstants.BETA;
//            final double g = OlciSnowAlbedoConstants.G;
//            final double lambda0 = OlciSnowAlbedoConstants.LAMBDA0;
//            final double f = (1.0 + g) / 2.0;
//            final double kappa = 1.0 - f;
//            final double tauMol =
//                    (p1 / Math.pow(wvl[i], 4.0)) * (1.0 + (p2 / Math.pow(wvl[i], 2.0) + (p3 / Math.pow(wvl[i], 4.0))));
//            final double tauA = beta * lambda0 / wvl[i];
//            final double trans = Math.exp(-(0.5 * tauMol + kappa * tauA) / mu0);
//            fLambda[i] = trans * solarSpectrum[i] / 1000;
//        }
//
//        return fLambda;
//    }

    public static double[] getFLambda(SolarSpectrumTable solarSpectrumTable, double mu0) {

//        final double[] wvl = solarSpectrumTable.getWvl();
        final double[] solarSpectrum = solarSpectrumTable.getSolarSpectrum();

        double[] fLambda = new double[solarSpectrum.length];

        for (int i = 0; i < fLambda.length; i++) {
            final double trans = 1.0;   // for the new solar spectrum
            fLambda[i] = trans * solarSpectrum[i] / 1000;
        }

        return fLambda;
    }

    public static double getFLambda(SolarSpectrumTable solarSpectrumTable, double mu0, double wvl) {
        final double[] fLambda = getFLambda(solarSpectrumTable, mu0);
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


}

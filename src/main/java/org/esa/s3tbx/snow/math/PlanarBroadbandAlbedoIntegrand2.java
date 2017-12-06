package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.esa.s3tbx.snow.RefractiveIndexTable;
import org.esa.s3tbx.snow.SnowUtils;
import org.esa.s3tbx.snow.SolarSpectrumTable;

/**
 *
 *
 * Created by Olaf on 03.10.2017.
 */
public class PlanarBroadbandAlbedoIntegrand2 implements UnivariateFunction {

    public static final int NUMERATOR = 1;
    public static final int DENOMINATOR = 2;

    private RefractiveIndexTable refractiveIndexInterpolatedTable;
    private SolarSpectrumTable solarSpectrumTable;
    private double mu_0;
    private double d;
    private int type;

    public PlanarBroadbandAlbedoIntegrand2(RefractiveIndexTable refractiveIndexInterpolatedTable,
                                    SolarSpectrumTable solarSpectrumTable,
                                    double mu_0, double d, int type) {
        this.refractiveIndexInterpolatedTable = refractiveIndexInterpolatedTable;
        this.solarSpectrumTable = solarSpectrumTable;
        this.mu_0 = mu_0;
        this.d = d;

        if (type != NUMERATOR && type != DENOMINATOR) {
            throw new IllegalArgumentException("Integrand must be the NUMERATOR (1) or DENOMINATOR (2)");
        }
        this.type = type;
    }

    public double value(double x) {
        double u;
        if (mu_0 == 1.0) {
            u = 1.0;
        } else {
            u = SnowUtils.computeU(mu_0);
        }
        final double chi = SnowUtils.getRefractiveIndex(refractiveIndexInterpolatedTable, x);
        final double k = 4.0*Math.PI*chi/x;

         final double planarSpectralAlbedo = Math.exp(-3.6*u*Math.sqrt(k*d));
         final double fLambda = SnowUtils.getFLambda(solarSpectrumTable, mu_0, x);

        if (type == NUMERATOR) {
//            System.out.println("wvl, planarSpectralAlbedo*flambda = " + x + ", " + planarSpectralAlbedo * fLambda);
            return planarSpectralAlbedo * fLambda;
        } else {
//            System.out.println("wvl, flambda, planarSpectralAlbedo = " +
//                                       x + ", " + fLambda + ", " + planarSpectralAlbedo);
            return fLambda;
        }
    }
}

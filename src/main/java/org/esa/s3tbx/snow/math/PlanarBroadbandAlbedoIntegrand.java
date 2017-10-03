package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.UnivariateFunction;

/**
 *
 *
 * Created by Olaf on 03.10.2017.
 */
class PlanarBroadbandAlbedoIntegrand implements UnivariateFunction {

    static final int NUMERATOR = 1;
    static final int DENOMINATOR = 2;

    private double mu_0;
    private double p;
    private int type;

    PlanarBroadbandAlbedoIntegrand(double mu_0, double p, int type) {
        this.mu_0 = mu_0;
        this.p = p;

        if (type != NUMERATOR && type != DENOMINATOR) {
            throw new IllegalArgumentException("Integrand must be the NUMERATOR (1) or DENOMINATOR (2)");
        }
        this.type = type;
    }

    public double value(double x) {
        final double a = -2.3;
        final double b = 3.3;
        final double lambda_0 = 1326.6E-3;
        final double d_lambda = 125.3E-3;
        final double xi_0 = 32.38;
        final double xi_1 = -160140.33;
        final double xi_2 = 7959.53;
        final double gamma_1 = 85.34E-3;
        final double gamma_2 = 401.79E-3;
        final double eta = 0.00877;
        final double nu = 4.05;
        final double beta = 1.3;
        final double f = 0.85;
        final double tau_a = 0.1;
        final double p_0 = 1013.25;

        final double factor2 = xi_0 + xi_1 * Math.exp(-x / gamma_1) + xi_2 * Math.exp(-x / gamma_2);
        final double factor3 = Math.exp(-(0.5 * (p / p_0) * eta * Math.pow(x, -nu) +
                (1.0 - f) * tau_a * Math.pow(lambda_0 / x, beta))) / mu_0;

        if (type == NUMERATOR) {
            final double factor1 = a + b / (1.0 + Math.exp((x - lambda_0) / d_lambda));
            return factor1 * factor2 * factor3;
        } else {
            return factor2 * factor3;
        }
    }
}

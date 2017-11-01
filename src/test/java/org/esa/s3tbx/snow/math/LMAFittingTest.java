package org.esa.s3tbx.snow.math;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.esa.s3tbx.snow.math.lma.LMA;
import org.esa.s3tbx.snow.math.lma.LMAFunction;
import org.esa.s3tbx.snow.math.lma.implementations.Polynomial;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

public class LMAFittingTest {

    @Test
    public void testPolynominalFit() throws Exception {

//        curveFitter.addObservedPoint(4.000E-001,  2.365E-011);
//        curveFitter.addObservedPoint(4.100E-001,  2.669E-011);
//        curveFitter.addObservedPoint(4.200E-001,  3.135E-011);
//        curveFitter.addObservedPoint(4.300E-001,  4.140E-011);
//        curveFitter.addObservedPoint(4.400E-001,  6.268E-011);
//        curveFitter.addObservedPoint(4.500E-001,  9.239E-011);
//        curveFitter.addObservedPoint(4.600E-001,  1.325E-010);
//        curveFitter.addObservedPoint(4.700E-001,  1.956E-010);
//        curveFitter.addObservedPoint(4.800E-001,  2.861E-010);
//        curveFitter.addObservedPoint(4.900E-001,  4.172E-010);
//        curveFitter.addObservedPoint(5.000E-001,  5.889E-010);
//        curveFitter.addObservedPoint(5.100E-001,  8.036E-010);
//        curveFitter.addObservedPoint(5.200E-001,  1.076E-009);

        LMA lma = new LMA(
                new Polynomial(),
                new double[] {1, 1,1., 1., 1.},
                new double[][] {
                        {0.4, .41, .42, .43, .44, .45, .46, .47, .48, .49, .5, .51, .52},
                        {2.365E-011, 2.669E-011, 3.135E-011, 4.140E-011, 6.268E-011, 9.239E-011,
                                1.325E-010, 1.956E-010, 2.861E-010, 4.172E-010, 5.889E-010, 8.036E-010, 1.076E-009}}
        );
        lma.fit();

        System.out.println("POLY RESULT PARAMETERS: " + Arrays.toString(lma.parameters));
    }

    @Test
    public void testSinFit() throws Exception {
        double[] x = {0.0, 0.1, 0.2, 0.3, 0.5, 0.7};//, 1.1, 1.4, 2.5, 6.4, 7.9, 10.4, 12.6};
        double[] a = {2.2, 0.4};
        double[][] data = {x, sin.generateData(x, a)};
        LMA lma = new LMA(sin, new double[] {0.1, 10}, data, null);
        lma.fit();
        System.out.println("SIN RESULT PARAMETERS: " + Arrays.toString(lma.parameters));
    }

    private static LMAFunction sin = new LMAFunction() {
        public double getY(double x, double[] a) {
            return a[0] * Math.sin(x / a[1]);
        }
        public double getPartialDerivate(double x, double[] a, int parameterIndex) {
            switch (parameterIndex) {
                case 0: return Math.sin(x / a[1]);
                case 1: return a[0] * Math.cos(x / a[1]) * (-x / (a[1] * a[1]));
            }
            throw new RuntimeException("No such fit parameter: " + parameterIndex);
        }
    };
}

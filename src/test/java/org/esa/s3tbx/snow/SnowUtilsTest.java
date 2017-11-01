package org.esa.s3tbx.snow;

import org.junit.Ignore;
import org.junit.Test;

import static junit.framework.Assert.assertNotNull;
import static junit.framework.TestCase.assertEquals;


public class SnowUtilsTest {

    @Test
    @Ignore
    public void testComputeKappa2() throws Exception {

        // todo: do own fit to 3rd order polynoms, compare with Table 1 which seems to be wrong (i.e third, fourth row)
        //
        double expected4 = 3.84987E-9 - 5.287E-9 + 2.2594E-9 + 1.38411E-9;
        System.out.println("expected4 = " + expected4);
        assertEquals(expected4, SnowUtils.computeKappa2(0.525, 0), 1.E-11);
        double expected5 = 3.84987E-9 - 0.5*5.287E-9 + 0.25*2.2594E-9 + 0.125*1.38411E-9;
        System.out.println("expected5 = " + expected5);
        assertEquals(expected5, SnowUtils.computeKappa2(0.4375, 0), 1.E-11);

        // results should be c_0[row] at lambda_0[row] from Table 1:
        double expected0 = 3.84987E-9;
        assertEquals(expected0, SnowUtils.computeKappa2(0.35, 0), 1.E-11);
        double expected1 = 2.20653E-9;
        assertEquals(expected1, SnowUtils.computeKappa2(0.525, 0), 1.E-11);
        assertEquals(expected1, SnowUtils.computeKappa2(0.525, 1), 1.E-11);
        double expected2 = 1.30689;
        assertEquals(expected2, SnowUtils.computeKappa2(0.7, 1), 1.E-12);
        assertEquals(expected2, SnowUtils.computeKappa2(0.7, 2), 1.E-12);
        double expected3 = 1.27525;
        assertEquals(expected3, SnowUtils.computeKappa2(1.98, 2), 1.E-3);
        assertEquals(expected3, SnowUtils.computeKappa2(1.98, 3), 1.E-3);

        System.out.println("kappa_2 at 400  = " + SnowUtils.computeKappa2(0.400, 0));
        System.out.println("kappa_2 at 490  = " + SnowUtils.computeKappa2(0.49, 0));
        System.out.println("kappa_2 at 865  = " + SnowUtils.computeKappa2(0.865, 0));
        System.out.println("kappa_2 at 1020 = " + SnowUtils.computeKappa2(1.02, 0));

//        final double[] kappa2 = SnowUtils.computeKappa2();
//        assertNotNull(kappa2);
//        assertEquals(21, kappa2.length);
//        assertEquals(2.365E-11, kappa2[0], 1.E-12);

    }

}

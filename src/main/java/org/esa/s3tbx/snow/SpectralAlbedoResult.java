package org.esa.s3tbx.snow;

/**
 * Holder for spectral albedo results
 *
 * @author olafd
 */
public class SpectralAlbedoResult {
    private final double[][] spectralAlbedos;
    private final double r0;
    private final double f;
    private final double l;
    private final double m;

    private final double r0RelErr;
    private final double fRelErr;
    private final double lRelErr;
    private final double mRelErr;

    SpectralAlbedoResult(double[][] spectralAlbedos, double r0, double f, double l, double m,
                         double r0RelErr, double fRelErr, double lRelErr, double mRelErr) {
        this.spectralAlbedos = spectralAlbedos;
        this.r0 = r0;
        this.f = f;
        this.l = l;
        this.m = m;
        this.r0RelErr = r0RelErr;
        this.fRelErr = fRelErr;
        this.lRelErr = lRelErr;
        this.mRelErr = mRelErr;
    }

    public double[][] getSpectralAlbedos() {
        return spectralAlbedos;
    }

    public double getR0() {
        return r0;
    }

    public double getF() {
        return f;
    }

    public double getL() {
        return l;
    }

    public double getM() {
        return m;
    }

    public double getR0RelErr() {
        return r0RelErr;
    }

    public double getfRelErr() {
        return fRelErr;
    }

    public double getlRelErr() {
        return lRelErr;
    }

    public double getmRelErr() {
        return mRelErr;
    }
}

package org.esa.s3tbx.snow;

/**
 * Holder for SICE snow properties:
 * - snow grain size
 * - snow specific area
 * - snow pollution load
 * - relative impurity load
 * - spectral albedo
 *
 * @author olafd
 */
public class SiceSnowPropertiesResult {
    private double effAbsLength;
    private double snowGrainSize;
    private double snowSpecificArea;
    private double relativeImpurityLoad;
    private double r0a1Thresh;
    private double[][] spectralAlbedos;

    SiceSnowPropertiesResult(double effAbsLength,
                             double snowGrainSize,
                             double snowSpecificArea,
                             double relativeImpurityLoad,
                             double r0a1Thresh,
                             double[][] spectralAlbedos) {
        this.effAbsLength = effAbsLength;
        this.snowGrainSize = snowGrainSize;
        this.snowSpecificArea = snowSpecificArea;
        this.relativeImpurityLoad = relativeImpurityLoad;
        this.r0a1Thresh = r0a1Thresh;
        this.spectralAlbedos = spectralAlbedos;
    }

    public void setSpectralAlbedos(double[][] spectralAlbedos) {
        this.spectralAlbedos = spectralAlbedos;
    }

    public double[][] getSpectralAlbedos() {
        return spectralAlbedos;
    }

    public double getR0a1Thresh() {
        return r0a1Thresh;
    }

    public void setR0a1Thresh(double r0a1Thresh) {
        this.r0a1Thresh = r0a1Thresh;
    }

    public double getEffAbsLength() {
        return effAbsLength;
    }

    public double getSnowGrainSize() {
        return snowGrainSize;
    }

    public double getSnowSpecificArea() {
        return snowSpecificArea;
    }

    public double getRelativeImpurityLoad() {
        return relativeImpurityLoad;
    }
}

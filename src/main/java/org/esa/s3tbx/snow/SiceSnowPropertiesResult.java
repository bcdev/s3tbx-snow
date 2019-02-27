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
    private double[][] spectralAlbedos;
    private double effAbsLength;
    private double snowGrainSize;
    private double snowSpecificArea;
    private double relativeImpurityLoad;

    SiceSnowPropertiesResult(double[][] spectralAlbedos,
                             double effAbsLength,
                             double snowGrainSize,
                             double snowSpecificArea,
                             double relativeImpurityLoad) {
        this.spectralAlbedos = spectralAlbedos;
        this.effAbsLength = effAbsLength;
        this.snowGrainSize = snowGrainSize;
        this.snowSpecificArea = snowSpecificArea;
        this.relativeImpurityLoad = relativeImpurityLoad;
    }

    public void setSpectralAlbedos(double[][] spectralAlbedos) {
        this.spectralAlbedos = spectralAlbedos;
    }

    public double[][] getSpectralAlbedos() {
        return spectralAlbedos;
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

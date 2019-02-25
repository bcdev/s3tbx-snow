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
    private final double[][] spectralAlbedos;
    private final double snowGrainSize;
    private final double snowSpecificArea;
    private final double relativeImpurityLoad;

    SiceSnowPropertiesResult(double[][] spectralAlbedos,
                             double snowGrainSize,
                             double snowSpecificArea,
                             double relativeImpurityLoad) {
        this.spectralAlbedos = spectralAlbedos;
        this.snowGrainSize = snowGrainSize;
        this.snowSpecificArea = snowSpecificArea;
        this.relativeImpurityLoad = relativeImpurityLoad;
    }

    public double[][] getSpectralAlbedos() {
        return spectralAlbedos;
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

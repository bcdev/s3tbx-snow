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
    private double[][] broadbandAlbedos;

    SiceSnowPropertiesResult(double effAbsLength,
                             double snowGrainSize,
                             double snowSpecificArea,
                             double relativeImpurityLoad,
                             double r0a1Thresh,
                             double[][] spectralAlbedos,
                             double[][] broadbandAlbedos) {
        this.effAbsLength = effAbsLength;
        this.snowGrainSize = snowGrainSize;
        this.snowSpecificArea = snowSpecificArea;
        this.relativeImpurityLoad = relativeImpurityLoad;
        this.r0a1Thresh = r0a1Thresh;
        this.spectralAlbedos = spectralAlbedos;
        this.broadbandAlbedos = broadbandAlbedos;
    }

    public void setSpectralAlbedos(double[][] spectralAlbedos) {
        this.spectralAlbedos = spectralAlbedos;
    }

    public double[][] getSpectralAlbedos() {
        return spectralAlbedos;
    }

    public void setBroadbandAlbedos(double[][] broadbandAlbedos) {
        this.broadbandAlbedos = broadbandAlbedos;
    }

    public void setPlanarBroadbandAlbedos(double[] planarBroadbandAlbedos) {
        if (this.broadbandAlbedos == null) {
            this.broadbandAlbedos = new double[2][];
        }
        this.broadbandAlbedos[0] = planarBroadbandAlbedos;
    }

    public void setSphericalBroadbandAlbedos(double[] planarBroadbandAlbedos) {
        if (this.broadbandAlbedos == null) {
            this.broadbandAlbedos = new double[2][];
        }
        this.broadbandAlbedos[1] = planarBroadbandAlbedos;
    }

    public double[][] getBroadbandAlbedos() {
        return broadbandAlbedos;
    }

    public double[] getPlanarBroadbandAlbedos() {
        return broadbandAlbedos[0];
    }

    public double[] getSphericalBroadbandAlbedos() {
        return broadbandAlbedos[1];
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

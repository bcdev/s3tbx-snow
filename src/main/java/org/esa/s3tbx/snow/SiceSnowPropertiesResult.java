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
    private double scatteringAngle;
    private double effAbsLength;
    private double snowGrainSize;
    private double snowSpecificArea;
    private SiceSnowImpurity snowImpurity;
    private double r0a1Thresh;
    private double[][] spectralAlbedos;
    private double[][] broadbandAlbedos;

    SiceSnowPropertiesResult(double effAbsLength,
                             double snowGrainSize,
                             double snowSpecificArea,
                             SiceSnowImpurity snowImpurity,
                             double scatteringAngle,
                             double r0a1Thresh,
                             double[][] spectralAlbedos,
                             double[][] broadbandAlbedos) {
        this.effAbsLength = effAbsLength;
        this.snowGrainSize = snowGrainSize;
        this.snowSpecificArea = snowSpecificArea;
        this.snowImpurity = snowImpurity;
        this.scatteringAngle = scatteringAngle;
        this.r0a1Thresh = r0a1Thresh;
        this.spectralAlbedos = spectralAlbedos;
        this.broadbandAlbedos = broadbandAlbedos;
    }

    public void setSphericalSpectralAlbedos(double[] sphericalSpectralAlbedos) {
        if (this.spectralAlbedos == null) {
            this.spectralAlbedos = new double[2][];
        }
        this.spectralAlbedos[0] = sphericalSpectralAlbedos;
    }

    public void setPlanarSpectralAlbedos(double[] planarSpectralAlbedos) {
        if (this.spectralAlbedos == null) {
            this.spectralAlbedos = new double[2][];
        }
        this.spectralAlbedos[1] = planarSpectralAlbedos;
    }

    public void setSphericalBroadbandAlbedos(double[] sphericalBroadbandAlbedos) {
        if (this.broadbandAlbedos == null) {
            this.broadbandAlbedos = new double[2][];
        }
        this.broadbandAlbedos[0] = sphericalBroadbandAlbedos;
    }

    public void setPlanarBroadbandAlbedos(double[] planarBroadbandAlbedos) {
        if (this.broadbandAlbedos == null) {
            this.broadbandAlbedos = new double[2][];
        }
        this.broadbandAlbedos[1] = planarBroadbandAlbedos;
    }

    public double[] getSphericalSpectralAlbedos() {
        return spectralAlbedos[0];
    }

    public double[] getPlanarSpectralAlbedos() {
        return spectralAlbedos[1];
    }

    public double[] getSphericalBroadbandAlbedos() {
        return broadbandAlbedos[0];
    }

    public double[] getPlanarBroadbandAlbedos() {
        return broadbandAlbedos[1];
    }

    public void setScatteringAngle(double scatteringAngle) {
        this.scatteringAngle = scatteringAngle;
    }

    public double getScatteringAngle() {
        return scatteringAngle;
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

    public SiceSnowImpurity getSnowImpurity() {
        return snowImpurity;
    }
}

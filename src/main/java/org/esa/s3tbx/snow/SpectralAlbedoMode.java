package org.esa.s3tbx.snow;

/**
 * Enumeration for mode (AK algorithm) for spectral albedo retrieval.
 *
 * @author olafd
 */
public enum SpectralAlbedoMode {
    SIGMOIDAL("SIGMOIDAL"),
    EXPONENTIAL_SQRT("EXPONENTIAL_SQRT"),
    EXPONENTIAL_3PARAM("EXPONENTIAL_3PARAM");

    private String name;

    SpectralAlbedoMode(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }
}

package org.esa.s3tbx.snow;

/**
 * Enumeration for mode (AK algorithm) for spectral albedo retrieval.
 *
 * @author olafd
 */
public enum SpectralAlbedoMode {
    SIGMOIDAL_FIT("SIGMOIDAL_FIT"),
    EXPONENTIAL_SQRT_FIT("EXPONENTIAL_SQRT_FIT"),
    POLYNOMINAL_FIT("POLYNOMINAL_FIT"),
    EXPONENTIAL_3PARAM_FIT("EXPONENTIAL_3PARAM_FIT");

    private String name;

    SpectralAlbedoMode(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }
}

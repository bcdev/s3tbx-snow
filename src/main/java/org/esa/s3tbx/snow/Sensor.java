package org.esa.s3tbx.snow;

import static org.esa.s3tbx.snow.SensorConstants.*;

/**
 * Enumeration for supported sensors for snow albedo retrieval (just OLCI so far).
 *
 * @author olafd
 */
public enum Sensor {

    OLCI("OLCI", OLCI_SZA_NAME, OLCI_VZA_NAME, OLCI_SAA_NAME, OLCI_VAA_NAME,
         OLCI_L1B_FLAGS_NAME, OLCI_INVALID_BIT, OLCI_LAND_BIT,
         OLCI_DEFAULT_REQUIRED_RADIANCE_BAND_NAMES,
         OLCI_DEFAULT_REQUIRED_REFL_BAND_NAMES,
         OLCI_DEFAULT_REQUIRED_BRR_BAND_NAMES,
         OLCI_TARGET_TPGS);

    private String name;
    private String szaName;
    private String vzaName;
    private String saaName;
    private String vaaName;
    private String l1bFlagsName;
    private int invalidBit;
    private int landBit;
    private String[] requiredRadianceBandNames;
    private String[] requiredReflBandNames;
    private String[] requiredBrrBandNames;
    private String[] targetTpgs;

    Sensor(String name, String szaName, String vzaName, String saaName, String vaaName,
           String l1bFlagsName, int invalidBit, int landBit, String requiredRadianceBandNames[],
           String requiredReflBandNames[], String[] requiredBrrBandNames, String[] targetTpgs) {
        this.name = name;
        this.szaName = szaName;
        this.vzaName = vzaName;
        this.saaName = saaName;
        this.vaaName = vaaName;
        this.l1bFlagsName = l1bFlagsName;
        this.invalidBit = invalidBit;
        this.landBit = landBit;
        this.requiredRadianceBandNames = requiredRadianceBandNames;
        this.requiredReflBandNames = requiredReflBandNames;
        this.requiredBrrBandNames = requiredBrrBandNames;
        this.targetTpgs = targetTpgs;
    }

    public String getName() {
        return name;
    }

    public String getSzaName() {
        return szaName;
    }

    public String getVzaName() {
        return vzaName;
    }

    public String getSaaName() {
        return saaName;
    }

    public String getVaaName() {
        return vaaName;
    }

    public String getL1bFlagsName() {
        return l1bFlagsName;
    }

    public int getInvalidBit() {
        return invalidBit;
    }

    public int getLandBit() {
        return landBit;
    }

    public String[] getRequiredRadianceBandNames() {
        return requiredRadianceBandNames;
    }

    public void setRequiredRadianceBandNames(String[] requiredRadianceBandNames) {
        this.requiredRadianceBandNames = requiredRadianceBandNames;
    }

    public String[] getRequiredReflBandNames() {
        return requiredReflBandNames;
    }

    public String[] getRequiredBrrBandNames() {
        return requiredBrrBandNames;
    }

    public void setRequiredBrrBandNames(String[] requiredBrrBandNames) {
        this.requiredBrrBandNames = requiredBrrBandNames;
    }

    public String[] getTargetTpgs() {
        return targetTpgs;
    }
}

package org.esa.s3tbx.snow;

import static org.esa.s3tbx.snow.SensorConstants.*;

/**
 * Enumeration for supported sensors for snow albedo retrieval (just OLCI so far).
 *
 * @author olafd
 */
public enum Sensor {

    OLCI("OLCI", OLCI_SZA_NAME, OLCI_VZA_NAME, OLCI_L1B_FLAGS_NAME, OLCI_INVALID_BIT, OLCI_TARGET_TPGS);

    private String name;
    private String szaName;
    private String vzaName;
    private String l1bFlagsName;
    private int invalidBit;
    private String[] targetTpgs;

    Sensor(String name, String szaName, String vzaName, String l1bFlagsName, int invalidBit, String[] targetTpgs) {
        this.name = name;
        this.szaName = szaName;
        this.vzaName = vzaName;
        this.l1bFlagsName = l1bFlagsName;
        this.invalidBit = invalidBit;
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

    public String getL1bFlagsName() {
        return l1bFlagsName;
    }

    public int getInvalidBit() {
        return invalidBit;
    }

    public String[] getTargetTpgs() {
        return targetTpgs;
    }
}

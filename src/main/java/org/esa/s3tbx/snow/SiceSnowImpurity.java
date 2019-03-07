package org.esa.s3tbx.snow;

/**
 * Holder for SICE snow impurity:
 * - relative impurity load (volumetric concentration of pollutants)
 * - absorption Angstroem exponent
 * - normalized absorption coefficient at 1000nm
 * - pollution type (flag: 1-soot, 2-dust,3-algae,0-uncertain)
 *
 * @author olafd
 */
public class SiceSnowImpurity {
    private double relativeImpurityLoad;
    private double absorptionAngstromExp;
    private double normalizedAbsCoeff;
    private int pollutionType;

    SiceSnowImpurity(double relativeImpurityLoad,
                     double absorptionAngstromExp,
                     double normalizedAbsCoeff,
                     int pollutionType) {
        this.relativeImpurityLoad = relativeImpurityLoad;
        this.absorptionAngstromExp = absorptionAngstromExp;
        this.normalizedAbsCoeff = normalizedAbsCoeff;
        this.pollutionType = pollutionType;
    }

    public double getRelativeImpurityLoad() {
        return relativeImpurityLoad;
    }

    public double getAbsorptionAngstromExp() {
        return absorptionAngstromExp;
    }

    public double getNormalizedAbsCoeff() {
        return normalizedAbsCoeff;
    }

    public int getPollutionType() {
        return pollutionType;
    }
}

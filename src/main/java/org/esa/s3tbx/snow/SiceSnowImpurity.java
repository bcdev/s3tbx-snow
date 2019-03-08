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
    private double concentrationOfPollutants;
    private double absorptionAngstromExp;
    private double normalizedAbsCoeff;
    private int pollutionType;

    SiceSnowImpurity(double concentrationOfPollutants,
                     double absorptionAngstromExp,
                     double normalizedAbsCoeff,
                     int pollutionType) {
        this.concentrationOfPollutants = concentrationOfPollutants;
        this.absorptionAngstromExp = absorptionAngstromExp;
        this.normalizedAbsCoeff = normalizedAbsCoeff;
        this.pollutionType = pollutionType;
    }

    public double getConcentrationOfPollutants() {
        return concentrationOfPollutants;
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

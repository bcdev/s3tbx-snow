/*
 *
 *  * Copyright (C) 2012 Brockmann Consult GmbH (info@brockmann-consult.de)
 *  *
 *  * This program is free software; you can redistribute it and/or modify it
 *  * under the terms of the GNU General Public License as published by the Free
 *  * Software Foundation; either version 3 of the License, or (at your option)
 *  * any later version.
 *  * This program is distributed in the hope that it will be useful, but WITHOUT
 *  * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 *  * more details.
 *  *
 *  * You should have received a copy of the GNU General Public License along
 *  * with this program; if not, see http://www.gnu.org/licenses/
 *
 */

package org.esa.s3tbx.snow;

import com.bc.ceres.core.ProgressMonitor;
import org.apache.commons.lang.ArrayUtils;
import org.esa.s3tbx.olci.radiometry.rayleigh.RayleighCorrectionOp;
import org.esa.snap.core.datamodel.Band;
import org.esa.snap.core.datamodel.FlagCoding;
import org.esa.snap.core.datamodel.Product;
import org.esa.snap.core.datamodel.ProductData;
import org.esa.snap.core.gpf.Operator;
import org.esa.snap.core.gpf.OperatorException;
import org.esa.snap.core.gpf.OperatorSpi;
import org.esa.snap.core.gpf.Tile;
import org.esa.snap.core.gpf.annotations.OperatorMetadata;
import org.esa.snap.core.gpf.annotations.Parameter;
import org.esa.snap.core.gpf.annotations.SourceProduct;
import org.esa.snap.core.util.ProductUtils;
import org.esa.snap.core.util.math.MathUtils;

import java.awt.*;
import java.util.Map;
import java.util.Set;

/**
 * Computes snow properties from OLCI L1b data products.
 *
 * @author olafd
 */
@OperatorMetadata(alias = "OLCI.SnowProperties",
        description = "Computes snow properties from OLCI L1b data products.",
        authors = "Olaf Danne (Brockmann Consult), Alexander Kokhanovsky (Vitrociset)",
        copyright = "(c) 2017, 2018 by ESA, Brockmann Consult",
        category = "Optical/Thematic Land Processing",
        version = "2.0.10-SNAPSHOT")

public class OlciSnowPropertiesOp extends Operator {

    private static final String ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX = "albedo_spectral_spherical_";
    private static final String ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX = "albedo_spectral_planar_";
    private static final String PPA_SPECTRAL_OUTPUT_PREFIX = "ppa_spectral_";

    private static final String ALBEDO_BB_OUTPUT_PREFIX = "albedo_bb_";

    private static final String ALBEDO_BROADBAND_VIS_SPHERICAL_BAND_NAME = "albedo_bb_spherical_vis";
    private static final String ALBEDO_BROADBAND_NIR_SPHERICAL_BAND_NAME = "albedo_bb_spherical_nir";
    private static final String ALBEDO_BROADBAND_SW_SPHERICAL_BAND_NAME = "albedo_bb_spherical_sw";

    private static final String ALBEDO_BROADBAND_VIS_PLANAR_BAND_NAME = "albedo_bb_planar_vis";
    private static final String ALBEDO_BROADBAND_NIR_PLANAR_BAND_NAME = "albedo_bb_planar_nir";
    private static final String ALBEDO_BROADBAND_SW_PLANAR_BAND_NAME = "albedo_bb_planar_sw";

    private static final String[] ALBEDO_BROADBAND_SUFFIXES = {"vis", "nir", "sw"};

    private static final String GRAIN_DIAMETER_BAND_NAME = "grain_diameter";
    private static final String SNOW_SPECIFIC_AREA_BAND_NAME = "snow_specific_area";
    private static final String NDBI_BAND_NAME = "ndbi";
    private static final String POLLUTION_MASK_BAND_NAME = "pollution_mask";
    private static final String POLLUTION_F_BAND_NAME = "f";
    private static final String POLLUTION_L_BAND_NAME = "l";
    private static final String POLLUTION_M_BAND_NAME = "m";
    private static final String POLLUTION_R0_BAND_NAME = "r_0";
    private static final String POLLUTION_F_REL_ERR_BAND_NAME = "f_rel_err";
    private static final String POLLUTION_L_REL_ERR_BAND_NAME = "l_rel_err";
    private static final String POLLUTION_M_REL_ERR_BAND_NAME = "m_rel_err";
    private static final String POLLUTION_R0_REL_ERR_BAND_NAME = "r_0_rel_err";
    private static final String NDSI_MASK_BAND_NAME = "ndsi_mask";
    private static final String NDSI_BAND_NAME = "ndsi";


    @Parameter(description = "The OLCI wavelengths for spectral snow quantities " +
            "(spherical and planar albedos, PPA) which will be written to the target product. " +
            "Note that processing may become slow if many wavelenghts are selected.",
            label = "Select OLCI wavelengths for spectral snow quantities",
            valueSet = {
                    "Oa01 (400 nm)", "Oa02 (412.5 nm)", "Oa03 (442.5 nm)", "Oa04 (490 nm)", "Oa05 (510 nm)",
                    "Oa06 (560 nm)", "Oa07 (620 nm)", "Oa08 (665 nm)", "Oa09 (673.75 nm)", "Oa10 (681.25 nm)",
                    "Oa11 (708.75 nm)", "Oa12 (753.75 nm)", "Oa13 (761.25 nm)", "Oa14 (764.375 nm)", "Oa15 (767.5 nm)",
                    "Oa16 (778.75 nm)", "Oa17 (865 nm)", "Oa18 (885 nm)", "Oa19 (900 nm)", "Oa20 (940 nm)",
                    "Oa21 (1020 nm)"
            },
            defaultValue = "")
    private String[] spectralAlbedoTargetBands;

    @Parameter(description = "Name of binary mask band in cloud mask product (if present)",
            label = "Name of binary mask band in cloud mask product (if present)",
            defaultValue = "cloud_over_snow")
    private String cloudMaskBandName;

    @Parameter(defaultValue = "false",
            label = "Consider NDSI snow mask",
            description =
                    "If selected, NDSI will be computed from 865 and 1020nm for snow identification. " +
                            "Then, only snow pixels will be considered for albedo computation.")
    private boolean considerNdsiSnowMask;

    @Parameter(defaultValue = "0.03",
            description = "NDSI threshold for snow identification",
            label = "NDSI threshold for snow identification")
    private double ndsiThresh;

    @Parameter(defaultValue = "true",
            description =
                    "If selected, polluted snow will be retrieved and specific algorithms for albedo retrieval will be used.")
    private boolean considerSnowPollution;

    @Parameter(defaultValue = "false",
            label = "Write additional parameters (f, l, m, R_0) in case of polluted snow (expert option)",
            description =
                    "If selected, additional parameters (f, l, m, R_0) will be written to the target product in case " +
                            "of polluted snow. Option 'Consider snow pollution' must be selected as well. Useful for experts only.")
    private boolean writeAdditionalSnowPollutionParms;

    @Parameter(defaultValue = "false",
            label = "Write uncertainties of additional parameters in case of polluted snow (expert option)",
            description =
                    "If selected, uncertainties of additional parameters will be written to the target product in case " +
                            "of polluted snow. Options 'Consider snow pollution' and 'Write additional parameters " +
                            "(f, l, m, R_0)...' must be selected as well. Useful for experts only.")
    private boolean writeUncertaintiesOfAdditionalSnowPollutionParms;

    @Parameter(defaultValue = "0.01",
            description = "Assumed uncertainty in Rayleigh corrected reflectances",
            label = "Assumed uncertainty of Rayleigh corrected reflectances")
    private double deltaBrr;

    @Parameter(defaultValue = "0.1",
            description = "Snow is regarded as polluted if snow reflectance at 400nm is smaller that R_0 - thresh. " +
                    "See algorithm descriptions for more details.",
            label = "Snow pollution threshold")
    private double pollutionDelta;

    @Parameter(defaultValue = "false",
            label = "Compute PPA",
            description =
                    "If selected, PPA (Probability of Photon Absorption) is computed at selected OLCI wavelengths" +
                            " and written to target product")
    private boolean computePPA;

    @Parameter(defaultValue = "false",
            description =
                    "If selected, Rayleigh corrected reflectances at selected OLCI wavelengths are written to target product")
    private boolean copyReflectanceBands;

    @Parameter(defaultValue = "1020.0",
            valueSet = {"1020.0", "865.0"},
            description = "OLCI reference wavelength used in computations of snow quantities",
            label = "OLCI reference wavelength (nm)")
    private double refWvl;

    @Parameter(defaultValue = "0.9798",
            description = "OLCI SVC gain for band 1 (default value as provided by Sentinel-3A Product Notice – " +
                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
            label = "OLCI SVC gain for band 1 (400nm)")
    private double olciGainBand1;

    @Parameter(defaultValue = "0.9892",
            description = "OLCI SVC gain for band 6 (default value as provided by Sentinel-3A Product Notice – " +
                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
            label = "OLCI SVC gain for band 6 (560nm)")
    private double olciGainBand5;

    @Parameter(defaultValue = "1.0",
            description = "OLCI SVC gain for band 17 (default value as provided by Sentinel-3A Product Notice – " +
                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
            label = "OLCI SVC gain for band 17 (865nm)")
    private double olciGainBand17;

    @Parameter(defaultValue = "0.914",
            description = "OLCI SVC gain for band 21 (default value as provided by Sentinel-3A Product Notice – " +
                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
            label = "OLCI SVC gain for band 21 (1020nm)")
    private double olciGainBand21;

    @SourceProduct(description = "OLCI L1b or Rayleigh corrected product",
            label = "OLCI L1b or Rayleigh corrected product")
    private Product sourceProduct;

    @SourceProduct(description = "Cloud over snow binary mask product",
            label = "Cloud mask product",
            optional = true)
    private Product cloudMaskProduct;

    private Sensor sensor = Sensor.OLCI;
    private double[] olciGains;

    private String[] requiredBrrBandNamesAlbedo;

    private String[] requiredRadianceBandNamesPPA;     // selected spectral bands
    private String[] requiredBrrBandNamesPPA;

    private Product targetProduct;


    private Product reflProduct;

    private int width;
    private int height;

    //    private SolarSpectrumTable solarSpectrumTable;
    private SolarSpectrumExtendedTable solarSpectrumExtendedTable;
    private RefractiveIndexTable refractiveIndexInterpolatedTable;


    @Override
    public void initialize() throws OperatorException {
        String[] requiredRadianceBandNamesAlbedo;

        olciGains = new double[4];
        olciGains[0] = olciGainBand1;
        olciGains[1] = olciGainBand5;
        olciGains[2] = olciGainBand17;
        olciGains[3] = olciGainBand21;
        requiredRadianceBandNamesAlbedo = new String[]{"Oa01_radiance", "Oa06_radiance", "Oa17_radiance", "Oa21_radiance"};
        requiredBrrBandNamesAlbedo = new String[]{"rBRR_01", "rBRR_06", "rBRR_17", "rBRR_21"};

        // get required radiance / BRR bands for PPA computation
        if (spectralAlbedoTargetBands != null) {
            requiredRadianceBandNamesPPA = new String[spectralAlbedoTargetBands.length];
            requiredBrrBandNamesPPA = new String[spectralAlbedoTargetBands.length];
            for (int i = 0; i < spectralAlbedoTargetBands.length; i++) {
                // 'Oa01 (400 nm)' --> 'Oa01_radiance'
                // 'Oa01 (400 nm)' --> 'rBRR_01'
                requiredRadianceBandNamesPPA[i] = "Oa" + spectralAlbedoTargetBands[i].substring(2, 4) + "_radiance";
                requiredBrrBandNamesPPA[i] = "rBRR_" + spectralAlbedoTargetBands[i].substring(2, 4);
            }
        }

        validateSourceProduct(sourceProduct);

        width = sourceProduct.getSceneRasterWidth();
        height = sourceProduct.getSceneRasterHeight();

        if (isValidL1bSourceProduct(sourceProduct)) {
            // apply Rayleigh correction
            RayleighCorrectionOp rayleighCorrectionOp = new RayleighCorrectionOp();
            rayleighCorrectionOp.setSourceProduct(sourceProduct);
            rayleighCorrectionOp.setParameterDefaultValues();
            rayleighCorrectionOp.setParameter("computeTaur", false);
            final String[] sourceBandNames = SnowUtils.setupRcSourceBands(requiredRadianceBandNamesAlbedo,
                                                                          requiredRadianceBandNamesPPA);
            rayleighCorrectionOp.setParameter("sourceBandNames", sourceBandNames);
            reflProduct = rayleighCorrectionOp.getTargetProduct();
        } else {
            reflProduct = sourceProduct;
        }

        // read auxiliary data:
        RefractiveIndexTable refractiveIndexTable = new RefractiveIndexTable();
        solarSpectrumExtendedTable = new SolarSpectrumExtendedTable();

        // interpolate input refractive indices (at 83 wavelengths) to full grid 0.3-1.02um from solar spectrum auxdata
        refractiveIndexInterpolatedTable = SnowUtils.getRefractiveIndexInterpolated(refractiveIndexTable,
                                                                                    solarSpectrumExtendedTable);

        if (cloudMaskProduct != null) {
            validateCloudMaskProduct();
        }

        createTargetProduct();
    }

    @Override
    public void computeTileStack(Map<Band, Tile> targetTiles, Rectangle targetRectangle, ProgressMonitor pm) throws OperatorException {
        if (targetRectangle.x == 0 && targetRectangle.y == 0) {
            System.out.println("calling computeTileStack: targetRectangle = " + targetRectangle);

            Set<Map.Entry<Band, Tile>> entries = targetTiles.entrySet();
            entries.forEach(targetTileStream -> {
                Band targetBand = targetTileStream.getKey();
                System.out.println("targetBand.getName() = " + targetBand.getName());
            });
        }
        try {
            Tile[] rhoToaTilesAlbedo = new Tile[requiredBrrBandNamesAlbedo.length];
            for (int i = 0; i < requiredBrrBandNamesAlbedo.length; i++) {
                final Band rhoToaBandAlbedo = reflProduct.getBand(requiredBrrBandNamesAlbedo[i]);
                rhoToaTilesAlbedo[i] = getSourceTile(rhoToaBandAlbedo, targetRectangle);
            }
            Tile[] rhoToaTilesPPA = null;
            if (computePPA && spectralAlbedoTargetBands != null) {
                rhoToaTilesPPA = new Tile[requiredBrrBandNamesPPA.length];
                for (int i = 0; i < requiredBrrBandNamesPPA.length; i++) {
                    final Band rhoToaBandPPA = reflProduct.getBand(requiredBrrBandNamesPPA[i]);
                    rhoToaTilesPPA[i] = getSourceTile(rhoToaBandPPA, targetRectangle);
                }
            }

            Tile szaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getSzaName()), targetRectangle);
            Tile saaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getSaaName()), targetRectangle);
            Tile vzaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getVzaName()), targetRectangle);
            Tile vaaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getVaaName()), targetRectangle);
            Tile l1FlagsTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getL1bFlagsName()), targetRectangle);

            Tile cloudMaskTile = null;
            if (cloudMaskProduct != null) {
                cloudMaskTile = getSourceTile(cloudMaskProduct.getRasterDataNode(cloudMaskBandName), targetRectangle);
            }

            for (int y = targetRectangle.y; y < targetRectangle.y + targetRectangle.height; y++) {
                checkForCancellation();
                for (int x = targetRectangle.x; x < targetRectangle.x + targetRectangle.width; x++) {

                    double[] rhoToaAlbedo = new double[requiredBrrBandNamesAlbedo.length];   // now always 01, 05, 17, 21
                    for (int i = 0; i < requiredBrrBandNamesAlbedo.length; i++) {
                        rhoToaAlbedo[i] = olciGains[i] * rhoToaTilesAlbedo[i].getSampleDouble(x, y);
                        rhoToaAlbedo[i] = Math.max(0.0, rhoToaAlbedo[i]);
                    }

                    double brr400 = rhoToaAlbedo[0];
                    final double brr865 = rhoToaAlbedo[2];
                    final double brr1020 = rhoToaAlbedo[3];

                    double ndbi = Double.NaN;
                    boolean isBareIce = false;
                    if (brr400 > 0.0 && brr1020 > 0.0) {
                        final double iceIndicator = brr400 / brr1020;
                        ndbi = (iceIndicator - 1) / (iceIndicator + 1);
                        ndbi = Math.max(-1.0, Math.min(ndbi, 1.0));
                        isBareIce = ndbi > OlciSnowPropertiesConstants.BARE_ICE_THRESH;
                    }

                    final boolean l1Valid = !l1FlagsTile.getSampleBit(x, y, sensor.getInvalidBit());
                    final boolean isLandOrBareIce = l1FlagsTile.getSampleBit(x, y, sensor.getLandBit()) || isBareIce;
                    final boolean isNotCloud = cloudMaskTile == null ||
                            (cloudMaskTile != null && cloudMaskTile.getSampleDouble(x, y) != 1.0);

                    final boolean pixelIsValid = l1Valid && isLandOrBareIce && isNotCloud;

                    if (pixelIsValid) {

                        final Band ndbiBand = targetProduct.getBand(NDBI_BAND_NAME);
                        if (!Double.isNaN(ndbi)) {
                            targetTiles.get(ndbiBand).setSample(x, y, SnowUtils.cutTo4DecimalPlaces(ndbi));
                        } else {
                            targetTiles.get(ndbiBand).setSample(x, y, ndbi);
                        }

                        double ndsi = (brr865 - brr1020) / (brr865 + brr1020);
                        boolean validNdsi = true;
                        if (considerNdsiSnowMask) {
                            if (ndsi <= ndsiThresh || brr400 <= 0.5) {
                                validNdsi = false;
                            }
                        }
                        final Band ndsiBand = targetProduct.getBand(NDSI_BAND_NAME);
                        if (!Double.isNaN(ndsi)) {
                            targetTiles.get(ndsiBand).setSample(x, y, SnowUtils.cutTo4DecimalPlaces(ndsi));
                        } else {
                            targetTiles.get(ndsiBand).setSample(x, y, ndsi);
                        }

                        if (validNdsi) {
                            double[] rhoToaPPA = null;
                            if (computePPA && spectralAlbedoTargetBands != null) {
                                rhoToaPPA = new double[requiredBrrBandNamesPPA.length];
                                for (int i = 0; i < requiredBrrBandNamesPPA.length; i++) {
                                    if (rhoToaTilesPPA != null) {
                                        rhoToaPPA[i] = rhoToaTilesPPA[i].getSampleDouble(x, y);
                                        rhoToaPPA[i] = Math.max(0.0, rhoToaPPA[i]);
                                    }
                                }
                            }

                            final double sza = szaTile.getSampleDouble(x, y);
                            final double vza = vzaTile.getSampleDouble(x, y);
                            final double mu_0 = Math.cos(sza * MathUtils.DTOR);

                            double[][] spectralAlbedos;

                            final double saa = saaTile.getSampleDouble(x, y);
                            final double vaa = vaaTile.getSampleDouble(x, y);
                            final double raa = SnowUtils.getRelAzi(saa, vaa);
                            final double r0Thresh = OlciSnowPropertiesAlgorithm.computeR0ReflectancePollutionThresh(sza, vza, raa);
                            final double pollutedSnowSeparationValue = r0Thresh - pollutionDelta;
                            final boolean isPollutedSnow = considerSnowPollution && brr400 < pollutedSnowSeparationValue;

                            SpectralAlbedoResult spectralAlbedoResult;

                            if (considerSnowPollution) {
                                if (isPollutedSnow) {
                                    spectralAlbedoResult =
                                            OlciSnowPropertiesAlgorithm.computeSpectralAlbedosPolluted(rhoToaAlbedo,
                                                                                                       sza, vza,
                                                                                                       deltaBrr);
                                } else {
                                    spectralAlbedoResult =
                                            OlciSnowPropertiesAlgorithm.computeSpectralAlbedos(rhoToaAlbedo, deltaBrr, sza, vza);
                                }
                            } else {
                                spectralAlbedoResult =
                                        OlciSnowPropertiesAlgorithm.computeSpectralAlbedos(rhoToaAlbedo, deltaBrr, sza, vza);
                            }

                            spectralAlbedos = spectralAlbedoResult.getSpectralAlbedos();

                            final double[] spectralSphericalAlbedos = spectralAlbedos[0];
                            final double[] spectralPlanarAlbedos = spectralAlbedos[1];

                            final double r0 = spectralAlbedoResult.getR0();
                            final double f = spectralAlbedoResult.getF();
                            final double l = spectralAlbedoResult.getL();
                            final double m = spectralAlbedoResult.getM();

                            setTargetTilesSpectralAlbedos(spectralSphericalAlbedos,
                                                          ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX, targetTiles, x, y);
                            setTargetTilesSpectralAlbedos(spectralPlanarAlbedos,
                                                          ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX, targetTiles, x, y);

                            final double refAlbedo = refWvl == 1020.0 ?
                                    spectralSphericalAlbedos[spectralSphericalAlbedos.length - 1] :
                                    spectralSphericalAlbedos[spectralSphericalAlbedos.length - 5];

                            double grainDiam = Double.NaN;
                            if (refAlbedo > 0.0) {
//                                    grainDiam = l / 13.08;
                                grainDiam = l / 11.38;      //  AK, October 2018. l and grainDiam are in mm !!
                            }
                            final double[] broadbandPlanarAlbedo =
                                    OlciSnowPropertiesAlgorithm.computeBroadbandAlbedo(mu_0,
                                                                                       rhoToaAlbedo,
                                                                                       isPollutedSnow,
                                                                                       refractiveIndexInterpolatedTable,
                                                                                       solarSpectrumExtendedTable,
                                                                                       sza, vza);
                            final double[] broadbandSphericalAlbedo =
                                    OlciSnowPropertiesAlgorithm.computeBroadbandAlbedo(1.0,
                                                                                       rhoToaAlbedo,
                                                                                       isPollutedSnow,
                                                                                       refractiveIndexInterpolatedTable,
                                                                                       solarSpectrumExtendedTable,
                                                                                       sza, vza);

                            setTargetTilesBroadbandAlbedos(broadbandPlanarAlbedo, targetTiles, "planar", x, y);
                            setTargetTilesBroadbandAlbedos(broadbandSphericalAlbedo, targetTiles, "spherical", x, y);

                            if (computePPA && spectralAlbedoTargetBands != null) {
                                final double[] spectralPPA = OlciSnowPropertiesAlgorithm.computeSpectralPPA(rhoToaPPA, sza, vza);
                                setTargetTilesSpectralPPA(spectralPPA, PPA_SPECTRAL_OUTPUT_PREFIX, targetTiles, x, y);
                            }

                            final Band grainDiameterBand = targetProduct.getBand(GRAIN_DIAMETER_BAND_NAME);
                            final Band snowSpecificAreaBand = targetProduct.getBand(SNOW_SPECIFIC_AREA_BAND_NAME);
                            if (!Double.isNaN(grainDiam)) {
                                final double grainDiamMetres = grainDiam / 1000.0;  // in m
                                targetTiles.get(grainDiameterBand).setSample(x, y, SnowUtils.cutTo4DecimalPlaces(grainDiam));

                                final double snowSpecificArea = 6.0 / (OlciSnowPropertiesConstants.RHO_ICE * grainDiamMetres);
                                targetTiles.get(snowSpecificAreaBand).setSample(x, y, SnowUtils.cutTo7DecimalPlaces(snowSpecificArea));
                            } else {
                                targetTiles.get(grainDiameterBand).setSample(x, y, Double.NaN);
                                targetTiles.get(snowSpecificAreaBand).setSample(x, y, Double.NaN);
                            }

                            if (considerSnowPollution) {
                                final Band pollutionMaskBand = targetProduct.getBand(POLLUTION_MASK_BAND_NAME);
                                if (brr400 > 0.0 && !Double.isNaN(r0)) {
                                    int pollutionMask = isPollutedSnow ? 1 : 0;
                                    targetTiles.get(pollutionMaskBand).setSample(x, y, pollutionMask);
                                } else {
                                    targetTiles.get(pollutionMaskBand).setSample(x, y, 0);
                                }
                                if (writeAdditionalSnowPollutionParms) {
                                    final Band pollutionFBand = targetProduct.getBand(POLLUTION_F_BAND_NAME);
                                    targetTiles.get(pollutionFBand).setSample(x, y, isPollutedSnow ? f : Float.NaN);
                                    final Band pollutionLBand = targetProduct.getBand(POLLUTION_L_BAND_NAME);
                                    targetTiles.get(pollutionLBand).setSample(x, y, isPollutedSnow ? l : Float.NaN);
                                    final Band pollutionMBand = targetProduct.getBand(POLLUTION_M_BAND_NAME);
                                    targetTiles.get(pollutionMBand).setSample(x, y, isPollutedSnow ? m : Float.NaN);
                                    final Band pollutionR0Band = targetProduct.getBand(POLLUTION_R0_BAND_NAME);
                                    targetTiles.get(pollutionR0Band).setSample(x, y, isPollutedSnow ? r0 : Float.NaN);

                                    if (writeUncertaintiesOfAdditionalSnowPollutionParms) {
                                        final double r0RelErr = spectralAlbedoResult.getR0RelErr();
                                        final double fRelErr = spectralAlbedoResult.getfRelErr();
                                        final double lRelErr = spectralAlbedoResult.getlRelErr();
                                        final double mRelErr = spectralAlbedoResult.getmRelErr();
                                        final Band pollutionFRelErrBand = targetProduct.getBand(POLLUTION_F_REL_ERR_BAND_NAME);
                                        targetTiles.get(pollutionFRelErrBand).setSample(x, y, isPollutedSnow ? fRelErr : Float.NaN);
                                        final Band pollutionLRelErrBand = targetProduct.getBand(POLLUTION_L_REL_ERR_BAND_NAME);
                                        targetTiles.get(pollutionLRelErrBand).setSample(x, y, isPollutedSnow ? lRelErr : Float.NaN);
                                        final Band pollutionMRelErrBand = targetProduct.getBand(POLLUTION_M_REL_ERR_BAND_NAME);
                                        targetTiles.get(pollutionMRelErrBand).setSample(x, y, isPollutedSnow ? mRelErr : Float.NaN);
                                        final Band pollutionR0RelErrBand = targetProduct.getBand(POLLUTION_R0_REL_ERR_BAND_NAME);
                                        targetTiles.get(pollutionR0RelErrBand).setSample(x, y, isPollutedSnow ? r0RelErr : Float.NaN);
                                    }
                                }
                            }
                            if (considerNdsiSnowMask) {
                                final Band ndsiMaskBand = targetProduct.getBand(NDSI_MASK_BAND_NAME);
                                int ndsiMask = validNdsi ? 1 : 0;
                                targetTiles.get(ndsiMaskBand).setSample(x, y, ndsiMask);
                            }
                        } else {
                            setTargetTilesInvalid(targetTiles, x, y);
                        }

                    } else {
                        setTargetTilesInvalid(targetTiles, x, y);
                    }
                }
            }
        } catch (Exception e) {
            throw new OperatorException(e);
        }
    }

    private void createTargetProduct() {
        targetProduct = new Product(sourceProduct.getName(), sourceProduct.getProductType(), width, height);

        targetProduct.addBand(ALBEDO_BROADBAND_VIS_SPHERICAL_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_NIR_SPHERICAL_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_SW_SPHERICAL_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_VIS_PLANAR_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_NIR_PLANAR_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_SW_PLANAR_BAND_NAME, ProductData.TYPE_FLOAT32);

        targetProduct.addBand(GRAIN_DIAMETER_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(SNOW_SPECIFIC_AREA_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(NDBI_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(NDSI_BAND_NAME, ProductData.TYPE_FLOAT32);

        if (considerSnowPollution) {
            targetProduct.addBand(POLLUTION_MASK_BAND_NAME, ProductData.TYPE_INT16);
            if (writeAdditionalSnowPollutionParms) {
                targetProduct.addBand(POLLUTION_F_BAND_NAME, ProductData.TYPE_FLOAT32);
                targetProduct.addBand(POLLUTION_L_BAND_NAME, ProductData.TYPE_FLOAT32);
                targetProduct.addBand(POLLUTION_M_BAND_NAME, ProductData.TYPE_FLOAT32);
                targetProduct.addBand(POLLUTION_R0_BAND_NAME, ProductData.TYPE_FLOAT32);

                if (writeUncertaintiesOfAdditionalSnowPollutionParms) {
                    targetProduct.addBand(POLLUTION_F_REL_ERR_BAND_NAME, ProductData.TYPE_FLOAT32);
                    targetProduct.addBand(POLLUTION_L_REL_ERR_BAND_NAME, ProductData.TYPE_FLOAT32);
                    targetProduct.addBand(POLLUTION_M_REL_ERR_BAND_NAME, ProductData.TYPE_FLOAT32);
                    targetProduct.addBand(POLLUTION_R0_REL_ERR_BAND_NAME, ProductData.TYPE_FLOAT32);
                }
            }
        }

        if (considerNdsiSnowMask) {
            targetProduct.addBand(NDSI_MASK_BAND_NAME, ProductData.TYPE_INT16);
        }

        if (spectralAlbedoTargetBands != null && spectralAlbedoTargetBands.length > 0) {
            for (final String targetBand : spectralAlbedoTargetBands) {
                final int spectralBandIndex = Integer.parseInt(targetBand.substring(2, 4));
                final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[spectralBandIndex - 1] * 1000.0;
                final Band sphericalBand =
                        targetProduct.addBand(ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX + (int) wvl, ProductData.TYPE_FLOAT32);
                sphericalBand.setSpectralWavelength((float) wvl);
                sphericalBand.setSpectralBandIndex(spectralBandIndex);
                final Band planarBand =
                        targetProduct.addBand(ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX + (int) wvl, ProductData.TYPE_FLOAT32);
                planarBand.setSpectralWavelength((float) wvl);
                planarBand.setSpectralBandIndex(spectralBandIndex);
                if (computePPA) {
                    final Band ppaBand =
                            targetProduct.addBand(PPA_SPECTRAL_OUTPUT_PREFIX + (int) wvl, ProductData.TYPE_FLOAT32);
                    ppaBand.setSpectralWavelength((float) wvl);
                    ppaBand.setSpectralBandIndex(spectralBandIndex);
                }
            }

            if (copyReflectanceBands) {
                final String[] allBrrBands = (String[]) ArrayUtils.addAll(requiredBrrBandNamesAlbedo,
                                                                          requiredBrrBandNamesPPA);
                for (Band band : reflProduct.getBands()) {
                    for (String brrBandName : allBrrBands) {
                        if (band.getName().equals(brrBandName) && !targetProduct.containsBand(brrBandName)) {
                            ProductUtils.copyBand(band.getName(), reflProduct, targetProduct, true);
                            ProductUtils.copyRasterDataNodeProperties(band, targetProduct.getBand(band.getName()));
                        }
                    }
                }
            }
        }

        for (Band band : targetProduct.getBands()) {
            band.setUnit("dl");
            band.setNoDataValue(Float.NaN);
            band.setNoDataValueUsed(true);
        }
        targetProduct.getBand(GRAIN_DIAMETER_BAND_NAME).setUnit("mm");
        targetProduct.getBand(SNOW_SPECIFIC_AREA_BAND_NAME).setUnit("m^2/kg");
        final Band pollutionLBand = targetProduct.getBand(POLLUTION_L_BAND_NAME);
        if (pollutionLBand != null) {
            pollutionLBand.setUnit("mm");
        }
        Band pollutionFBand = targetProduct.getBand(POLLUTION_F_BAND_NAME);
        if (pollutionFBand != null) {
            pollutionFBand.setUnit("1/mm");
        }

        if (considerSnowPollution) {
            targetProduct.getBand(POLLUTION_MASK_BAND_NAME).setNoDataValueUsed(false);
        }

        if (cloudMaskProduct != null) {
            ProductUtils.copyBand(cloudMaskBandName, cloudMaskProduct, targetProduct, true);
            ProductUtils.copyFlagBands(cloudMaskProduct, targetProduct, true);
            FlagCoding idepixFlagCoding =
                    cloudMaskProduct.getFlagCodingGroup().get(OlciSnowPropertiesConstants.IDEPIX_CLASSIF_BAND_NAME);
            ProductUtils.copyFlagCoding(idepixFlagCoding, targetProduct);
            ProductUtils.copyMasks(cloudMaskProduct, targetProduct);
        }

        for (String tpg : sensor.getTargetTpgs()) {
            if (!targetProduct.containsTiePointGrid(tpg)) {
                ProductUtils.copyTiePointGrid(tpg, sourceProduct, targetProduct);
            }
        }

        ProductUtils.copyMetadata(sourceProduct, targetProduct);
        ProductUtils.copyMasks(sourceProduct, targetProduct);
        ProductUtils.copyFlagBands(sourceProduct, targetProduct, true);
        ProductUtils.copyGeoCoding(sourceProduct, targetProduct);
        targetProduct.setAutoGrouping("rBRR:albedo_spectral_spherical:albedo_spectral_planar:ppa_spectral");

        setTargetProduct(targetProduct);
    }

    private void setTargetTilesSpectralAlbedos(double[] spectralAlbedos,
                                               String prefix,
                                               Map<Band, Tile> targetTiles,
                                               int x, int y) {
        if (spectralAlbedoTargetBands != null && spectralAlbedoTargetBands.length > 0) {
            for (final String targetBand : spectralAlbedoTargetBands) {
                final int spectralBandIndex = Integer.parseInt(targetBand.substring(2, 4));
                final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[spectralBandIndex - 1] * 1000.0;
                final Band spectralAlbedoBand = targetProduct.getBand(prefix + (int) wvl);
                final double result = spectralAlbedos[spectralBandIndex - 1];
                targetTiles.get(spectralAlbedoBand).setSample(x, y, SnowUtils.cutTo4DecimalPlaces(result));
            }
        }
    }

    private void setTargetTilesSpectralPPA(double[] spectralPPA,
                                           String prefix,
                                           Map<Band, Tile> targetTiles,
                                           int x, int y) {
        int spectralPPABandIndex = 0;
        if (spectralAlbedoTargetBands != null && spectralAlbedoTargetBands.length > 0) {
            for (final String targetBand : spectralAlbedoTargetBands) {
                final int spectralBandIndex = Integer.parseInt(targetBand.substring(2, 4));
                final double wvl = OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI[spectralBandIndex - 1] * 1000.0;
                final Band spectralAlbedoBand = targetProduct.getBand(prefix + (int) wvl);
                targetTiles.get(spectralAlbedoBand).setSample(x, y, spectralPPA[spectralPPABandIndex++]);
            }
        }
    }

    private void setTargetTilesBroadbandAlbedos(double[] bbAlbedos, Map<Band, Tile> targetTiles,
                                                String bbMode, int x, int y) {
        int index = 0;
        for (final String bbSuffix : ALBEDO_BROADBAND_SUFFIXES) {
            final Band targetBand = targetProduct.getBand(ALBEDO_BB_OUTPUT_PREFIX + bbMode + "_" + bbSuffix);
            final double result = bbAlbedos[index++];
            targetTiles.get(targetBand).setSample(x, y, SnowUtils.cutTo4DecimalPlaces(result));
        }
    }

    private void setTargetTilesInvalid(Map<Band, Tile> targetTiles, int x, int y) {
        for (Tile tile : targetTiles.values()) {
            tile.setSample(x, y, Float.NaN);
        }
    }

    private void validateSourceProduct(Product sourceProduct) {
        boolean isOlci = isValidL1bSourceProduct(sourceProduct);
        if (!isOlci) {
            isOlci = isValidRayleighCorrectedSourceProduct(sourceProduct);
            if (!isOlci) {
                throw new OperatorException("Source product not applicable to this operator.\n" +
                                                    "Only OLCI L1b or Rayleigh corrected products are currently supported");
            }
        }
    }

    private void validateCloudMaskProduct() {
        Band cloudMaskBand = cloudMaskProduct.getBand(cloudMaskBandName);
        if (cloudMaskBand == null) {
            throw new OperatorException("Specified cloud mask product does not contain a band named '" +
                                                cloudMaskBandName + "'. Please check.");
        }
        if (cloudMaskProduct.getSceneRasterWidth() != width || cloudMaskProduct.getSceneRasterHeight() != height) {
            throw new OperatorException("Dimensions of cloud mask product differ from source product. Please check.");
        }
    }

    private boolean isValidL1bSourceProduct(Product sourceProduct) {
        for (int i = 0; i < OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length; i++) {
            final String bandName = "Oa" + String.format("%02d", i + 1) + "_radiance";
            if (!sourceProduct.containsBand(bandName)) {
                return false;
            }
        }
        return true;
    }

    private boolean isValidRayleighCorrectedSourceProduct(Product sourceProduct) {
        final String[] allRequiredBrrBands = (String[]) ArrayUtils.addAll(requiredBrrBandNamesAlbedo,
                                                                          requiredBrrBandNamesPPA);
        for (String bandName : allRequiredBrrBands) {
            if (!sourceProduct.containsBand(bandName)) {
                if (!sourceProduct.containsBand(bandName)) {
                    throw new OperatorException("Source product is not a valid L1b product and cannot be handled as " +
                                                        "Rayleigh corrected product either, as it does not contain " +
                                                        "mandatory band '" + bandName + "'. \n Mandatory bands are " +
                                                        "'rBRR_*' for indices 1, 5, 17, 21 " +
                                                        "(400nm, 510, 865 and 1020nm), and in addition for " +
                                                        "all manually selected wavelengths.");
                }
            }
        }
        return true;
    }


    /**
     * The Service Provider Interface (SPI) for the operator.
     * It provides operator meta-data and is a factory for new operator instances.
     */
    public static class Spi extends OperatorSpi {

        public Spi() {
            super(OlciSnowPropertiesOp.class);
        }
    }
}

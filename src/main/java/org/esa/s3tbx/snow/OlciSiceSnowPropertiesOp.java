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
import org.esa.snap.core.util.math.RsMathUtils;

import java.awt.*;
import java.util.Map;

/**
 * Computes snow properties from OLCI L1b data products.
 *
 * @author olafd
 */
@OperatorMetadata(alias = "OLCI.SnowProperties.SICE",
        description = "Computes snow properties from OLCI L1b data products, using new SICE algorithm.",
        authors = "Olaf Danne (Brockmann Consult), Alexander Kokhanovsky (Vitrociset)",
        copyright = "(c) 2019 by ESA, Brockmann Consult",
        category = "Optical/Thematic Land Processing",
        version = "3.0-SNAPSHOT")

public class OlciSiceSnowPropertiesOp extends Operator {

    private static final String SICE_POLLUTION_TYPE_FLAG_BAND_NAME = "sice_pollution_type_flags";
    private static final String SICE_GROUND_TYPE_FLAG_BAND_NAME = "sice_ground_type_flags";
    private static final String ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX = "albedo_spectral_spherical_";
    private static final String ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX = "albedo_spectral_planar_";

    private static final String ALBEDO_BB_OUTPUT_PREFIX = "albedo_bb_";

    private static final String ALBEDO_BROADBAND_VIS_SPHERICAL_BAND_NAME = "albedo_bb_spherical_vis";
    private static final String ALBEDO_BROADBAND_NIR_SPHERICAL_BAND_NAME = "albedo_bb_spherical_nir";
    private static final String ALBEDO_BROADBAND_SW_SPHERICAL_BAND_NAME = "albedo_bb_spherical_sw";

    private static final String ALBEDO_BROADBAND_VIS_PLANAR_BAND_NAME = "albedo_bb_planar_vis";
    private static final String ALBEDO_BROADBAND_NIR_PLANAR_BAND_NAME = "albedo_bb_planar_nir";
    private static final String ALBEDO_BROADBAND_SW_PLANAR_BAND_NAME = "albedo_bb_planar_sw";

    private static final String[] ALBEDO_BROADBAND_SUFFIXES = {"vis", "nir", "sw"};

    private static final String GRAIN_DIAMETER_BAND_NAME = "grain_diameter";
    private static final String SCATTERING_ANGLE_BAND_NAME = "scattering_angle";
    private static final String SNOW_SPECIFIC_AREA_BAND_NAME = "snow_specific_area";
    private static final String CONCENTRATION_OF_POLLUTANTS_BAND_NAME = "concentration_of_pollutants";
    private static final String EFFECTIVE_ABSORPTION_LENGTH_BAND_NAME = "efefctive_absorption_length";
    private static final String ABSORPTION_ANGSTROEM_EXPONENT_BAND_NAME = "absorption_angstroem_exponent";
    private static final String NORNALIZED_ABSORPTION_COEFFICIENT_BAND_NAME = "normalized_absorption_coefficient";
    private static final String NDBI_BAND_NAME = "ndbi";
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

    @Parameter(defaultValue = "false",
            label = "Consider NDSI snow mask",
            description = "If selected, NDSI will be computed from 865 and 1020nm for snow identification. " +
                    "Then, only 'NDSI snow' pixels will be considered for snow properties retrieval.")
    private boolean considerNdsiSnowMask;

    @Parameter(defaultValue = "0.03",
            description = "NDSI threshold for snow identification",
            label = "NDSI threshold for snow identification")
    private double ndsiThresh;

    // todo: discuss more options
//    @Parameter(defaultValue = "false",
//            label = "Write absorption Angstroem exponent",
//            description =
//                    "If selected, absorption Angstroem exponent will be written to the target product.")
//    private boolean writeAbsorptionAngstroemExponent;
    private boolean writeAbsorptionAngstroemExponent = false;

    //    @Parameter(defaultValue = "false",
//            label = "Write normalized absorption coefficient",
//            description =
//                    "If selected, normalized absorption coefficient will be written to the target product.")
//    private boolean writeNormalizedAbsorptionCoefficient;
    private boolean writeNormalizedAbsorptionCoefficient = false;

    //    @Parameter(defaultValue = "false",
//            label = "Write effective absorption length",
//            description =
//                    "If selected, effective absorption length will be written to the target product.")
//    private boolean writeEffectiveAbsorptionLength;
    private boolean writeEffectiveAbsorptionLength = false;

    //    @Parameter(defaultValue = "false",
//            label = "Copy Bottom-of-Atmosphere reflectance bands",
//            description =
//                    "If selected, Bottom-of-Atmosphere reflectance bands at selected OLCI wavelengths are written to target product")
//    private boolean copyBrrBands;
    private boolean copyBrrBands = false;

    @SourceProduct(description = "OLCI L1b radiance product",
            label = "OLCI L1b radiance product")
    private Product sourceProduct;

    @SourceProduct(description = "OLCI Rayleigh corrected product",
            label = "OLCI Rayleigh corrected product")
    private Product brrProduct;

    @SourceProduct(description = "Cloud over snow binary mask product",
            label = "Cloud mask product",
            optional = true)
    private Product cloudMaskProduct;

    private Sensor sensor = Sensor.OLCI;

    private String[] requiredBrrBandNames;

    private Product targetProduct;

    private int width;
    private int height;

    private RefractiveIndexTable refractiveIndexTable;


    @Override
    public void initialize() throws OperatorException {
        System.out.println("entering initialize...");
        requiredBrrBandNames = new String[]{"rBRR_01", "rBRR_06", "rBRR_10", "rBRR_11", "rBRR_17", "rBRR_21"};

        if (!isValidL1bSourceProduct(sourceProduct)) {
            throw new OperatorException("Source product is not a valid OLCI L1b radiance product");
        }

        validateRayleighCorrectedSourceProduct(brrProduct);

        width = sourceProduct.getSceneRasterWidth();
        height = sourceProduct.getSceneRasterHeight();

        // read auxiliary data:
        refractiveIndexTable = new RefractiveIndexTable();

        if (cloudMaskProduct != null) {
            validateCloudMaskProduct();
        }

        createTargetProduct();
        System.out.println("done initialize...");
    }

    @Override
    public void computeTileStack(Map<Band, Tile> targetTiles, Rectangle targetRectangle, ProgressMonitor pm) throws OperatorException {
        System.out.println("entering computeTileStack...: " + targetRectangle);
        try {
            Tile[] radianceTiles = new Tile[OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length];
            Tile[] fluxTiles = new Tile[OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length];
            for (int i = 0; i < OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length; i++) {
                String radianceBandName = "Oa" + String.format("%02d", i + 1) + "_radiance";
                String fluxBandName = "solar_flux_band_" + (i + 1);
                final Band radianceBand = sourceProduct.getBand(radianceBandName);
                final Band fluxBand = sourceProduct.getBand(fluxBandName);
                radianceTiles[i] = getSourceTile(radianceBand, targetRectangle);
                fluxTiles[i] = getSourceTile(fluxBand, targetRectangle);
            }

            Tile[] brrTiles = new Tile[requiredBrrBandNames.length];
            for (int i = 0; i < requiredBrrBandNames.length; i++) {
                final Band brrBand = brrProduct.getBand(requiredBrrBandNames[i]);
                brrTiles[i] = getSourceTile(brrBand, targetRectangle);
            }

            Tile szaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getSzaName()), targetRectangle);
            Tile saaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getSaaName()), targetRectangle);
            Tile vzaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getVzaName()), targetRectangle);
            Tile vaaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getVaaName()), targetRectangle);
            Tile l1FlagsTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getL1bFlagsName()), targetRectangle);

            Tile idepixClassifTile = null;
            if (cloudMaskProduct != null) {
                idepixClassifTile = getSourceTile(cloudMaskProduct.getRasterDataNode
                        (OlciSnowPropertiesConstants.IDEPIX_CLASSIF_BAND_NAME), targetRectangle);
            }

            final Band sicePollutionFlagBand = targetProduct.getBand(SICE_POLLUTION_TYPE_FLAG_BAND_NAME);
            final Tile sicePollutionFlagTile = targetTiles.get(sicePollutionFlagBand);
            final Band siceGroundFlagBand = targetProduct.getBand(SICE_GROUND_TYPE_FLAG_BAND_NAME);
            final Tile siceGroundFlagTile = targetTiles.get(siceGroundFlagBand);

            for (int y = targetRectangle.y; y < targetRectangle.y + targetRectangle.height; y++) {
                checkForCancellation();
                for (int x = targetRectangle.x; x < targetRectangle.x + targetRectangle.width; x++) {

                    double[] rhoToa = new double[OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length];
                    final double sza = szaTile.getSampleDouble(x, y);
                    for (int i = 0; i < OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length; i++) {
                        final float rad = radianceTiles[i].getSampleFloat(x, y);
                        final float flux = fluxTiles[i].getSampleFloat(x, y);
                        rhoToa[i] = RsMathUtils.radianceToReflectance(rad, (float) sza, flux);
                        rhoToa[i] = Math.max(0.0, rhoToa[i]);
                    }

                    final double rtoa400 = rhoToa[0];
                    final double rtoa865 = rhoToa[4];
                    final double rtoa1020 = rhoToa[5];

                    final double brr400 = brrTiles[0].getSampleDouble(x, y);
                    final double brr560 = brrTiles[1].getSampleDouble(x, y);
                    final double brr681 = brrTiles[2].getSampleDouble(x, y);
                    final double brr709 = brrTiles[3].getSampleDouble(x, y);
                    final double brr865 = brrTiles[4].getSampleDouble(x, y);
                    final double brr1020 = brrTiles[5].getSampleDouble(x, y);

                    final boolean l1Valid = !l1FlagsTile.getSampleBit(x, y, sensor.getInvalidBit());

                    boolean isCloud = false;
                    if (idepixClassifTile != null) {
                        isCloud = idepixClassifTile.getSampleBit(x, y, OlciSnowPropertiesConstants.IDEPIX_CLOUD) ||
                                idepixClassifTile.getSampleBit(x, y, OlciSnowPropertiesConstants.IDEPIX_CLOUD_BUFFER) ||
                                idepixClassifTile.getSampleBit(x, y, OlciSnowPropertiesConstants.IDEPIX_CLOUD_SHADOW);
                    }

                    // 20181207: do not exclude high SZA, but just raise a flag
                    final boolean pixelIsValid = l1Valid && !isCloud;

                    if (pixelIsValid) {
                        double ndsi = (rtoa865 - rtoa1020) / (rtoa865 + rtoa1020);
                        boolean validNdsi = true;
                        if (considerNdsiSnowMask) {
                            if (ndsi <= ndsiThresh || rtoa400 <= 0.5) {
                                validNdsi = false;
                            }
                        }
                        final Band ndsiBand = targetProduct.getBand(NDSI_BAND_NAME);
                        if (!Double.isNaN(ndsi)) {
                            targetTiles.get(ndsiBand).setSample(x, y, SnowUtils.cutTo4DecimalPlaces(ndsi));
                        } else {
                            targetTiles.get(ndsiBand).setSample(x, y, ndsi);
                        }

                        double ndbi = (rtoa400 - rtoa1020) / (rtoa400 + rtoa1020);
                        final Band ndbiBand = targetProduct.getBand(NDBI_BAND_NAME);
                        if (!Double.isNaN(ndbi)) {
                            targetTiles.get(ndbiBand).setSample(x, y, SnowUtils.cutTo4DecimalPlaces(ndbi));
                        } else {
                            targetTiles.get(ndbiBand).setSample(x, y, ndbi);
                        }

                        if (validNdsi) {
                            // all default computations:
                            final double vza = vzaTile.getSampleDouble(x, y);
                            final double saa = saaTile.getSampleDouble(x, y);
                            final double vaa = vaaTile.getSampleDouble(x, y);
                            final double raa = SnowUtils.getRelAziSice(saa, vaa);

                            double r0 = OlciSiceSnowPropertiesAlgorithm.computeR0(brr865, brr1020);
                            double xx = OlciSiceSnowPropertiesAlgorithm.computeXX(r0, sza, vza);
                            SiceSnowPropertiesResult siceSnowProperties =
                                    OlciSiceSnowPropertiesAlgorithm.computeGeneralSnowProperties
                                            (brr400, brr560, brr681, brr709, brr1020, r0, xx);
                            OlciSiceSnowPropertiesAlgorithm.computeSpectralAlbedos(siceSnowProperties,
                                                                                   rhoToa, brr400, sza, vza, raa);
                            OlciSiceSnowPropertiesAlgorithm.computeBroadbandAlbedos(siceSnowProperties,
                                                                                    brr400, sza,
                                                                                    refractiveIndexTable);

                            setTargetTilesSpectralAlbedos(siceSnowProperties.getSphericalSpectralAlbedos(),
                                                          ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX, targetTiles, x, y);
                            setTargetTilesSpectralAlbedos(siceSnowProperties.getPlanarSpectralAlbedos(),
                                                          ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX, targetTiles, x, y);

                            setTargetTilesBroadbandAlbedos(siceSnowProperties.getSphericalBroadbandAlbedos(),
                                                           targetTiles, "spherical", x, y);
                            setTargetTilesBroadbandAlbedos(siceSnowProperties.getPlanarBroadbandAlbedos(),
                                                           targetTiles, "planar", x, y);

                            final double snowGrainSize = siceSnowProperties.getSnowGrainSize();
                            final Band snowGrainSizeBand = targetProduct.getBand(GRAIN_DIAMETER_BAND_NAME);
                            targetTiles.get(snowGrainSizeBand).setSample(x, y,
                                                                         SnowUtils.cutTo4DecimalPlaces(snowGrainSize));

                            final double snowSpecificArea = siceSnowProperties.getSnowSpecificArea();
                            final Band snowSpecificAreaBand = targetProduct.getBand(SNOW_SPECIFIC_AREA_BAND_NAME);
                            targetTiles.get(snowSpecificAreaBand).setSample(x, y,
                                                                            SnowUtils.cutTo4DecimalPlaces(snowSpecificArea));

                            final double scatteringAngle = siceSnowProperties.getScatteringAngle();
                            final Band scatteringAngleBand = targetProduct.getBand(SCATTERING_ANGLE_BAND_NAME);
                            targetTiles.get(scatteringAngleBand).setSample(x, y,
                                                                           SnowUtils.cutTo4DecimalPlaces(scatteringAngle));

                            final double concentrationOfPollutants =
                                    siceSnowProperties.getSnowImpurity().getConcentrationOfPollutants();
                            final Band concOfPollutantsBand = targetProduct.getBand(CONCENTRATION_OF_POLLUTANTS_BAND_NAME);
                            targetTiles.get(concOfPollutantsBand).setSample(x, y,
                                                                            SnowUtils.cutTo4DecimalPlaces(concentrationOfPollutants));

                            // set flags:
                            final int pollutionTypeFlag =
                                    OlciSiceSnowPropertiesAlgorithm.computePollutionTypeFlag(siceSnowProperties, ndbi);
                            sicePollutionFlagTile.setSample(x, y, pollutionTypeFlag);
                            final int groundTypeFlag =
                                    OlciSiceSnowPropertiesAlgorithm.computeGroundTypeFlag(siceSnowProperties,
                                                                                          rtoa400, rtoa1020, ndsi, ndbi);
                            siceGroundFlagTile.setSample(x, y, groundTypeFlag);

                            // write optional output:
                            if (considerNdsiSnowMask) {
                                final Band ndsiMaskBand = targetProduct.getBand(NDSI_MASK_BAND_NAME);
                                int ndsiMask = validNdsi ? 1 : 0;
                                targetTiles.get(ndsiMaskBand).setSample(x, y, ndsiMask);
                            }

                            if (writeAbsorptionAngstroemExponent) {
                                final Band band = targetProduct.getBand(ABSORPTION_ANGSTROEM_EXPONENT_BAND_NAME);
                                final double absorptionAngstromExp =
                                        siceSnowProperties.getSnowImpurity().getAbsorptionAngstromExp();
                                targetTiles.get(band).setSample(x, y, absorptionAngstromExp);
                            }

                            if (writeNormalizedAbsorptionCoefficient) {
                                final Band band = targetProduct.getBand(NORNALIZED_ABSORPTION_COEFFICIENT_BAND_NAME);
                                final double normalizedAbsCoeff =
                                        siceSnowProperties.getSnowImpurity().getNormalizedAbsCoeff();
                                targetTiles.get(band).setSample(x, y, normalizedAbsCoeff);
                            }

                            if (writeEffectiveAbsorptionLength) {
                                final Band band = targetProduct.getBand(EFFECTIVE_ABSORPTION_LENGTH_BAND_NAME);
                                final double effAbsLength = siceSnowProperties.getEffAbsLength();
                                targetTiles.get(band).setSample(x, y, effAbsLength);
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
        targetProduct.addBand(SCATTERING_ANGLE_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(CONCENTRATION_OF_POLLUTANTS_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(NDBI_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(NDSI_BAND_NAME, ProductData.TYPE_FLOAT32);

        final Band sicePollutionFlagBand = targetProduct.addBand(SICE_POLLUTION_TYPE_FLAG_BAND_NAME, ProductData.TYPE_INT8);
        FlagCoding pollutionTypeFlagCoding =
                SnowUtils.createSicePollutionTypeFlagCoding(SICE_POLLUTION_TYPE_FLAG_BAND_NAME);
        sicePollutionFlagBand.setSampleCoding(pollutionTypeFlagCoding);
        targetProduct.getFlagCodingGroup().add(pollutionTypeFlagCoding);
        SnowUtils.setupSicePollutionTypeBitmask(targetProduct);

        final Band siceGroundFlagBand = targetProduct.addBand(SICE_GROUND_TYPE_FLAG_BAND_NAME, ProductData.TYPE_INT8);
        FlagCoding groundTypeFlagCoding =
                SnowUtils.createSiceGroundTypeFlagCoding(SICE_GROUND_TYPE_FLAG_BAND_NAME);
        siceGroundFlagBand.setSampleCoding(groundTypeFlagCoding);
        targetProduct.getFlagCodingGroup().add(groundTypeFlagCoding);
        SnowUtils.setupSiceGroundTypeBitmask(targetProduct);

        if (considerNdsiSnowMask) {
            targetProduct.addBand(NDSI_MASK_BAND_NAME, ProductData.TYPE_INT16);
        }

        if (writeAbsorptionAngstroemExponent) {
            targetProduct.addBand(ABSORPTION_ANGSTROEM_EXPONENT_BAND_NAME, ProductData.TYPE_FLOAT32);
        }

        if (writeNormalizedAbsorptionCoefficient) {
            targetProduct.addBand(NORNALIZED_ABSORPTION_COEFFICIENT_BAND_NAME, ProductData.TYPE_FLOAT32);
        }

        if (writeEffectiveAbsorptionLength) {
            targetProduct.addBand(EFFECTIVE_ABSORPTION_LENGTH_BAND_NAME, ProductData.TYPE_FLOAT32);
            targetProduct.getBand(EFFECTIVE_ABSORPTION_LENGTH_BAND_NAME).setUnit("mm");
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
            }

            if (copyBrrBands) {
                final String[] allBrrBands = (String[]) ArrayUtils.addAll(requiredBrrBandNames, null);
                for (Band band : brrProduct.getBands()) {
                    for (String brrBandName : allBrrBands) {
                        if (band.getName().equals(brrBandName) && !targetProduct.containsBand(brrBandName)) {
                            ProductUtils.copyBand(band.getName(), brrProduct, targetProduct, true);
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

        if (cloudMaskProduct != null) {
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
        targetProduct.setAutoGrouping("rBRR:albedo_spectral_spherical:albedo_spectral_planar");

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

    private void validateCloudMaskProduct() {
        Band pixelClassifBand = cloudMaskProduct.getBand(OlciSnowPropertiesConstants.IDEPIX_CLASSIF_BAND_NAME);
        if (pixelClassifBand == null) {
            throw new OperatorException("Specified cloud mask product does not contain the IdePix flag band '" +
                                                OlciSnowPropertiesConstants.IDEPIX_CLASSIF_BAND_NAME +
                                                "'. Please check.");
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

    private void validateRayleighCorrectedSourceProduct(Product sourceProduct) {
        for (String bandName : requiredBrrBandNames) {
            if (!sourceProduct.containsBand(bandName)) {
                if (!sourceProduct.containsBand(bandName)) {
                    throw new OperatorException("Source product is not a valid L1b product and cannot be handled as " +
                                                        "Rayleigh corrected product either, as it does not contain " +
                                                        "mandatory band '" + bandName + "'. \n Mandatory bands are " +
                                                        "'rBRR_*' for indices 1, 6, 10, 11, 17, 21 " +
                                                        "(400nm, 560, 681, 709, 865 and 1020nm), and in addition for " +
                                                        "all manually selected wavelengths.");
                }
            }
        }
    }

    /**
     * The Service Provider Interface (SPI) for the operator.
     * It provides operator meta-data and is a factory for new operator instances.
     */
    public static class Spi extends OperatorSpi {

        public Spi() {
            super(OlciSiceSnowPropertiesOp.class);
        }
    }
}

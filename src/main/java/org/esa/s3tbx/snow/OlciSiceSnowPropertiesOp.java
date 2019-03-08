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
import org.esa.s3tbx.processor.rad2refl.Rad2ReflOp;
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
        version = "0.8")

public class OlciSiceSnowPropertiesOp extends Operator {

    public static final String S3_SNOW_FLAG_BAND_NAME = "s3snow_flags";
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

    @Parameter(description = "Name of binary mask band in cloud mask product (if present)",
            label = "Name of binary mask band in cloud mask product (if present)",
            defaultValue = "cloud_over_snow")
    private String cloudMaskBandName;

    @Parameter(defaultValue = "false",
            label = "Consider NDSI snow mask",
            description = "If selected, NDSI will be computed from 865 and 1020nm for snow identification. " +
                    "Then, only 'NDSI snow' pixels will be considered for snow properties retrieval.")
    private boolean considerNdsiSnowMask;

    @Parameter(defaultValue = "0.03",
            description = "NDSI threshold for snow identification",
            label = "NDSI threshold for snow identification")
    private double ndsiThresh;

    @Parameter(defaultValue = "false",
            label = "Write absorption Angstroem exponent",
            description =
                    "If selected, absorption Angstroem exponent will be written to the target product.")
    private boolean writeAbsorptionAngstroemExponent;

    @Parameter(defaultValue = "false",
            label = "Write normalized absorption coefficient",
            description =
                    "If selected, normalized absorption coefficient will be written to the target product.")
    private boolean writeNormalizedAbsorptionCoefficient;

    @Parameter(defaultValue = "false",
            label = "Write effective absorption length",
            description =
                    "If selected, effective absorption length will be written to the target product.")
    private boolean writeEffectiveAbsorptionLength;

    @Parameter(defaultValue = "false",
            description =
                    "If selected, Rayleigh corrected reflectances at selected OLCI wavelengths are written to target product")
    private boolean copyBrrBands;

    // todo: check if gains are still requested
//    @Parameter(defaultValue = "0.9798",
//            description = "OLCI SVC gain for band 1 (default value as provided by Sentinel-3A Product Notice – " +
//                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
//            label = "OLCI SVC gain for band 1 (400nm)")
//    private double olciGainBand1;
    private double olciGainBand1 = 1.0;

    //    @Parameter(defaultValue = "0.9892",
//            description = "OLCI SVC gain for band 6 (default value as provided by Sentinel-3A Product Notice – " +
//                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
//            label = "OLCI SVC gain for band 6 (560nm)")
//    private double olciGainBand6;
    private double olciGainBand6 = 1.0;

    //    @Parameter(defaultValue = "0.9962",
//            description = "OLCI SVC gain for band 6 (default value as provided by Sentinel-3A Product Notice – " +
//                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
//            label = "OLCI SVC gain for band 10 (681nm)")
//    private double olciGainBand10;
    private double olciGainBand10 = 1.0;

    //    @Parameter(defaultValue = "0.996",
//            description = "OLCI SVC gain for band 6 (default value as provided by Sentinel-3A Product Notice – " +
//                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
//            label = "OLCI SVC gain for band 11 (709nm)")
//    private double olciGainBand11;
    private double olciGainBand11 = 1.0;

    //    @Parameter(defaultValue = "0.914",
//            description = "OLCI SVC gain for band 21 (default value as provided by Sentinel-3A Product Notice – " +
//                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
//            label = "OLCI SVC gain for band 21 (1020nm)")
//    private double olciGainBand21;
    private double olciGainBand21 = 1.0;

    @SourceProduct(description = "OLCI L1b or Rayleigh corrected product",
            label = "OLCI L1b or Rayleigh corrected product")
    private Product sourceProduct;

    @SourceProduct(description = "Cloud over snow binary mask product",
            label = "Cloud mask product",
            optional = true)
    private Product cloudMaskProduct;

    private Sensor sensor = Sensor.OLCI;
    private double[] olciGains;

    private String[] requiredBrrBandNames;

    private Product targetProduct;


    private Product rhoToaProduct;
    private Product brrProduct;

    private double[] wvlFullGrid;

    private int width;
    private int height;

    //    private SolarSpectrumTable solarSpectrumTable;
    private SolarSpectrumExtendedTable solarSpectrumExtendedTable;
    private RefractiveIndexTable refractiveIndexInterpolatedTable;
    private boolean validL1bSourceProduct;


    @Override
    public void initialize() throws OperatorException {
        String[] requiredRadianceBandNamesAlbedo;

        olciGains = new double[5];
        olciGains[0] = olciGainBand1;
        olciGains[1] = olciGainBand6;
        olciGains[2] = olciGainBand10;
        olciGains[3] = olciGainBand11;
        olciGains[4] = olciGainBand21;
        requiredRadianceBandNamesAlbedo =
                new String[]{"Oa01_radiance", "Oa06_radiance", "Oa10_radiance", "Oa11_radiance", "Oa17_radiance", "Oa21_radiance"};
        requiredBrrBandNames = new String[]{"rBRR_01", "rBRR_06", "rBRR_10", "rBRR_11", "rBRR_17", "rBRR_21"};

        validateSourceProduct(sourceProduct);

        width = sourceProduct.getSceneRasterWidth();
        height = sourceProduct.getSceneRasterHeight();

        validL1bSourceProduct = isValidL1bSourceProduct(sourceProduct);
        if (validL1bSourceProduct) {
            Rad2ReflOp rad2ReflOp = new Rad2ReflOp();
            rad2ReflOp.setSourceProduct(sourceProduct);
            rad2ReflOp.setParameterDefaultValues();
            rhoToaProduct = rad2ReflOp.getTargetProduct();      // band names: Oa%2d_reflectance

            // apply Rayleigh correction for 5 required bands
            RayleighCorrectionOp rayleighCorrectionOp = new RayleighCorrectionOp();
            rayleighCorrectionOp.setSourceProduct(sourceProduct);
            rayleighCorrectionOp.setParameterDefaultValues();
            rayleighCorrectionOp.setParameter("computeTaur", false);
            final String[] sourceBandNames = SnowUtils.setupRcSourceBands(requiredRadianceBandNamesAlbedo, null);
            rayleighCorrectionOp.setParameter("sourceBandNames", sourceBandNames);
            brrProduct = rayleighCorrectionOp.getTargetProduct();  // band names: rBRR_%2d
        } else {
            // in this case the source product MUST be a Rayleigh corrected product with all rBRR and rtoa bands
            // todo: check
            rhoToaProduct = sourceProduct;     // band names: rtoa_%2d   !!!
            brrProduct = sourceProduct;        // band names: rBRR_%2d
        }

        // set up fine resolution wavelength grid (5nm step)
        int numWvl = (int) ((OlciSnowPropertiesConstants.BB_WVL_3 - OlciSnowPropertiesConstants.BB_WVL_1) / 0.005 + 1);
        wvlFullGrid = new double[numWvl];
        for (int i = 0; i < numWvl; i++) {
            wvlFullGrid[i] = OlciSnowPropertiesConstants.BB_WVL_1 + i * 0.005;
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
        try {
            Tile[] rhoToaTiles = new Tile[OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length];
            for (int i = 0; i < OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length; i++) {
                // todo
                String rToaBandName;
                if (validL1bSourceProduct) {
                    // band names: Oa%2d_reflectance
                } else {
                    // band names: rtoa_%2d
                }
                final Band rhoToaBand = rhoToaProduct.getBand(requiredBrrBandNames[i]);
                rhoToaTiles[i] = getSourceTile(rhoToaBand, targetRectangle);
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

            Tile cloudMaskTile = null;
            if (cloudMaskProduct != null) {
                cloudMaskTile = getSourceTile(cloudMaskProduct.getRasterDataNode(cloudMaskBandName), targetRectangle);
            }

            for (int y = targetRectangle.y; y < targetRectangle.y + targetRectangle.height; y++) {
                checkForCancellation();
                for (int x = targetRectangle.x; x < targetRectangle.x + targetRectangle.width; x++) {

                    double[] rhoToa = new double[OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length];
                    for (int i = 0; i < OlciSnowPropertiesConstants.WAVELENGTH_GRID_OLCI.length; i++) {
                        rhoToa[i] = olciGains[i] * brrTiles[i].getSampleDouble(x, y);
                        rhoToa[i] = Math.max(0.0, rhoToa[i]);
                    }

                    final double brr400 = brrTiles[0].getSampleDouble(x, y);
                    final double brr560 = brrTiles[1].getSampleDouble(x, y);
                    final double brr681 = brrTiles[2].getSampleDouble(x, y);
                    final double brr709 = brrTiles[3].getSampleDouble(x, y);
                    final double brr865 = brrTiles[4].getSampleDouble(x, y);
                    final double brr1020 = brrTiles[5].getSampleDouble(x, y);

                    final boolean l1Valid = !l1FlagsTile.getSampleBit(x, y, sensor.getInvalidBit());
                    final boolean isNotCloud = cloudMaskTile == null || cloudMaskTile.getSampleDouble(x, y) != 1.0;

                    final double sza = szaTile.getSampleDouble(x, y);

                    final boolean szaIsInvalid = sza > 75.0;
                    // 20181207: do not exclude high SZA, but just raise a flag
                    final boolean pixelIsValid = l1Valid && isNotCloud;
                    final Band s3SnowFlagBand = targetProduct.getBand(S3_SNOW_FLAG_BAND_NAME);
                    targetTiles.get(s3SnowFlagBand).setSample(x, y,
                                                              OlciSnowPropertiesConstants.S3_SNOW_SZA_HIGH,
                                                              szaIsInvalid);

                    if (pixelIsValid) {
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
                            final double vza = vzaTile.getSampleDouble(x, y);
                            final double mu_0 = Math.cos(sza * MathUtils.DTOR);

                            double[][] spectralAlbedos;

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
                                                                                    brr400, sza, wvlFullGrid);

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

        final Band cloudFlagBand = targetProduct.addBand("s3snow_flags", ProductData.TYPE_INT8);
        FlagCoding flagCoding = SnowUtils.createS3SnowFlagCoding(S3_SNOW_FLAG_BAND_NAME);
        cloudFlagBand.setSampleCoding(flagCoding);
        targetProduct.getFlagCodingGroup().add(flagCoding);
        SnowUtils.setupS3SnowBitmask(targetProduct);

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
        return true;
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

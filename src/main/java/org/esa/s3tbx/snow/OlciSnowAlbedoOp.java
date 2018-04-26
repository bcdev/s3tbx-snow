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
 * Computes snow albedo quantities from OLCI L1b data products.
 *
 * @author olafd
 */
@OperatorMetadata(alias = "OLCI.SnowAlbedo",
        description = "Computes snow albedo quantities from OLCI L1b data products.",
        authors = "Alexander Kokhanovsky (EUMETSAT),  Olaf Danne (Brockmann Consult)",
        copyright = "(c) 2017, 2018 by ESA, EUMETSAT, Brockmann Consult",
        category = "Optical/Thematic Land Processing",
        version = "1.5-SNAPSHOT")

public class OlciSnowAlbedoOp extends Operator {

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
    private static final String ICE_FLAG_BAND_NAME = "ice_indicator";
    private static final String POLLUTION_MASK_BAND_NAME = "pollution_mask";
    private static final String NDSI_MASK_BAND_NAME = "ndsi_mask";
    private static final String NDSI_BAND_NAME = "ndsi";

    // AK, 20171127: no longer a user option, simple approx is best
    private SpectralAlbedoMode spectralAlbedoComputationMode = SpectralAlbedoMode.SIMPLE_APPROXIMATION;

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
            description =
                    "If selected, NDSI will be computed from 865 and 1020nm for snow identification. " +
                            "Then, only snow pixels will be considered for albedo computation.")
    private boolean considerNdsiSnowMask;

    @Parameter(defaultValue = "0.03",
            description = "NDSI threshold for snow identification",
            label = "NDSI threshold for snow identification")
    private double ndsiThresh;

    @Parameter(defaultValue = "false",
            description =
                    "If selected, polluted snow will be retrieved and specific algorithms for albedo retrieval will be used.")
    private boolean considerSnowPollution;

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
            description = "OLCI SVC gain for band 5 (default value as provided by Sentinel-3A Product Notice – " +
                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
            label = "OLCI SVC gain for band 5 (560nm)")
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

    @Parameter(defaultValue = "false",
            label = "Use new algorithm for spectral albedo (AK, 20180404)",
            description = "If selected, new algorithm for spectral albedo (provided by AK, 20180404) is used.")
    private boolean useAlgoApril2018;


    private Sensor sensor = Sensor.OLCI;
    private double[] olciGains;

    private String[] requiredBrrBandNamesAlbedo;

    private String[] requiredRadianceBandNamesPPA;     // selected spectral bands
    private String[] requiredBrrBandNamesPPA;

    private Product targetProduct;


    private Product reflProduct;

    private int width;
    private int height;

    private SolarSpectrumTable solarSpectrumTable;
    private RefractiveIndexTable refractiveIndexInterpolatedTable;


    @Override
    public void initialize() throws OperatorException {
        String[] requiredRadianceBandNamesAlbedo;

        olciGains = new double[4];
        olciGains[0] = olciGainBand1;
        olciGains[1] = olciGainBand5;
        olciGains[2] = olciGainBand17;
        olciGains[3] = olciGainBand21;
        requiredRadianceBandNamesAlbedo = new String[]{"Oa01_radiance", "Oa05_radiance", "Oa17_radiance", "Oa21_radiance"};
        requiredBrrBandNamesAlbedo = new String[]{"rBRR_01", "rBRR_05", "rBRR_17", "rBRR_21"};

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
        solarSpectrumTable = new SolarSpectrumTable();

        // interpolate input refractive indices (at 83 wavelengths) to full grid 0.3-1.02um from solar spectrum auxdata
        refractiveIndexInterpolatedTable = SnowUtils.getRefractiveIndexInterpolated(refractiveIndexTable,
                                                                                    solarSpectrumTable);

        if (cloudMaskProduct != null) {
            validateCloudMaskProduct();
        }

        createTargetProduct();
    }

    @Override
    public void computeTileStack(Map<Band, Tile> targetTiles, Rectangle targetRectangle, ProgressMonitor pm) throws OperatorException {
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
                    final boolean pixelIsValid = !l1FlagsTile.getSampleBit(x, y, sensor.getInvalidBit()) &&
                            (cloudMaskTile == null ||
                                    (cloudMaskTile != null && cloudMaskTile.getSampleDouble(x, y) != 1.0));
                    if (pixelIsValid) {
                        double[] rhoToaAlbedo = new double[requiredBrrBandNamesAlbedo.length];   // now always 01, 05, 17, 21
                        for (int i = 0; i < requiredBrrBandNamesAlbedo.length; i++) {
                            rhoToaAlbedo[i] = olciGains[i] * rhoToaTilesAlbedo[i].getSampleDouble(x, y);
                            rhoToaAlbedo[i] = Math.max(0.0, rhoToaAlbedo[i]);
                        }

                        double brr400 = rhoToaAlbedo[0];
                        final double brr560 = rhoToaAlbedo[1];
                        final double brr865 = rhoToaAlbedo[2];
                        final double brr1020 = rhoToaAlbedo[3];
                        double ndsi = Double.NaN;
                        boolean validNdsi = true;
                        if (considerNdsiSnowMask) {
                            ndsi = (brr865 - brr1020)/(brr865 + brr1020);
                            if (ndsi <=ndsiThresh || brr400 <= 0.5) {
                                validNdsi = false;
                            }
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

                            // Sigma site in
                            // subset_0_of_S3A_OL_1_EFR____20170529T004035_20170529T004335_20170529T030013_0179_018_145_1260_SVL_O_NR_002_rayleigh.dim
//                        if (x == 54 && y == 49) {
//                            System.out.println("x = " + x);
//                        }

                            double[][] spectralAlbedos;
                            double r0 = Double.NaN;
                            boolean isPollutedSnow = false;
                            if (considerSnowPollution) {
                                final double saa = saaTile.getSampleDouble(x, y);
                                final double vaa = vaaTile.getSampleDouble(x, y);
                                final double raa = SnowUtils.getRelAzi(saa, vaa);

                                r0 = OlciSnowAlbedoAlgorithm.computeR0PollutionThresh(sza, vza, raa);
                                final double pollutedSnowSeparationValue = r0 - pollutionDelta;
                                isPollutedSnow = brr400 < pollutedSnowSeparationValue;

                                if (isPollutedSnow) {
                                    final double[] pollutedSnowParams =
                                            OlciSnowAlbedoAlgorithm.computePollutedSnowParams(brr400, brr1020, sza, vza, raa);
                                    spectralAlbedos =
                                            OlciSnowAlbedoAlgorithm.computeSpectralAlbedosPolluted(pollutedSnowParams,
                                                                                                   sza, vza,
                                                                                                   useAlgoApril2018);
                                } else {
                                    spectralAlbedos =
                                            OlciSnowAlbedoAlgorithm.computeSpectralAlbedos(rhoToaAlbedo, sza, vza,
                                                                                           refWvl,
                                                                                           spectralAlbedoComputationMode,
                                                                                           useAlgoApril2018);
                                }
                            } else {
                                spectralAlbedos =
                                        OlciSnowAlbedoAlgorithm.computeSpectralAlbedos(rhoToaAlbedo, sza, vza,
                                                                                       refWvl,
                                                                                       spectralAlbedoComputationMode,
                                                                                       useAlgoApril2018);
                            }
                            final double[] spectralSphericalAlbedos = spectralAlbedos[0];
                            final double[] spectralPlanarAlbedos = spectralAlbedos[1];

                            setTargetTilesSpectralAlbedos(spectralSphericalAlbedos,
                                                          ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX, targetTiles, x, y);
                            setTargetTilesSpectralAlbedos(spectralPlanarAlbedos,
                                                          ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX, targetTiles, x, y);

                            final double refAlbedo = refWvl == 1020.0 ?
                                    spectralSphericalAlbedos[spectralSphericalAlbedos.length - 1] :
                                    spectralSphericalAlbedos[spectralSphericalAlbedos.length - 5];

                            double grainDiam = Double.NaN;
                            if (refAlbedo > 0.0) {
                                if (useAlgoApril2018) {
                                    final double l = OlciSnowAlbedoAlgorithm.computeLFromTwoWavelengths(rhoToaAlbedo,
                                                                                                        sza, vza);
                                    grainDiam = l/13.08;
                                    // todo: consider also polluted snow here
                                } else {
                                    grainDiam = OlciSnowAlbedoAlgorithm.computeGrainDiameter(refAlbedo, refWvl);
                                }
                            }
                            // todo: this is a test with 'manual' summation rather than Simpson integration.
                            // Check why Simpson is so slow!
                            final double[] broadbandPlanarAlbedo =
                                    OlciSnowAlbedoAlgorithm.computeBroadbandAlbedo(mu_0,
                                                                                   grainDiam,
                                                                                   refractiveIndexInterpolatedTable,
                                                                                   solarSpectrumTable);
                            final double[] broadbandSphericalAlbedo =
                                    OlciSnowAlbedoAlgorithm.computeBroadbandAlbedo(1.0,
                                                                                   grainDiam,
                                                                                   refractiveIndexInterpolatedTable,
                                                                                   solarSpectrumTable);

                            setTargetTilesBroadbandAlbedos(broadbandPlanarAlbedo, targetTiles, "planar", x, y);
                            setTargetTilesBroadbandAlbedos(broadbandSphericalAlbedo, targetTiles, "spherical", x, y);

                            if (computePPA && spectralAlbedoTargetBands != null) {
                                final double[] spectralPPA = OlciSnowAlbedoAlgorithm.computeSpectralPPA(rhoToaPPA, sza, vza);
                                setTargetTilesSpectralPPA(spectralPPA, PPA_SPECTRAL_OUTPUT_PREFIX, targetTiles, x, y);
                            }

                            final Band grainDiameterBand = targetProduct.getBand(GRAIN_DIAMETER_BAND_NAME);
                            final Band snowSpecificAreaBand = targetProduct.getBand(SNOW_SPECIFIC_AREA_BAND_NAME);
                            if (!Double.isNaN(grainDiam)) {
                                final double grainDiamMillim = grainDiam / 1000.0;  // in mm
                                targetTiles.get(grainDiameterBand).setSample(x, y, SnowUtils.cutTo4DecimalPlaces(grainDiamMillim));

                                final double snowSpecificArea = 6.0 / (OlciSnowAlbedoConstants.RHO_ICE * grainDiamMillim);
                                targetTiles.get(snowSpecificAreaBand).setSample(x, y, SnowUtils.cutTo7DecimalPlaces(snowSpecificArea));
                            } else {
                                targetTiles.get(grainDiameterBand).setSample(x, y, Double.NaN);
                                targetTiles.get(snowSpecificAreaBand).setSample(x, y, Double.NaN);
                            }

                            final Band iceFlagBand = targetProduct.getBand(ICE_FLAG_BAND_NAME);
                            if (brr400 > 0.0 && brr1020 > 0.0) {
                                double iceFlag = brr400 / brr1020;
                                targetTiles.get(iceFlagBand).setSample(x, y, SnowUtils.cutTo4DecimalPlaces(iceFlag));
                            } else {
                                targetTiles.get(iceFlagBand).setSample(x, y, Double.NaN);
                            }

                            if (considerSnowPollution) {
                                final Band pollutionMaskBand = targetProduct.getBand(POLLUTION_MASK_BAND_NAME);
                                if (brr400 > 0.0 && !Double.isNaN(r0)) {
                                    int pollutionMask = isPollutedSnow ? 1 : 0;
                                    targetTiles.get(pollutionMaskBand).setSample(x, y, pollutionMask);
                                } else {
                                    targetTiles.get(pollutionMaskBand).setSample(x, y, 0);
                                }
                            }
                            if (considerNdsiSnowMask) {
                                final Band ndsiMaskBand = targetProduct.getBand(NDSI_MASK_BAND_NAME);
                                final Band ndsiBand = targetProduct.getBand(NDSI_BAND_NAME);
                                int ndsiMask = validNdsi ? 1 : 0;
                                targetTiles.get(ndsiMaskBand).setSample(x, y, ndsiMask);
                                targetTiles.get(ndsiBand).setSample(x, y, ndsi);
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
        targetProduct.addBand(ICE_FLAG_BAND_NAME, ProductData.TYPE_FLOAT32);

        if (considerSnowPollution) {
            targetProduct.addBand(POLLUTION_MASK_BAND_NAME, ProductData.TYPE_INT16);
        }

        if (considerNdsiSnowMask) {
            targetProduct.addBand(NDSI_MASK_BAND_NAME, ProductData.TYPE_INT16);
            targetProduct.addBand(NDSI_BAND_NAME, ProductData.TYPE_FLOAT32);
        }

        if (spectralAlbedoTargetBands != null && spectralAlbedoTargetBands.length > 0) {
            for (final String targetBand : spectralAlbedoTargetBands) {
                final int spectralBandIndex = Integer.parseInt(targetBand.substring(2, 4));
                final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[spectralBandIndex - 1] * 1000.0;
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
        targetProduct.getBand(SNOW_SPECIFIC_AREA_BAND_NAME).setUnit("mm");

        if (considerSnowPollution) {
            targetProduct.getBand(POLLUTION_MASK_BAND_NAME).setNoDataValueUsed(false);
        }

        if (cloudMaskProduct != null) {
            ProductUtils.copyBand(cloudMaskBandName, cloudMaskProduct, targetProduct, true);
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
                final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[spectralBandIndex - 1] * 1000.0;
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
                final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[spectralBandIndex - 1] * 1000.0;
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
        for (int i = 0; i < OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI.length; i++) {
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
                                                        "'rBRR_*' for 440nm, 865 or 1020nm, and for all manually " +
                                                        "selected wavelengths.");
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
            super(OlciSnowAlbedoOp.class);
        }
    }
}

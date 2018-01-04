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
        copyright = "(c) 2017 by EUMETSAT, Brockmann Consult",
        category = "Optical/Thematic Land Processing",
        version = "1.0")

public class OlciSnowAlbedoOp extends Operator {

    private static final String ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX = "albedo_spectral_spherical_";
    private static final String ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX = "albedo_spectral_planar_";

    private static final String ALBEDO_BB_SPHERICAL_OUTPUT_PREFIX = "albedo_bb_spherical_";
    private static final String ALBEDO_BB_PLANAR_OUTPUT_PREFIX = "albedo_bb_planar_";
    private static final String ALBEDO_BB_OUTPUT_PREFIX = "albedo_bb_";

    private static final String ALBEDO_BROADBAND_VIS_SPHERICAL_BAND_NAME = "albedo_bb_spherical_vis";
    private static final String ALBEDO_BROADBAND_NIR_SPHERICAL_BAND_NAME = "albedo_bb_spherical_nir";
    private static final String ALBEDO_BROADBAND_SW_SPHERICAL_BAND_NAME = "albedo_bb_spherical_sw";

    private static final String ALBEDO_BROADBAND_VIS_PLANAR_BAND_NAME = "albedo_bb_planar_vis";
    private static final String ALBEDO_BROADBAND_NIR_PLANAR_BAND_NAME = "albedo_bb_planar_nir";
    private static final String ALBEDO_BROADBAND_SW_PLANAR_BAND_NAME = "albedo_bb_planar_sw";

    private static final String ALBEDO_BROADBAND_SPHERICAL_PREFIX = "albedo_bb_spherical_";
    private static final String ALBEDO_BROADBAND_PLANAR_PREFIX = "albedo_bb_planar_";
    private static final String[] ALBEDO_BROADBAND_SUFFIXES = {"vis", "nir", "sw"};

    private static final String GRAIN_DIAMETER_BAND_NAME = "grain_diameter";

//    @Parameter(label = "Spectral albedo computation mode", defaultValue = "SIMPLE_APPROXIMATION",
//            description = "Spectral albedo computation mode (i.e. suitable way of curve fitting)")
//    private SpectralAlbedoMode spectralAlbedoComputationMode;

    // AK, 20171127: no longer a user option, simple approx is best
    private SpectralAlbedoMode spectralAlbedoComputationMode = SpectralAlbedoMode.SIMPLE_APPROXIMATION;

    @Parameter(description = "The OLCI wavelengths for spectral spherical and planar albedos which will " +
            "be written to the target product.",
            label = "Select OLCI wavelengths for spectral albedos",
            valueSet = {
                    "Oa01 (400 nm)", "Oa02 (412.5 nm)", "Oa03 (442.5 nm)", "Oa04 (490 nm)", "Oa05 (510 nm)",
                    "Oa06 (560 nm)", "Oa07 (620 nm)", "Oa08 (665 nm)", "Oa09 (673.75 nm)", "Oa10 (681.25 nm)",
                    "Oa11 (708.75 nm)", "Oa12 (753.75 nm)", "Oa13 (761.25 nm)", "Oa14 (764.375 nm)", "Oa15 (767.5 nm)",
                    "Oa16 (778.75 nm)", "Oa17 (865 nm)", "Oa18 (885 nm)", "Oa19 (900 nm)", "Oa20 (940 nm)",
                    "Oa21 (1020 nm)"
            },
            defaultValue = "")
    String[] spectralAlbedoTargetBands;

    @Parameter(defaultValue = "false",
            description = "If selected, Rayleigh corrected reflectances are written to target product")
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

    private Sensor sensor = Sensor.OLCI;
    private double[] olciGains;

    private Product targetProduct;


    private int reflType;

    private Product reflProduct;

    RefractiveIndexTable refractiveIndexTable;
    SolarSpectrumTable solarSpectrumTable;
    private RefractiveIndexTable refractiveIndexInterpolatedTable;


    @Override
    public void initialize() throws OperatorException {
        checkSensorType(sourceProduct, sensor);

        // for SIMPLE_APPROXIMATION we need only OLCI gains for bands (1, 21) or (1, 17)
        olciGains = new double[2];
        olciGains[0] = olciGainBand1;
        olciGains[1] = olciGainBand21;

        if (refWvl == 865.0) {
            sensor.setRequiredRadianceBandNames(new String[]{"Oa01_radiance", "Oa17_radiance"});
            sensor.setRequiredBrrBandNames(new String[]{"rBRR_01", "rBRR_17"});
            olciGains[1] = olciGainBand17;
        }

        if (isValidRayleighCorrectedSourceProduct(sourceProduct, sensor)) {
            reflProduct = sourceProduct;
            reflType = SensorConstants.REFL_TYPE_BRR;
        } else if (isValidL1bSourceProduct(sourceProduct, sensor)) {
            // apply Rayleigh correction
            // todo: compute Rayleigh correction for all selected spectral albedo/ppa wavelengths + bands 1, 17, 21
            // make separate brr subsets for albedo (only 1 + 17/21) and ppa (the selected bands)
            RayleighCorrectionOp rayleighCorrectionOp = new RayleighCorrectionOp();
            rayleighCorrectionOp.setSourceProduct(sourceProduct);
            rayleighCorrectionOp.setParameterDefaultValues();
            rayleighCorrectionOp.setParameter("computeTaur", false);
            // for SIMPLE_APPROXIMATION we need only OLCI BRR bands 1, 21
            rayleighCorrectionOp.setParameter("sourceBandNames", sensor.getRequiredRadianceBandNames());
            reflProduct = rayleighCorrectionOp.getTargetProduct();
            reflType = SensorConstants.REFL_TYPE_TOA;
        } else {
            throw new OperatorException
                    ("Input product not supported - must be " + Sensor.OLCI.getName() +
                             " L1b or Rayleigh corrected BRR product");
        }

        // read auxiliary data:
        refractiveIndexTable = new RefractiveIndexTable();
        solarSpectrumTable = new SolarSpectrumTable();

        // interpolate input refractive indices (at 83 wavelengths) to full grid 0.3-1.02um from solar spectrum auxdata
        refractiveIndexInterpolatedTable = SnowUtils.getRefractiveIndexInterpolated(refractiveIndexTable,
                                                                                    solarSpectrumTable);

        createTargetProduct();
    }

    @Override
    public void computeTileStack(Map<Band, Tile> targetTiles, Rectangle targetRectangle, ProgressMonitor pm) throws OperatorException {
        try {
            final String[] reflBandNames = SnowUtils.getReflectanceTypeBandNames(sensor, reflType);
            Tile[] rhoToaTiles = new Tile[reflBandNames.length];
            for (int i = 0; i < reflBandNames.length; i++) {
                final Band rhoToaBand = reflProduct.getBand(reflBandNames[i]);
                rhoToaTiles[i] = getSourceTile(rhoToaBand, targetRectangle);
            }

            Tile szaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getSzaName()), targetRectangle);
            Tile vzaTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getVzaName()), targetRectangle);
            Tile l1FlagsTile = getSourceTile(sourceProduct.getRasterDataNode(sensor.getL1bFlagsName()), targetRectangle);

            for (int y = targetRectangle.y; y < targetRectangle.y + targetRectangle.height; y++) {
                checkForCancellation();
                for (int x = targetRectangle.x; x < targetRectangle.x + targetRectangle.width; x++) {
                    if (!l1FlagsTile.getSampleBit(x, y, sensor.getInvalidBit())) {

                        // we compute snow albedo over land only
                        double[] rhoToa = new double[reflBandNames.length];
                        for (int i = 0; i < reflBandNames.length; i++) {
                            rhoToa[i] = olciGains[i] * rhoToaTiles[i].getSampleDouble(x, y);
                            rhoToa[i] = Math.max(0.0, rhoToa[i]);
                        }

                        final double vza = vzaTile.getSampleDouble(x, y);
                        final double sza = szaTile.getSampleDouble(x, y);
                        final double mu_0 = Math.cos(sza * MathUtils.DTOR);

                        // Sigma site in subset_0_of_S3A_OL_1_EFR____20170529T004035_20170529T004335_20170529T030013_0179_018_145_1260_SVL_O_NR_002_rayleigh.dim
//                        if (x == 54 && y == 49) {
//                            System.out.println("x = " + x);
//                        }

                        // default is actually SIMPLE_APPROXIMATION (latest approach from AK, 20171120):
                        final double[][] sphericalAlbedos =
                                OlciSnowAlbedoAlgorithm.computeSphericalAlbedos(rhoToa, sza, vza,
                                                                                refWvl,
                                                                                spectralAlbedoComputationMode);
                        final double[] spectralSphericalAlbedos = sphericalAlbedos[0];
                        final double[] spectralPlanarAlbedos = sphericalAlbedos[1];

                        setTargetTilesSpectralAlbedos(spectralSphericalAlbedos,
                                                      ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX, targetTiles, x, y);
                        setTargetTilesSpectralAlbedos(spectralPlanarAlbedos,
                                                      ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX, targetTiles, x, y);

                        final double refAlbedo = refWvl == 1020.0 ?
                                spectralSphericalAlbedos[spectralSphericalAlbedos.length - 1] :
                                spectralSphericalAlbedos[spectralSphericalAlbedos.length - 5];
                        final double grainDiam = OlciSnowAlbedoAlgorithm.computeGrainDiameter(refAlbedo, refWvl);
                        // todo: this is a test with 'manual' summation rather than Simpson integration.
                        // Check why Simpson is so slow!
                        final double[] broadbandPlanarAlbedo =
                                OlciSnowAlbedoAlgorithm.computeBroadbandAlbedo_test(mu_0,
                                                                                    grainDiam,
                                                                                    refractiveIndexInterpolatedTable,
                                                                                    solarSpectrumTable);
                        final double[] broadbandSphericalAlbedo =
                                OlciSnowAlbedoAlgorithm.computeBroadbandAlbedo_test(1.0,
                                                                                    grainDiam,
                                                                                    refractiveIndexInterpolatedTable,
                                                                                    solarSpectrumTable);

                        setTargetTilesBroadbandAlbedos(broadbandPlanarAlbedo, targetTiles, "planar", x, y);
                        setTargetTilesBroadbandAlbedos(broadbandSphericalAlbedo, targetTiles, "spherical", x, y);

                        final Band grainDiameterBand = targetProduct.getBand(GRAIN_DIAMETER_BAND_NAME);
                        targetTiles.get(grainDiameterBand).setSample(x, y, grainDiam / 1000.0);  // in mm
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
        targetProduct = new Product(sourceProduct.getName(), sourceProduct.getProductType(),
                                    sourceProduct.getSceneRasterWidth(), sourceProduct.getSceneRasterHeight());

        if (copyReflectanceBands) {
            for (Band band : reflProduct.getBands()) {
                if (band.getName().startsWith(SensorConstants.OLCI_BRR_BAND_PREFIX)) {
                    ProductUtils.copyBand(band.getName(), reflProduct, targetProduct, true);
                    ProductUtils.copyRasterDataNodeProperties(band, targetProduct.getBand(band.getName()));
                }
            }
        }

        targetProduct.addBand(ALBEDO_BROADBAND_VIS_SPHERICAL_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_NIR_SPHERICAL_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_SW_SPHERICAL_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_VIS_PLANAR_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_NIR_PLANAR_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_SW_PLANAR_BAND_NAME, ProductData.TYPE_FLOAT32);

        targetProduct.addBand(GRAIN_DIAMETER_BAND_NAME, ProductData.TYPE_FLOAT32);

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
            }
        }

        for (Band band : targetProduct.getBands()) {
            band.setUnit("dl");
            band.setNoDataValue(Float.NaN);
            band.setNoDataValueUsed(true);
        }
        targetProduct.getBand(GRAIN_DIAMETER_BAND_NAME).setUnit("mm");

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

    private void setTargetTilesSpectralAlbedos(double[] spectralAlbedos, String prefix, Map<Band, Tile> targetTiles, int x, int y) {
        if (spectralAlbedoTargetBands != null && spectralAlbedoTargetBands.length > 0) {
            for (final String targetBand : spectralAlbedoTargetBands) {
                final int spectralBandIndex = Integer.parseInt(targetBand.substring(2, 4));
                final double wvl = OlciSnowAlbedoConstants.WAVELENGTH_GRID_OLCI[spectralBandIndex - 1] * 1000.0;
                final Band spectralAlbedoBand = targetProduct.getBand(prefix + (int) wvl);
                targetTiles.get(spectralAlbedoBand).setSample(x, y, spectralAlbedos[spectralBandIndex - 1]);
            }
        }
    }

    private void setTargetTilesBroadbandAlbedos(double[] bbAlbedos, Map<Band, Tile> targetTiles,
                                                String bbMode, int x, int y) {
        int index = 0;
        for (final String bbSuffix : ALBEDO_BROADBAND_SUFFIXES) {
            final Band targetBand = targetProduct.getBand(ALBEDO_BB_OUTPUT_PREFIX + bbMode + "_" + bbSuffix);
            targetTiles.get(targetBand).setSample(x, y, bbAlbedos[index++]);
        }
    }

    private void setTargetTilesInvalid(Map<Band, Tile> targetTiles, int x, int y) {
        for (Tile tile : targetTiles.values()) {
            tile.setSample(x, y, Float.NaN);
        }
    }

    private static void checkSensorType(Product sourceProduct, Sensor sensor) {
        boolean isOlci = isValidL1bSourceProduct(sourceProduct, sensor);
        if (!isOlci) {
            isOlci = isValidRayleighCorrectedSourceProduct(sourceProduct, sensor);
            if (!isOlci) {
                throw new OperatorException("Source product not applicable to this operator.\n" +
                                                    "Only OLCI is currently supported");
            }
        }
    }

    private static boolean isValidL1bSourceProduct(Product sourceProduct, Sensor sensor) {
        for (String bandName : sensor.getRequiredRadianceBandNames()) {
            if (!sourceProduct.containsBand(bandName)) {
                return false;
            }
        }
        return true;
    }

    private static boolean isValidRayleighCorrectedSourceProduct(Product sourceProduct, Sensor sensor) {
        for (String bandName : sensor.getRequiredBrrBandNames()) {
            if (!sourceProduct.containsBand(bandName)) {
                return false;
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

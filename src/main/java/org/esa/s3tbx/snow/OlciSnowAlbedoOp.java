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
    private static final String ALBEDO_BROADBAND_SPHERICAL_BAND_NAME = "albedo_broadband_spherical";
    private static final String ALBEDO_BROADBAND_PLANAR_BAND_NAME = "albedo_broadband_planar";
    private static final String GRAIN_DIAMETER_BAND_NAME = "grain_diameter";


    @Parameter(label = "Spectral albedo computation mode", defaultValue = "EXPONENTIAL_4PARAM_FIT",
            description = "Spectral albedo computation mode (i.e. suitable way of curve fitting)")
    private SpectralAlbedoMode spectralAlbedoComputationMode;

    @Parameter(defaultValue = "false",
            description = "If selected, Rayleigh corrected reflectances are written to target product")
    private boolean copyReflectanceBands;

    // we need gains only for OLCI BRR bands 1, 2, 3, 4, 5, 12, 17, 21
    @Parameter(defaultValue = "0.9798, 0.9718, 0.9747, 0.9781, 0.9827, 1.003, 1.0, 0.914",
            description = "OLCI SVC gains for bands 1-5, 12, 17, 21 (default values as provided by Sentinel-3A Product Notice – " +
                    "OLCI Level-2 Ocean Colour, July 5th, 2017",
            label = "OLCI SVC gains (bands 1-5, 12, 17, 21)")
    private double[] olciGains;


    @SourceProduct(description = "OLCI L1b or Rayleigh corrected product",
            label = "OLCI L1b or Rayleigh corrected product")
    private Product sourceProduct;

    private Sensor sensor = Sensor.OLCI;

    private Product targetProduct;


    private int reflType;

    private Product reflProduct;


    @Override
    public void initialize() throws OperatorException {
        checkSensorType(sourceProduct);

        if (isValidRayleighCorrectedSourceProduct(sourceProduct, Sensor.OLCI)) {
            reflProduct = sourceProduct;
            reflType = SensorConstants.REFL_TYPE_BRR;
        } else if (isValidL1bSourceProduct(sourceProduct, Sensor.OLCI)) {

            // apply Rayleigh correction
            RayleighCorrectionOp rayleighCorrectionOp = new RayleighCorrectionOp();
            rayleighCorrectionOp.setSourceProduct(sourceProduct);
            rayleighCorrectionOp.setParameterDefaultValues();
            rayleighCorrectionOp.setParameter("computeTaur", false);
            // we need OLCI BRR bands 1, 2, 3, 4, 5, 12, 17, 21
            rayleighCorrectionOp.setParameter("sourceBandNames", Sensor.OLCI.getRequiredRadianceBandNames());
            reflProduct = rayleighCorrectionOp.getTargetProduct();
            reflType = SensorConstants.REFL_TYPE_TOA;
        } else {
            throw new OperatorException
                    ("Input product not supported - must be " + Sensor.OLCI.getName() +
                             " L1b or Rayleigh corrected BRR product");
        }

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

            Tile szaTile = getSourceTile(sourceProduct.getRasterDataNode(Sensor.OLCI.getSzaName()), targetRectangle);
            Tile vzaTile = getSourceTile(sourceProduct.getRasterDataNode(Sensor.OLCI.getVzaName()), targetRectangle);
            Tile l1FlagsTile = getSourceTile(sourceProduct.getRasterDataNode(Sensor.OLCI.getL1bFlagsName()), targetRectangle);

            for (int y = targetRectangle.y; y < targetRectangle.y + targetRectangle.height; y++) {
                checkForCancellation();
                for (int x = targetRectangle.x; x < targetRectangle.x + targetRectangle.width; x++) {
                    if (!l1FlagsTile.getSampleBit(x, y, Sensor.OLCI.getInvalidBit())) {

                        // we compute snow albedo over land only
                        double[] rhoToa = new double[reflBandNames.length];
                        for (int i = 0; i < reflBandNames.length; i++) {
                            rhoToa[i] = olciGains[i] * rhoToaTiles[i].getSampleDouble(x, y);
                            rhoToa[i] = Math.max(0.0, rhoToa[i]);
                        }

                        final double sza = szaTile.getSampleDouble(x, y);
                        final double vza = vzaTile.getSampleDouble(x, y);

                        if (x == 54 && y == 49) {
                            System.out.println("x = " + x);
                        }

                        // actually done with latest approach from AK, 20171120:
                        final double[][] sphericalAlbedos =
                                OlciSnowAlbedoAlgorithm.computeSphericalAlbedos_nov20(rhoToa, sza, vza);
                        final double[] spectralSphericalAlbedos = sphericalAlbedos[0];
                        final double[] spectralPlanarAlbedos = sphericalAlbedos[1];

                        // actually done with latest approach from AK, 20170929
//                        final double[] spectralSphericalAlbedos =
//                                OlciSnowAlbedoAlgorithm.computeSpectralSphericalAlbedos(rhoToa,
//                                                                                        sza, vza,
//                                                                                        spectralAlbedoComputationMode);

                        // This is a test using LMA fitting library, which helped to find the negative brr input,
                        // which obviously leads to infinite iterations in the apache-commons polynominal
                        // or sigmoidal fits. todo: check and identify that more clearly
//                            final double[] spectralSphericalAlbedos =
//                                    OlciSnowAlbedoAlgorithm.computeSpectralSphericalAlbedos_test_lma(rhoToa, sza, vza);

                        setTargetTilesSpectralAlbedos(spectralSphericalAlbedos,
                                                      ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX, targetTiles, x, y);

//                        final double[] spectralPlanarAlbedos =
//                                OlciSnowAlbedoAlgorithm.computePlanarFromSphericalAlbedos(spectralSphericalAlbedos, sza);
//                        setTargetTilesSpectralAlbedos(spectralPlanarAlbedos,
//                                                      ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX, targetTiles, x, y);

                        final OlciSnowAlbedoAlgorithm.SphericalBroadbandAlbedo sbbaTerms =
                                OlciSnowAlbedoAlgorithm.computeSphericalBroadbandAlbedoTerms(spectralSphericalAlbedos);

                        final double sbba = sbbaTerms.getR_b1() + sbbaTerms.getR_b2();
                        final double planarBroadbandAlbedo =
                                OlciSnowAlbedoAlgorithm.computePlanarFromSphericalAlbedo(sbba, sza);
                        final Band sphericalBBABand = targetProduct.getBand(ALBEDO_BROADBAND_SPHERICAL_BAND_NAME);
                        final Band planarBBABand = targetProduct.getBand(ALBEDO_BROADBAND_PLANAR_BAND_NAME);
                        final Band grainDiameterBand = targetProduct.getBand(GRAIN_DIAMETER_BAND_NAME);
                        targetTiles.get(sphericalBBABand).setSample(x, y, sbba);
                        targetTiles.get(planarBBABand).setSample(x, y, planarBroadbandAlbedo);
                        final double grainDiameterInMillimeter = sbbaTerms.getGrainDiameter() * 0.001;
                        targetTiles.get(grainDiameterBand).setSample(x, y, grainDiameterInMillimeter);
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
//                if (band.getName().endsWith(SensorConstants.OLCI_REFL_BAND_SUFFIX)) {
                    ProductUtils.copyBand(band.getName(), reflProduct, targetProduct, true);
                    ProductUtils.copyRasterDataNodeProperties(band, targetProduct.getBand(band.getName()));
                }
            }
        }

        targetProduct.addBand(ALBEDO_BROADBAND_SPHERICAL_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(ALBEDO_BROADBAND_PLANAR_BAND_NAME, ProductData.TYPE_FLOAT32);
        targetProduct.addBand(GRAIN_DIAMETER_BAND_NAME, ProductData.TYPE_FLOAT32);

        for (int i = 0; i < OlciSnowAlbedoConstants.SPECTRAL_ALBEDO_OUTPUT_WAVELENGTHS.length; i++) {
            final int wvl = OlciSnowAlbedoConstants.SPECTRAL_ALBEDO_OUTPUT_WAVELENGTHS[i];
            final Band sphericalBand =
                    targetProduct.addBand(ALBEDO_SPECTRAL_SPHERICAL_OUTPUT_PREFIX + wvl, ProductData.TYPE_FLOAT32);
            sphericalBand.setSpectralWavelength(wvl);
            sphericalBand.setSpectralBandIndex(OlciSnowAlbedoConstants.SPECTRAL_ALBEDO_OUTPUT_WAVELENGTH_INDICES[i]);
            final Band planarBand =
                    targetProduct.addBand(ALBEDO_SPECTRAL_PLANAR_OUTPUT_PREFIX + wvl, ProductData.TYPE_FLOAT32);
            planarBand.setSpectralWavelength(wvl);
            planarBand.setSpectralBandIndex(OlciSnowAlbedoConstants.SPECTRAL_ALBEDO_OUTPUT_WAVELENGTH_INDICES[i]);
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
//        ProductUtils.copyTiePointGrids(sourceProduct, targetProduct);
        ProductUtils.copyMasks(sourceProduct, targetProduct);
        ProductUtils.copyFlagBands(sourceProduct, targetProduct, true);
        ProductUtils.copyGeoCoding(sourceProduct, targetProduct);
        targetProduct.setAutoGrouping("rBRR:albedo_spectral_spherical:albedo_spectral_planar");

        setTargetProduct(targetProduct);
    }

    private void setTargetTilesSpectralAlbedos(double[] spectralAlbedos, String prefix, Map<Band, Tile> targetTiles, int x, int y) {
        for (int i = 0; i < OlciSnowAlbedoConstants.SPECTRAL_ALBEDO_OUTPUT_WAVELENGTHS.length; i++) {
            final int wvl = OlciSnowAlbedoConstants.SPECTRAL_ALBEDO_OUTPUT_WAVELENGTHS[i];
            final Band spectralAlbedoBand = targetProduct.getBand(prefix + wvl);
            int spectralIndex = OlciSnowAlbedoConstants.SPECTRAL_ALBEDO_OUTPUT_WAVELENGTH_INDICES[i];
            targetTiles.get(spectralAlbedoBand).setSample(x, y, spectralAlbedos[spectralIndex]);
        }
    }

    private void setTargetTilesInvalid(Map<Band, Tile> targetTiles, int x, int y) {
        for (Tile tile : targetTiles.values()) {
            tile.setSample(x, y, Float.NaN);
        }
    }

    private static void checkSensorType(Product sourceProduct) {
        boolean isOlci = isValidL1bSourceProduct(sourceProduct, Sensor.OLCI);
        if (!isOlci) {
            isOlci = isValidRayleighCorrectedSourceProduct(sourceProduct, Sensor.OLCI);
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

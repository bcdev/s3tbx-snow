package org.esa.s3tbx.snow;

import org.esa.snap.core.gpf.OperatorException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

/**
 * todo: add comment
 * To change this template use File | Settings | File Templates.
 * Date: 04.12.2017
 * Time: 12:30
 *
 * @author olafd
 */
public class SnowAuxdata {

    public static final int SOLAR_SPECTRUM_TABLE_LENGTH = 7923;
    public static final String SOLAR_SPECTRUM_FILE_NAME = "solar_spectrum.txt";

    public static final int REFRACTIVE_INDEX_TABLE_LENGTH = 162;
    public static final String REFRACTIVE_INDEX_FILE_NAME = "refractive_index.txt";

    private static SnowAuxdata instance;
    private SolarSpectrumTable solarSpectrumTable;

    public static SnowAuxdata getInstance() {
        if (instance == null) {
            instance = new SnowAuxdata();
        }

        return instance;
    }


    public SolarSpectrumTable createSolarSpectrumTable() throws IOException {
        final InputStream inputStream = getClass().getResourceAsStream(SOLAR_SPECTRUM_FILE_NAME);
        SolarSpectrumTable spectrumTable = new SolarSpectrumTable(SOLAR_SPECTRUM_TABLE_LENGTH);

        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
        StringTokenizer st;
        try {
            int i = 0;
            String line;
            while ((line = bufferedReader.readLine()) != null && i < SOLAR_SPECTRUM_TABLE_LENGTH) {
                line = line.trim();
                st = new StringTokenizer(line, ";", false);

                if (st.hasMoreTokens()) {
                    // x (whatever that is)
                    spectrumTable.setWvl(i, Double.parseDouble(st.nextToken()));
                }
                if (st.hasMoreTokens()) {
                    // y
                    spectrumTable.setSolarSpectrum(i, Double.parseDouble(st.nextToken()));
                }
                i++;
            }
        } catch (IOException | NumberFormatException e) {
            throw new OperatorException("Failed to load SolarSpectrumTable: \n" + e.getMessage(), e);
        } finally {
            inputStream.close();
        }
        return spectrumTable;
    }

    public RefractiveIndexTable createRefractiveIndexTable() throws IOException {
        final InputStream inputStream = getClass().getResourceAsStream(REFRACTIVE_INDEX_FILE_NAME);
        RefractiveIndexTable refractiveIndexTable = new RefractiveIndexTable(REFRACTIVE_INDEX_TABLE_LENGTH);

        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
        StringTokenizer st;
        try {
            int i = 0;
            String line;
            while ((line = bufferedReader.readLine()) != null && i < REFRACTIVE_INDEX_TABLE_LENGTH) {
                line = line.trim();
                st = new StringTokenizer(line, ";", false);

                if (st.hasMoreTokens()) {
                    // x (whatever that is)
                    refractiveIndexTable.setWvl(i, Double.parseDouble(st.nextToken()));
                }
                if (st.hasMoreTokens()) {
                    // real part
                    refractiveIndexTable.setRefractiveIndexReal(i, Double.parseDouble(st.nextToken()));
                }
                if (st.hasMoreTokens()) {
                    // imaginary part
                    refractiveIndexTable.setRefractiveIndexImag(i, Double.parseDouble(st.nextToken()));
                }
                i++;
            }
        } catch (IOException | NumberFormatException e) {
            throw new OperatorException("Failed to load SolarSpectrumTable: \n" + e.getMessage(), e);
        } finally {
            inputStream.close();
        }
        return refractiveIndexTable;
    }

}

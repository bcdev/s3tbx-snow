package org.esa.s3tbx.snow;

import org.esa.snap.core.gpf.OperatorException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

/**
 * Holder for new Solar Spectrum Table as provided from text file.
 * This table contains fluxes also for SZA = 0, 15, 30, 45 and 75deg and allows interpolation rather than
 * assuming a constant SZA of 60deg.
 *
 * @author olafd
 */
public class SolarSpectrumExtendedTable {

    private static final int SOLAR_SPECTRUM_TABLE_DEFAULT_LENGTH = 411;
    private static final String SOLAR_SPECTRUM_DEFAULT_FILE_NAME = "final_table_fluxes_angles.txt";
//    private static final String SOLAR_SPECTRUM_DEFAULT_FILE_NAME = "final_table_fluxes_angles_oct2018.txt";

    private double[] wvl;
    private double[][] solarSpectrum;
    private int length;
    private String filename;

    SolarSpectrumExtendedTable() {
        length = SOLAR_SPECTRUM_TABLE_DEFAULT_LENGTH;
        filename = SOLAR_SPECTRUM_DEFAULT_FILE_NAME;
        wvl = new double[length];
        solarSpectrum = new double[6][length];  // for SZA = 0, 15, 30, 45, 60, 75deg
//        solarSpectrum = new double[89][length];  // for SZA = 0, ... , 88deg
        readTableFromFile();
    }

    void readTableFromFile() {
        final InputStream inputStream = getClass().getResourceAsStream(filename);

        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
        StringTokenizer st;
        try {
            int rowIndex = 0;
            String line;
            while ((line = bufferedReader.readLine()) != null && rowIndex < length) {
                line = line.trim();
                st = new StringTokenizer(line, ";", false);

                if (st.hasMoreTokens()) {
                    setWvl(rowIndex, Double.parseDouble(st.nextToken()));  // this table is in microns!
                }

                int szaIndex = 0;
                while (st.hasMoreTokens()) {
                    setSolarSpectrum(szaIndex++, rowIndex, Double.parseDouble(st.nextToken()));
                }

                rowIndex++;
            }
        } catch (IOException | NumberFormatException e) {
            throw new OperatorException("Failed to load SolarSpectrumTable: \n" + e.getMessage(), e);
        } finally {
            try {
                inputStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    private void setWvl(int rowIndex, double wvl) {
        this.wvl[rowIndex] = wvl;
    }

    public double getWvl(int rowIndex) {
        return wvl[rowIndex];
    }

    public double[] getWvl() {
        return wvl;
    }

    public double[][] getSolarSpectrum() {
        return solarSpectrum;
    }

    public double getSolarSpectrum(int szaIndex, int rowIndex) {
        return solarSpectrum[szaIndex][rowIndex];
    }

    private void setSolarSpectrum(int szaIndex, int rowIndex, double spectrum) {
        this.solarSpectrum[szaIndex][rowIndex] = spectrum;
    }

}

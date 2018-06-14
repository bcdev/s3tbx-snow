package org.esa.s3tbx.snow;

import org.esa.snap.core.gpf.OperatorException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

/**
 * Holder for Solar Spectrum Table as provided from text file.
 *
 * @author olafd
 */
public class SolarSpectrumTable {

//    public static final int SOLAR_SPECTRUM_TABLE_DEFAULT_LENGTH = 7923;
//    public static final String SOLAR_SPECTRUM_DEFAULT_FILE_NAME = "solar_spectrum.txt";
    private static final int SOLAR_SPECTRUM_TABLE_DEFAULT_LENGTH = 411;
    private static final String SOLAR_SPECTRUM_DEFAULT_FILE_NAME = "solar_spectrum_2.txt";

    private double[] wvl;
    private double[] solarSpectrum;
    private int length;
    private String filename;

    SolarSpectrumTable() {
        length = SOLAR_SPECTRUM_TABLE_DEFAULT_LENGTH;
        filename = SOLAR_SPECTRUM_DEFAULT_FILE_NAME;
        wvl = new double[length];
        solarSpectrum = new double[length];
        readTableFromFile();
    }

    private void readTableFromFile() {
        final InputStream inputStream = getClass().getResourceAsStream(filename);

        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
        StringTokenizer st;
        try {
            int i = 0;
            String line;
            while ((line = bufferedReader.readLine()) != null && i < length) {
                line = line.trim();
                st = new StringTokenizer(line, ";", false);

                if (st.hasMoreTokens()) {
                    String nextToken = st.nextToken();
                    setWvl(i, Double.parseDouble(nextToken)/1000.0);  // in microns! For the moment we assume table has nm.
                }
                if (st.hasMoreTokens()) {
                    setSolarSpectrum(i, Double.parseDouble(st.nextToken()));
                }
                i++;
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

    private void setWvl(int index, double wvl) {
        this.wvl[index] = wvl;
    }

    public double getWvl(int index) {
        return wvl[index];
    }

    public double[] getSolarSpectrum() {
        return solarSpectrum;
    }

    public double getSolarSpectrum(int index) {
        return solarSpectrum[index];
    }

    private void setSolarSpectrum(int index, double spectrum) {
        this.solarSpectrum[index] = spectrum;
    }

}

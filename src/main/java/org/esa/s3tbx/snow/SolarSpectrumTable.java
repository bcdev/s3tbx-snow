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
 * Date: 05.12.2017
 * Time: 17:06
 *
 * @author olafd
 */
public class SolarSpectrumTable {

//    public static final int SOLAR_SPECTRUM_TABLE_DEFAULT_LENGTH = 7923;
//    public static final String SOLAR_SPECTRUM_DEFAULT_FILE_NAME = "solar_spectrum.txt";
    public static final int SOLAR_SPECTRUM_TABLE_DEFAULT_LENGTH = 411;
    public static final String SOLAR_SPECTRUM_DEFAULT_FILE_NAME = "solar_spectrum_2.txt";

    private double[] wvl;
    private double[] solarSpectrum;
    private int length;
    private String filename;

    public SolarSpectrumTable() {
        length = SOLAR_SPECTRUM_TABLE_DEFAULT_LENGTH;
        filename = SOLAR_SPECTRUM_DEFAULT_FILE_NAME;
        wvl = new double[length];
        solarSpectrum = new double[length];
        readTableFromFile();
    }

    public SolarSpectrumTable(int length) {
        this.length = length;
    }

    public SolarSpectrumTable(int length, String filename) {
        wvl = new double[length];
        solarSpectrum = new double[length];
        length = length;
        filename = filename;
    }

    public void readTableFromFile() {
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
                    setWvl(i, Double.parseDouble(st.nextToken())/1000.0);  // in microns! For the moment we assume table has nm.
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

    public void setWvl(double[] wvl) {
        this.wvl = wvl;
    }

    public void setSolarSpectrum(double[] solarSpectrum) {
        this.solarSpectrum = solarSpectrum;
    }

    public void setWvl(int index, double wvl) {
        this.wvl[index] = wvl;
    }

    public double getWvl(int index) {
        return wvl[index];
    }

    public double[] getWvl() {
        return wvl;
    }

    public double[] getSolarSpectrum() {
        return solarSpectrum;
    }

    public double getSolarSpectrum(int index) {
        return solarSpectrum[index];
    }

    public void setSolarSpectrum(int index, double spectrum) {
        this.solarSpectrum[index] = spectrum;
    }

}

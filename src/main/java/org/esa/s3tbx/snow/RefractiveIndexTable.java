package org.esa.s3tbx.snow;

import org.esa.snap.core.gpf.OperatorException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

/**
 * Holder for Refractive Index Table as provided from text file.
 *
 * @author olafd
 */
public class RefractiveIndexTable {
    public static final int REFRACTIVE_INDEX_TABLE_LENGTH = 162;
    public static final String REFRACTIVE_INDEX_FILE_NAME = "refractive_index.txt";

    private double[] wvl;
    private double[] refractiveIndexReal;
    private double[] refractiveIndexImag;
    private int length;
    private String filename;

    public RefractiveIndexTable() {
        length = REFRACTIVE_INDEX_TABLE_LENGTH;
        filename = REFRACTIVE_INDEX_FILE_NAME;
        wvl = new double[length];
        refractiveIndexReal = new double[length];
        refractiveIndexImag = new double[length];
        readTableFromFile();
    }

    public RefractiveIndexTable(int length) {
        this.length = length;
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
                    // x (whatever that is)
                    setWvl(i, Double.parseDouble(st.nextToken()));
                }
                if (st.hasMoreTokens()) {
                    // real part
                    setRefractiveIndexReal(i, Double.parseDouble(st.nextToken()));
                }
                if (st.hasMoreTokens()) {
                    // imaginary part
                    setRefractiveIndexImag(i, Double.parseDouble(st.nextToken()));
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

    public void setRefractiveIndexImag(double[] refractiveIndexImag) {
        this.refractiveIndexImag = refractiveIndexImag;
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

    public void setRefractiveIndexReal(int index, double refractiveIndexReal) {
        this.refractiveIndexReal[index] = refractiveIndexReal;
    }

    public double[] getRefractiveIndexImag() {
        return refractiveIndexImag;
    }

    public double getRefractiveIndexImag(int index) {
        return refractiveIndexImag[index];
    }

    public void setRefractiveIndexImag(int index, double refractiveIndexImag) {
        this.refractiveIndexImag[index] = refractiveIndexImag;
    }
}

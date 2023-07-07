/*
 * Scagnostics
 *
 * Leland Wilkinson and Anushka Anand (University of Illinois at Chicago)
 * This program accompanies the following paper:

 * Wilkinson L., Anand, A., and Grossman, R. (2006). High-Dimensional visual analytics:
 *   Interactive exploration guided by pairwise views of point distributions.
 *   IEEE Transactions on Visualization and Computer Graphics, November/December 2006 (Vol. 12, No. 6) pp. 1363-1372.
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software.
 * Supporting documentation must also include a citation of
 * the abovementioned article.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, THE AUTHORS MAKE NO
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */
package RScag.scagnostics.src.main.java;

import java.io.File;
import java.io.FileFilter;
import java.util.*;

public class Main {

    public static void main(String[] argv) {
        int numBins = 50; // user setting for number of bins
        int maxBins = 1000; // user setting for maximum number of nonempty bins allowed (maxBins >=
                            // numBins*numBins)

        String fold_path = "Rscag//data//";
        File[] files = getFileList(fold_path, ".csv");
        for (File file : files) {
            System.out.println(file.getName());
            double[][] points = getData(file);
            double[][] scagnostics = computeScagnostics(points, numBins, maxBins, file.getName());
            System.out.println(Arrays.toString(scagnostics[0]));
        }
    }

    private static File[] getFileList(String path, String suffix) {
        suffix = suffix.toLowerCase();
        if (!suffix.startsWith(".")) {
            suffix = "." + suffix;
        }
        final String suffixFinal = suffix;
        File directory = new File(path);
        File[] files = directory.listFiles(new FileFilter() {
            public boolean accept(File pathname) {
                if (pathname.getName().toLowerCase().endsWith(suffixFinal))
                    return true;
                return false;
            }
        });
        return files;
    }

    private static double[][] getData(File fname) {
        java.io.BufferedReader fin;
        try {
            fin = new java.io.BufferedReader(new java.io.FileReader(fname));
        } catch (java.io.FileNotFoundException fe) {
            javax.swing.JOptionPane.showMessageDialog(null, "File not found!", "Alert",
                    javax.swing.JOptionPane.ERROR_MESSAGE);
            return null;
        }
        try {
            String record;
            record = fin.readLine();
            // Get the scagnosticsLabels
            record = replaceSeparatorsWithBlanks(record);
            StringTokenizer st = new StringTokenizer(record, " ");
            int col = 0;
            int numVars = st.countTokens();
            String[] variableLabels = new String[numVars];

            while (st.hasMoreTokens()) {
                variableLabels[col] = st.nextToken();
                col++;
            }
            // Count the number of rows
            int numRows = 0;
            // record = fin.readLine();//ignore the label
            while (record != null) {
                record = fin.readLine();
                numRows++;
            }
            fin.close();
            System.out.println("Number of rows, cols " + numRows + " " + numVars);

            // Read in the data
            fin = new java.io.BufferedReader(new java.io.FileReader(fname));
            double[][] data = new double[numVars][numRows];
            // record = fin.readLine(); //ignore line with scagnosticsLabels
            record = fin.readLine();
            int j = 0;
            while (record != null) {
                record = replaceSeparatorsWithBlanks(record);
                st = new StringTokenizer(record, " ");
                for (int i = 0; i < numVars; i++) {
                    try {
                        String tmp = st.nextToken();
                        data[i][j] = Double.parseDouble(tmp);
                    } catch (Exception ie) {
                        javax.swing.JOptionPane.showMessageDialog(null, "Error reading from the file", "Alert",
                                javax.swing.JOptionPane.ERROR_MESSAGE);
                        return null;
                    }
                }
                record = fin.readLine();
                j++;
            }
            fin.close();

            return data;
        } catch (java.io.IOException ie) {
            javax.swing.JOptionPane.showMessageDialog(null, "Error reading from the file", "Alert",
                    javax.swing.JOptionPane.ERROR_MESSAGE);
            return null;
        }
    }

    private static String replaceSeparatorsWithBlanks(String record) {
        record = replaceAll(record, ",,", ",~,");
        record = replaceAll(record, "\t\t", "\t~\t");
        record = replaceAll(record, ",", " ");
        record = replaceAll(record, "\t", " ");
        return record;
    }

    private static String replaceAll(String source, String toReplace, String replacement) {
        int idx = source.lastIndexOf(toReplace);
        if (idx != -1) {
            StringBuffer sb = new StringBuffer(source);
            sb.replace(idx, idx + toReplace.length(), replacement);
            while ((idx = source.lastIndexOf(toReplace, idx - 1)) != -1) {
                sb.replace(idx, idx + toReplace.length(), replacement);
            }
            source = sb.toString();
        }

        return source;
    }

    private static double[][] computeScagnostics(double[][] points, int numBins, int maxBins, String filename) {
        normalizePoints(points);
        int nDim = points.length;
        int numCells = nDim * (nDim - 1) / 2;
        double[][] scagnostics = new double[numCells][Scagnostics.getNumScagnostics()];
        int k = 0;
        for (int i = 1; i < nDim; i++) {
            for (int j = 0; j < i; j++) {
                Scagnostics s = new Scagnostics(points[j], points[i], numBins, maxBins);
                scagnostics[k] = s.compute();
                s.outputResult("scag_results/" + filename.split(".csv")[0] + ".scag.csv", scagnostics[k]);
                k++;
            }
        }
        return scagnostics;
    }

    private static void normalizePoints(double[][] points) {
        double[] min = new double[points.length];
        double[] max = new double[points.length];
        for (int i = 0; i < points.length; i++)
            for (int j = 0; j < points[0].length; j++) {
                if (j == 0)
                    max[i] = min[i] = points[i][0];
                else if (min[i] > points[i][j])
                    min[i] = points[i][j];
                else if (max[i] < points[i][j])
                    max[i] = points[i][j];
            }
        for (int i = 0; i < points.length; i++)
            for (int j = 0; j < points[0].length; j++)
                points[i][j] = (points[i][j] - min[i]) / (max[i] - min[i]);
    }
}

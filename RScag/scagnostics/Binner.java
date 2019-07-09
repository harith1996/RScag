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
package RScag.scagnostics;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Vector;

public class Binner {
    private int maxBins;

    public Binner(int maxBins) {
        this.maxBins = maxBins;
    }

    public final BinnedData binHex(double[] x, double[] y, int nBins, int scale) {

        int n = x.length;

        int mBins = nBins+20;
        int nBin = mBins*mBins;

        double[] count = new double[nBin];
        double[] xbin = new double[nBin];
        double[] ybin = new double[nBin];
        Vector<Double>[] x_inbin = new Vector[nBin];
        Vector<Double>[] y_inbin = new Vector[nBin];

        //ini
        for(int i=0; i<nBin; i++)
        {
            x_inbin[i] = new Vector<Double>();
            y_inbin[i] = new Vector<Double>();
        }
        for (int i=0; i<n; i++)
        {
            if (Double.isNaN(x[i])) continue;
            if (Double.isNaN(y[i])) continue;
            if(x[i]==1){
            if(y[i]==1){
                x_inbin[nBin-1].add(x[i]);
                y_inbin[nBin-1].add(y[i]);
            }else {
                int location = mBins * ((int) Math.floor(y[i] * mBins) + 1) - 1;
                x_inbin[location].add(x[i]);
                y_inbin[location].add(y[i]);
            }
        }
        else if(y[i]==1)
        {
            int location = mBins*((int)Math.floor(x[i]*mBins)+1)-1;
            x_inbin[location].add(x[i]);
            y_inbin[location].add(y[i]);
        }
        else {
            int xlocat = (int) Math.floor(x[i] * mBins);
            int ylocat = (int) Math.floor(y[i] * mBins);
            int location = ylocat * mBins + xlocat;
            x_inbin[location].add(x[i]);
            y_inbin[location].add(y[i]);
        }
    }
    int m = 0;
    for(int i=0 ; i<nBin ; i++)
    {
        if(x_inbin[i].size()>0)
        {
            if(x_inbin[i].size()==1)//keep one node
            {
                xbin[m] = x_inbin[i].get(0);
                ybin[m] = y_inbin[i].get(0);
                count[m] = 1;
                m++;
            }else{
                int scount = x_inbin[i].size()/scale;
                if(scount==0)//keep one node
                {
                    xbin[m] = x_inbin[i].get(randSelect(x_inbin[i].size(),1)[0]);
                    ybin[m] = y_inbin[i].get(randSelect(y_inbin[i].size(),1)[0]);
                    count[m] = x_inbin[i].size();
                    m++;
                }
                else{//keep scount nodes
                    int[] select = randSelect(x_inbin[i].size(), scount);
                    for(int j=0 ; j<scount ; j++)
                    {
                        xbin[m] = x_inbin[i].get(select[j]);
                        ybin[m] = y_inbin[i].get(select[j]);
                        count[m] = (double)(x_inbin[i].size())/scount;
                        m++;
                    }
                }
            }
        }
    }
        nBin = deleteEmptyBins(count, xbin, ybin);

        if (nBin > maxBins) {
            nBins = 2 * nBins / 3;
            scale = (int)(1.5*scale);
            return binHex(x, y, nBins, scale);
        }
        double[] tcount = new double[nBin];
        double[] xtbin = new double[nBin];
        double[] ytbin = new double[nBin];
        System.arraycopy(count, 0, tcount, 0, nBin);//copy count to tcount
        System.arraycopy(xbin, 0, xtbin, 0, nBin);
        System.arraycopy(ybin, 0, ytbin, 0, nBin);
        return new BinnedData(xtbin, ytbin, tcount);
    }
    private int[] randSelect (int n, int m)
    {
        int[] rel = new int[m];
        for(int i = 0; i< m ; i++)
        {
            int tmp = (int)(Math.random()*n);
            while (isInclude(rel ,tmp))
            {
                tmp = (int)(Math.random()*n);
            }
            rel[i] = tmp;
        }
        return rel;
    }
    private boolean isInclude (int[] arr, int target)
    {
        boolean isIn = false;
        for(int i=0 ; i<arr.length ; i++)
        {
            if(arr[i]==target)
            {
                isIn = true;
            }
        }
        return isIn;
    }
    private int deleteEmptyBins(double[] count, double[] xbin, double[] ybin) {

        int k = 0;
        for (int i = 0; i < count.length; i++) {
            if (count[i] > 0) {
                count[k] = count[i];
                xbin[k] = xbin[i];
                ybin[k] = ybin[i];
                k++;
            }
        }
        return k;
    }
}
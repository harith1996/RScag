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

import javax.swing.plaf.basic.BasicScrollPaneUI;
import javax.swing.text.html.HTMLDocument;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.TransferQueue;

public class Scagnostics {
    private BinnedData bdata;
    private List nodes;        // nodes set
    private List edges;        // edges set
    private List triangles;    // triangles set
    private List mstEdges;     // minimum spanning tree set
    private Edge hullStart;   // entering edge of convex hull
    private Edge actE;
    private int totalPeeledCount;
    private int totalPeeledCountBackup;
    private int totalCount;
    private double alphaArea = 1, alphaPerimeter = 1, hullArea = 1, hullPerimeter = 1;
    private double totalOriginalMSTLengths;
    private double totalMSTOutlierLengths;
    private double[] sortedOriginalMSTLengths;
    private static int numScagnostics = 9;
    private final static int OUTLYING = 0, SKEWED = 1, CLUMPY = 2, SPARSE = 3,
            STRIATED = 4, CONVEX = 5, SKINNY = 6, STRINGY = 7, MONOTONIC = 8;
    private final static String[] scagnosticsLabels = {"Outlying", "Skewed", "Clumpy", "Sparse",
            "Striated", "Convex", "Skinny", "Stringy", "Monotonic"};
    private double[] counts;
    private int[] px, py;
    private boolean[] isOutlier;
    private boolean[] isOutlierbackup;
    private List<Edge> outlierEdges;
    private List<Edge> runtsEdges;
    private List<Edge> maxEdges;
    private double rel_outlying = 0;
    private double rel_clumpy = 0;
    private List<Double> maxValues;

    private Vector<Node>[] subCluster;
    private double FUZZ = .999;

    public Scagnostics(double[] x, double[] y, int numBins, int maxBins) {
        nodes = new ArrayList();
        edges = new ArrayList();
        triangles = new ArrayList();
        mstEdges = new ArrayList();
        Binner b = new Binner(maxBins);
        bdata = b.binHex(x, y, numBins, 3);
    }

    public double[] compute() {
        px = bdata.getXData();
        py = bdata.getYData();
        if (px.length < 3)
            return null;
        int xx = px[0];
        int yy = py[0];
        boolean isXConstant = true;
        boolean isYConstant = true;
        for (int i = 1; i < px.length; i++) {
            if (px[i] != xx) isXConstant = false;
            if (py[i] != yy) isYConstant = false;
        }
        if (isXConstant || isYConstant)
            return null;


        findOutliers(bdata);
        double[] result = new double[numScagnostics];
        result[OUTLYING] = rel_outlying;
        result[CLUMPY] = rel_clumpy;
        for (int i = 0; i < px.length; i++) {
            if (subCluster[i].size() > 3) {
                allOutlierExcept(subCluster[i]);
                double weight = (double) subCluster[i].size() / px.length;
                clear();
                computeDT(px, py);
                computeMST();
                computeAlphaGraph();
                computeTotalCount();
                computeAlphaArea();
                computeAlphaPerimeter();
                computeHullArea();
                computeHullPerimeter();
                double[] tmpResult = computeMeasures();
                result[SKEWED] += tmpResult[SKEWED] * weight;
                result[CONVEX] = tmpResult[CONVEX] * weight;
                result[SKINNY] = tmpResult[SKINNY] * weight;
                result[STRINGY] = tmpResult[STRINGY] * weight;
                result[STRIATED] = tmpResult[STRIATED] * weight;
                result[SPARSE] = tmpResult[SPARSE] * weight;
                result[MONOTONIC] = tmpResult[MONOTONIC] * weight;
            }
        }
        return result;
    }

    public static int getNumScagnostics() {
        return scagnosticsLabels.length;
    }

    public static String[] getScagnosticsLabels() {
        return scagnosticsLabels;
    }

    public static boolean[] computeScagnosticsExemplars(double[][] pts) {
        int nPts = pts.length;
        if (nPts < 2)
            return null;
        Cluster c = new Cluster(0, 0);
        int[] exemp = c.compute(pts);
        boolean[] exemplars = new boolean[nPts];
        for (int i = 0; i < exemp.length; i++)
            exemplars[exemp[i]] = true;
        return exemplars;
    }

    public static boolean[] computeScagnosticsOutliers(double[][] pts) {

        // Prim's algorithm

        int nPts = pts.length;     // p*(p-1)/2 points representing pairwise scatterplots
        int nVar = pts[0].length;  // number of scagnostics (9)
        if (nPts < 2)
            return null;
        int[][] edges = new int[nPts - 1][2];
        int[] list = new int[nPts];
        int[] degrees = new int[nPts];
        double[] cost = new double[nPts];
        double[] lengths = new double[nPts - 1];

        list[0] = 0;
        cost[0] = Double.POSITIVE_INFINITY;
        int cheapest = 0;

        for (int i = 1; i < nPts; i++) {
            for (int j = 0; j < nVar; j++) {
                double d = pts[i][j] - pts[0][j];
                cost[i] += d * d;
            }
            if (cost[i] < cost[cheapest])
                cheapest = i;
        }
        for (int j = 1; j < nPts; j++) {
            int end = list[cheapest];
            int jp = j - 1;
            edges[jp][0] = cheapest;
            edges[jp][1] = end;
            lengths[jp] = cost[cheapest];
            degrees[cheapest]++;
            degrees[end]++;
            cost[cheapest] = Double.POSITIVE_INFINITY;
            end = cheapest;

            for (int i = 1; i < nPts; i++) {
                if (cost[i] != Double.POSITIVE_INFINITY) {
                    double dist = 0.;
                    for (int k = 0; k < nVar; k++) {
                        double d = pts[i][k] - pts[end][k];
                        dist += d * d;
                    }
                    if (dist < cost[i]) {
                        list[i] = end;
                        cost[i] = dist;
                    }
                    if (cost[i] < cost[cheapest]) cheapest = i;
                }
            }
        }
        double cutoff = findCutoff(lengths);
        boolean[] outliers = new boolean[nPts];
        for (int i = 0; i < nPts; i++)
            outliers[i] = true;
        for (int i = 0; i < nPts - 1; i++) {
            if (lengths[i] < cutoff) {
                for (int k = 0; k < 2; k++) {
                    int node = edges[i][k];
                    outliers[node] = false;
                }
            }
        }
        return outliers;
    }

    private void clear() {
        nodes.clear();
        edges.clear();
        triangles.clear();
        mstEdges.clear();
    }

    private void findOutliers(BinnedData bdata) {
        this.counts = bdata.getCounts();
        isOutlier = new boolean[px.length];
        isOutlierbackup = new boolean[px.length];
        outlierEdges = new ArrayList<>();
        runtsEdges = new ArrayList<>();
        maxEdges = new ArrayList<>();
        maxValues = new ArrayList<Double>();

        subCluster = new Vector[px.length];
        for (int i = 0; i < px.length; i++) {
            subCluster[i] = new Vector<Node>();
        }

        computeDT(px, py);
        computeMST();
        sortedOriginalMSTLengths = getSortedMSTEdgeLengths();
        computeTotalOriginalMSTLengths();
        for (int i = 0; i < nodes.size(); i++) {
            Node cur_nd = (Node) nodes.get(i);
            subCluster[0].add(cur_nd);
        }
        //ini
        departNodes(0);
        //re
        isOutlier = deepClone(isOutlierbackup);

        clear();
        computeDT(px, py);
        computeMST();

        sortedOriginalMSTLengths = getSortedMSTEdgeLengths();

        if (subCluster[1].size() == 0) {
            System.out.println("single cluster.");
            getRuntsAndMaxEdgeInOneClu();
            rel_clumpy = maxValues.get(0);
        } else {
            int runt_sz = runtsEdges.size();
            int sum = 0;
            Vector<Integer> Counts = new Vector<Integer>();
            for (int i = 0; i < runt_sz; i++) {
                runtsEdges.get(i).onMST = false;
            }
            for (int i = 0; i < runt_sz; i++) {
                Node p1 = runtsEdges.get(i).p1;
                Node p2 = runtsEdges.get(i).p2;
                int index1 = findClosestClu(p1);
                int index2 = findClosestClu(p2);
                int count1 = subCluster[index1].size();
                int count2 = subCluster[index2].size();
                Counts.add((count1 + count2));
                sum += (count1 + count2);
            }
            for (int i = 0; i < runt_sz; i++) {
                rel_clumpy += (maxValues.get(i) * Counts.get(i)) / sum;
            }
        }
    }

    private int findClosestClu(Node p) {
        int index = 0;
        //find the nearest node
        int nod_sz = nodes.size();
        Node nearest_nd = (Node) nodes.get(0);
        double nearDis = Double.MAX_VALUE;
        for (int i = 0; i < nod_sz; i++) {
            Node cur_nd = (Node) nodes.get(i);
            double distance = cur_nd.distToNode(p.x, p.y);
            if (distance < nearDis) {
                nearest_nd = cur_nd;
                nearDis = distance;
            }
        }
        int clu_sz = subCluster.length;
        for (int i = 0; i < clu_sz; i++) {
            int cur_clu_sz = subCluster[i].size();
            for (int j = 0; j < cur_clu_sz; j++) {
                Node cur_nd = subCluster[i].get(j);
                if (cur_nd.distToNode(nearest_nd.x, nearest_nd.y) == 0) {
                    index = i;
                    break;
                }
            }
        }
        return index;
    }

    private boolean[] deepClone(boolean[] target) {
        boolean[] rel = new boolean[px.length];
        for (int i = 0; i < px.length; i++) {
            rel[i] = target[i];
        }
        return rel;
    }

    private void departNodes(int index) {
        if (checkSingleClu(index)) {
            return;
        }
        allOutlierExcept(subCluster[index]);
        subCluster[index].clear();
        clear();
        computeDT(px, py);
        computeMST();

        double[] sortedmst = getSortedMSTEdgeLengths();
        double cutoff = computeCutoff(sortedmst);
        boolean foundOutliers = computeMSTOutliers(cutoff);
        double[] sortedPeeledMSTLengths;
        while (foundOutliers) {
            clear();
            computeDT(px, py);
            computeMST();
            sortedPeeledMSTLengths = getSortedMSTEdgeLengths();
            cutoff = computeCutoff(sortedPeeledMSTLengths);
            foundOutliers = computeMSTOutliers(cutoff);
        }
        computeClusterMeasure(cutoff);

    }

    public void outputResult(String filename, double[] scagnosticsRel) {
        String strtemp;
        File file = new File(filename);
        for (int i = 0; i < numScagnostics; i++) {
            strtemp = scagnosticsLabels[i] + "," + scagnosticsRel[i];
            writeToFile(file, filename, strtemp);
        }
    }

    public void outputNode(String filename, int[] x, int[] y) {
        String strtemp;
        File file = new File(filename);
        for (int i = 0; i < x.length; i++) {
            strtemp = x[i] + "," + y[i];
            writeToFile(file, filename, strtemp);
        }
    }

    public void outputMSTedge(String filename, List<Edge> targetEdges) {
        Node p1, p2;
        Edge temp;
        String strtemp;
        File file = new File(filename);
        for (int i = 0; i < targetEdges.size(); i++) {
            temp = (Edge) targetEdges.get(i);
            p1 = temp.p1;
            p2 = temp.p2;
            strtemp = p1.x + "," + p1.y + "," + p2.x + "," + p2.y;
            writeToFile(file, filename, strtemp);
        }
    }

    public void writeToFile(File file, String filename, String content) {
        try {
            if (!file.exists()) {
                file.createNewFile();
                System.out.println("Create file successful!");
                writeFileContent(filename, content);
            } else {
                writeFileContent(filename, content);
            }
        } catch (Exception exc) {
            exc.printStackTrace();
        }
    }

    public static boolean writeFileContent(String filepath, String newstr)
            throws IOException {
        Boolean bool = false;
        String filein = newstr + "\r\n";
        String temp = "";

        FileInputStream fis = null;
        InputStreamReader isr = null;
        BufferedReader br = null;
        FileOutputStream fos = null;
        PrintWriter pw = null;
        try {
            File file = new File(filepath);
            fis = new FileInputStream(file);
            isr = new InputStreamReader(fis);
            br = new BufferedReader(isr);
            StringBuffer buffer = new StringBuffer();


            for (int i = 0; (temp = br.readLine()) != null; i++) {
                buffer.append(temp);

                buffer = buffer.append(System.getProperty("line.separator"));
            }
            buffer.append(filein);

            fos = new FileOutputStream(file);
            pw = new PrintWriter(fos);
            pw.write(buffer.toString().toCharArray());
            pw.flush();
            bool = true;
        } catch (Exception e) {
            // TODO: handle exception
            e.printStackTrace();
        } finally {

            if (pw != null) {
                pw.close();
            }
            if (fos != null) {
                fos.close();
            }
            if (br != null) {
                br.close();
            }
            if (isr != null) {
                isr.close();
            }
            if (fis != null) {
                fis.close();
            }
        }
        return bool;
    }

    private boolean checkSingleClu(int index) {
        boolean isSingleClu = true;
        allOutlierExcept(subCluster[index]);
        clear();
        computeDT(px, py);
        computeMST();
        double[] sortedPeeledMSTLengths = getSortedMSTEdgeLengths();
        double cutoff = computeCutoff(sortedPeeledMSTLengths);
        int tmp_size = sortedPeeledMSTLengths.length;
        for (int j = 0; j < tmp_size; j++) {
            if (sortedPeeledMSTLengths[j] > cutoff) {
                isSingleClu = false;
            }
        }
        return isSingleClu;
    }

    private void allOutlier() {
        for (int i = 0; i < px.length; i++) {
            isOutlier[i] = true;
        }
    }

    private void allOutlierExcept(Vector<Node> nods) {
        int tmp_size = nods.size();
        allOutlier();
        for (int j = 0; j < tmp_size; j++) {
            isOutlier[nods.get(j).pointID] = false;
        }
    }

    private int addNodesToSubClu(Vector<Node> nods) {
        int tmp_size = nods.size();
        int index = 0;
        for (int i = 0; i < px.length; i++) {
            if (subCluster[i].size() == 0) {
                index = i;
                for (int j = 0; j < tmp_size; j++) {

                    if (!isOutlierbackup[nods.get(j).pointID]) {
                        subCluster[i].add(nods.get(j));
                    }
                }
                break;
            }
        }
        return index;
    }

    private boolean isInSubClu(Node n) {
        boolean isIn = false;
        for (int i = 0; i < px.length; i++) {
            int tmp_size = subCluster[i].size();
            if (tmp_size != 0) {
                for (int j = 0; j < tmp_size; j++) {
                    if (subCluster[i].get(j).pointID == n.pointID) {
                        isIn = true;
                    }
                }
            }
        }
        return isIn;
    }


    private void computeTotalCount() {
        for (int i = 0; i < counts.length; i++) {
            totalCount += counts[i];
        }
    }

    // compute measures except outlying and clumpy
    private double[] computeMeasures() {
        double[] results = new double[numScagnostics];
        // Do not change order of these calls!
        results[SKEWED] = computeMSTEdgeLengthSkewnessMeasure();
        results[CONVEX] = computeConvexityMeasure();
        results[SKINNY] = computeSkinnyMeasure();
        results[STRINGY] = computeStringyMeasure();
        results[STRIATED] = computeStriationMeasure();
        results[SPARSE] = computeSparsenessMeasure();
        results[MONOTONIC] = computeMonotonicityMeasure();
        return results;
    }

    private void computeDT(int[] px, int[] py) {
        totalPeeledCount = 0;
        Random r = new Random(13579);

        for (int i = 0; i < px.length; i++) {
            int x = px[i] + (int) (8 * (r.nextDouble() - .5)); // perturb to prevent singularities
            int y = py[i] + (int) (8 * (r.nextDouble() - .5));
            double count = counts[i];
            if (!isOutlier[i]) {
                insert(x, y, count, i);
                totalPeeledCount += count;
            }
        }
        setNeighbors();
        markHull();
    }


    private void computeMST() {
        if (nodes.size() > 1) {
            List mstNodes = new ArrayList();
            Node mstNode = (Node) nodes.get(0);
            updateMSTNodes(mstNode, mstNodes);
            int count = 1;
            while (count < nodes.size()) {
                Edge addEdge = null;
                double wmin = Double.MAX_VALUE;
                Node nmin = null;
                Iterator mstIterator = mstNodes.iterator();
                while (mstIterator.hasNext()) {
                    mstNode = (Node) mstIterator.next();
                    Edge candidateEdge = mstNode.shortestEdge(false);
                    if (candidateEdge != null) {
                        double wt = candidateEdge.weight;
                        if (wt < wmin) {
                            wmin = wt;
                            nmin = mstNode;
                            addEdge = candidateEdge;
                        }
                    }
                }
                if (addEdge != null) {
                    Node addNode = addEdge.otherNode(nmin);
                    updateMSTNodes(addNode, mstNodes);
                    updateMSTEdges(addEdge, mstEdges);
                }
                count++;
            }
        }
    }

    private static double findCutoff(double[] distances) {
        int[] index = Sorts.indexedDoubleArraySort(distances, 0, 0);
        int n50 = distances.length / 2;
        int n25 = n50 / 2;
        int n75 = n50 + n50 / 2;
        return distances[index[n75]] + 1.5 * (distances[index[n75]] - distances[index[n25]]);
    }

    private boolean computeMSTOutliers(double omega) {
        boolean found = false;
        Iterator it = nodes.iterator();
        while (it.hasNext()) {
            Node n = (Node) it.next();
            if (n.getMstDegree() != 1) {
                continue;
            }
            Iterator ie = n.neighbors.iterator();
            while (ie.hasNext()) {
                Edge e = (Edge) ie.next();
                if (e.onMST && e.weight > omega && !e.onOutlier) {
                    found = true;
                    totalMSTOutlierLengths += e.weight;
                    e.onOutlier = true;
                    isOutlier[n.pointID] = true;
                    isOutlierbackup[n.pointID] = true;
                    e.onMST = false;
                    outlierEdges.add(e);
                    rel_outlying += e.weight / totalOriginalMSTLengths;
                }
            }
        }
        return found;
    }

    private double computeCutoff(double[] lengths) {
        if (lengths.length == 0) return 0;
        int n50 = lengths.length / 2;
        int n25 = n50 / 2;
        int n75 = n50 + n25;
        return lengths[n75] + 1.5 * (lengths[n75] - lengths[n25]);
    }

    private double computeAlphaValue() {
        int length = sortedOriginalMSTLengths.length;
        if (length == 0) return 100.;
        int n90 = (9 * length) / 10;
        double alpha = sortedOriginalMSTLengths[n90];
        return Math.min(alpha, 100.);
    }


    private double getMeanValue(int start, int end, double[] array) {
        double mean;
        double sum = 0;
        for (int i = start; i < end; i++) {
            sum += array[i];
        }
        mean = sum / (end - start);
        return mean;
    }

    private double computeMSTEdgeLengthSkewnessMeasure() {
        if (sortedOriginalMSTLengths.length == 0)
            return 0;
        int n = sortedOriginalMSTLengths.length;
        int n80 = n * 4 / 5;
        int n20 = n / 5;
        int n50 = n / 2;
        double Mean80to100 = getMeanValue(n80, n, sortedOriginalMSTLengths);
        double Mean0to20 = getMeanValue(0, n20, sortedOriginalMSTLengths);
        double skewness = (Mean80to100 - sortedOriginalMSTLengths[n50]) / (Mean80to100 - Mean0to20);
        double t = (double) totalCount / 500;
        double correction = .7 + .3 / (1 + t * t);
        return 1 - correction * (1 - skewness);
    }

    private void updateMSTEdges(Edge addEdge, List mstEdges) {
        mstEdges.add(addEdge);
        addEdge.onMST = true;
        addEdge.p1.mstDegree++;
        addEdge.p2.mstDegree++;
    }

    private void updateMSTNodes(Node addNode, List mstNodes) {
        mstNodes.add(addNode);
        addNode.onMST = true;
    }

    private double[] getSortedMSTEdgeLengths() {
        double[] lengths = computeEdgeLengths(mstEdges.iterator(), mstEdges.size());
        Sorts.doubleArraySort(lengths, 0, 0);
        return lengths;
    }

    private void computeTotalOriginalMSTLengths() {
        for (int i = 0; i < sortedOriginalMSTLengths.length; i++)
            totalOriginalMSTLengths += sortedOriginalMSTLengths[i];
    }

    private double[] computeEdgeLengths(Iterator graph, int n) {
        double[] lengths = new double[n];
        int i = 0;
        while (graph.hasNext()) {
            Edge e = (Edge) graph.next();
            lengths[i] = e.weight;
            i++;
        }
        return lengths;
    }

    private boolean pointsInCircle(Node n, double xc, double yc, double radius) {
        double r = FUZZ * radius;
        Iterator i = n.neighbors.iterator();
        while (i.hasNext()) {
            Edge e = (Edge) i.next();
            Node no = e.otherNode(n);
            double dist = no.distToNode(xc, yc);
            if (dist < r)
                return true;
        }
        return false;
    }

    private void computeAlphaGraph() { // requires initializing SEdge.onShape = false
        boolean deleted;
        double alpha = computeAlphaValue();
        do {
            Iterator i = edges.iterator();
            deleted = false;
            while (i.hasNext()) {
                Edge e = (Edge) i.next();
                if (e.inT.onComplex) {
                    if (alpha < e.weight / 2) {
                        e.inT.onComplex = false;
                        deleted = true;
                    } else {
                        if (e.invE != null)
                            if (e.invE.inT.onComplex)
                                continue;
                        if (!edgeIsExposed(alpha, e)) {
                            e.inT.onComplex = false;
                            deleted = true;
                        }
                    }
                }
            }
        } while (deleted);
        markShape();
    }

    private void markShape() {
        Iterator i = edges.iterator();
        while (i.hasNext()) {
            Edge e = (Edge) i.next();
            e.onShape = false;
            if (e.inT.onComplex) {
                if (e.invE == null) {
                    e.onShape = true;
                } else if (!e.invE.inT.onComplex)
                    e.onShape = true;
            }
        }
    }

    private boolean edgeIsExposed(double alpha, Edge e) {
        double x1 = e.p1.x;
        double x2 = e.p2.x;
        double y1 = e.p1.y;
        double y2 = e.p2.y;
        double xe = (x1 + x2) / 2;
        double ye = (y1 + y2) / 2;
        double d = Math.sqrt(alpha * alpha - e.weight * e.weight / 4);
        double xt = d * (y2 - y1) / e.weight;
        double yt = d * (x2 - x1) / e.weight;
        double xc1 = xe + xt;
        double yc1 = ye - yt;
        double xc2 = xe - xt;
        double yc2 = ye + yt;
        boolean pointsInCircle1 = pointsInCircle(e.p1, xc1, yc1, alpha) ||
                pointsInCircle(e.p2, xc1, yc1, alpha);
        boolean pointsInCircle2 = pointsInCircle(e.p1, xc2, yc2, alpha) ||
                pointsInCircle(e.p2, xc2, yc2, alpha);
        return !(pointsInCircle1 && pointsInCircle2);
    }

    private double computeStringyMeasure() {
        int count1 = 0;
        int count2 = 0;
        Iterator it = nodes.iterator();
        while (it.hasNext()) {
            Node n = (Node) it.next();
            if (n.mstDegree == 1)
                count1++;
            if (n.mstDegree == 2)
                count2++;
        }
        double result = (double) count2 / (double) (nodes.size() - count1);
        return result * result * result;
    }

    private void computeClusterMeasure(double cutoff) {

        Iterator it0 = mstEdges.iterator();
        boolean allSmaller = true;
        while (it0.hasNext()) {
            Edge cur_edge = (Edge) it0.next();
            if (cur_edge.weight > 1.6 * cutoff) {
                allSmaller = false;
                break;
            }
        }
        if (allSmaller) {

            int nod_sz = nodes.size();
            Vector<Node> nods = new Vector<>();
            for (int i = 0; i < nod_sz; i++) {
                nods.add((Node) nodes.get(i));
            }
            if (!isInSubClu(nods.get(0))) {
                addNodesToSubClu(nods);
            }

            return;
        }


        Edge runt_edge = getRuntsAndMaxEdge(cutoff);
        Node p1 = runt_edge.p1;
        Node p2 = runt_edge.p2;
        Vector<Node> node_p1 = new Vector<>();
        Vector<Node> node_p2 = new Vector<>();

        runt_edge.onMST = false;
        clearVisits();
        int count1 = p1.getMSTCount(node_p1);
        int count2 = p2.getMSTCount(node_p2);

        if (!isInSubClu(p1)) {
            int p1_index = addNodesToSubClu(node_p1);
            departNodes(p1_index);
        }
        if (!isInSubClu(p2)) {
            int p2_index = addNodesToSubClu(node_p2);
            departNodes(p2_index);
        }
    }

    private Edge getRuntsAndMaxEdge(double cutoff) {
        Iterator it = mstEdges.iterator();
        double[] maxLength = new double[1];
        double maxValue = 0;

        Edge[] maxEdge = new Edge[2];
        Node p1_tmp = new Node(0, 0, 1, 0);
        Node p2_tmp = new Node(0, 0, 1, 1);
        Edge tmpEdge1 = new Edge(p1_tmp, p2_tmp);
        Edge tmpEdge2 = new Edge(p1_tmp, p2_tmp);
        Edge runt_edge = new Edge(p1_tmp, p2_tmp);
        double tmp_value = 0;
        while (it.hasNext()) {
            Edge e = (Edge) it.next();
            clearVisits();
            e.onMST = false;  // break MST at this edge
            int runts = e.getRunts(maxLength, maxEdge);
            e.onMST = true;   // restore this edge to MST
            if (e.weight > 1.6 * cutoff && maxLength[0] > 0) {
                double value = runts * (1 - maxLength[0] / e.weight);
//                double value = 1 - maxLength[0] / e.weight;
                if (value > maxValue) {
                    maxValue = value;
                    tmp_value = (1 - maxLength[0] / e.weight);
                    runt_edge = e;
                    tmpEdge1 = maxEdge[0];//edge in smaller cluster
                    tmpEdge2 = maxEdge[1];//edge in bigger cluster
                }
            }
        }
        runtsEdges.add(runt_edge);
        maxEdges.add(tmpEdge1);
        maxEdges.add(tmpEdge2);
        maxValues.add(tmp_value);
        return runt_edge;
    }

    private void getRuntsAndMaxEdgeInOneClu() {
        Iterator it = mstEdges.iterator();
        double[] maxLength = new double[1];
        double maxValue = 0;

        Edge[] maxEdge = new Edge[2];
        Node p1_tmp = new Node(0, 0, 1, 0);
        Node p2_tmp = new Node(0, 0, 1, 1);
        Edge tmpEdge1 = new Edge(p1_tmp, p2_tmp);
        Edge tmpEdge2 = new Edge(p1_tmp, p2_tmp);
        Edge runt_edge = new Edge(p1_tmp, p2_tmp);
        double tmp_value = 0;
        while (it.hasNext()) {
            Edge e = (Edge) it.next();
            clearVisits();
            e.onMST = false;  // break MST at this edge
            int runts = e.getRunts(maxLength, maxEdge);
            e.onMST = true;   // restore this edge to MST
            if (maxLength[0] > 0 && runts > 1) {
                double value = runts * (1 - maxLength[0] / e.weight);
                if (value > maxValue) {
                    maxValue = value;
                    tmp_value = (1 - maxLength[0] / e.weight);
                    runt_edge = e;
                    tmpEdge1 = maxEdge[0];//edge in smaller cluster
                    tmpEdge2 = maxEdge[1];//edge in bigger cluster
                }
            }
        }
        runtsEdges.add(runt_edge);
        maxEdges.add(tmpEdge1);
        maxEdges.add(tmpEdge2);
        maxValues.add(tmp_value);
    }

    private void clearVisits() {
        Iterator it = nodes.iterator();
        while (it.hasNext()) {
            Node n = (Node) it.next();
            n.isVisited = false;
        }
    }

    private double computeMonotonicityMeasure() {
        int n = counts.length;
        double[] ax = new double[n];
        double[] ay = new double[n];
        double[] weights = new double[n];
        for (int i = 0; i < n; i++) {
            ax[i] = px[i];
            ay[i] = py[i];
            weights[i] = counts[i];
        }
        double[] rx = Sorts.rank(ax);
        double[] ry = Sorts.rank(ay);
        double s = computePearson(rx, ry, weights);
        return s * s;
    }

    private double computePearson(double[] x, double[] y, double[] weights) {
        int n = x.length;
        double xmean = 0;
        double ymean = 0;
        double xx = 0;
        double yy = 0;
        double xy = 0;
        double sumwt = 0;
        for (int i = 0; i < n; i++) {
            double wt = weights[i];
            if (wt > 0 && !isOutlier[i]) {
                sumwt += wt;
                xx += (x[i] - xmean) * wt * (x[i] - xmean);
                yy += (y[i] - ymean) * wt * (y[i] - ymean);
                xy += (x[i] - xmean) * wt * (y[i] - ymean);
                xmean += (x[i] - xmean) * wt / sumwt;
                ymean += (y[i] - ymean) * wt / sumwt;
            }
        }
        xy = xy / Math.sqrt(xx * yy);
        return xy;
    }

    private double computeSparsenessMeasure() {
        int n = sortedOriginalMSTLengths.length;
        int n90 = (9 * n) / 10;
        double sparse = Math.min(sortedOriginalMSTLengths[n90] / 1000, 1);
        double t = (double) totalCount / 500;
        double correction = .7 + .3 / (1 + t * t);
        return correction * sparse;
    }

    private double computeStriationMeasure() {
        double numEdges = 0;
        Iterator it = mstEdges.iterator();
        while (it.hasNext()) {
            Edge e = (Edge) it.next();
            Node n1 = e.p1;
            Node n2 = e.p2;
            if (n1.mstDegree == 2 && n2.mstDegree == 2) {
                Edge e1 = getAdjacentMSTEdge(n1, e);
                Edge e2 = getAdjacentMSTEdge(n2, e);
                if (cosineOfAdjacentEdges(e, e1, n1) < -.7 && cosineOfAdjacentEdges(e, e2, n2) < -.7)
                    numEdges++;
            }
        }
        return numEdges / (double) mstEdges.size();
    }

    private Edge getAdjacentMSTEdge(Node n, Edge e) {
        Iterator nt = n.neighbors.iterator();
        while (nt.hasNext()) {
            Edge et = (Edge) nt.next();
            if (et.onMST && !e.equals(et)) {
                return et;
            }
        }
        return null;
    }

    private double cosineOfAdjacentEdges(Edge e1, Edge e2, Node n) {
        double v1x = e1.otherNode(n).x - n.x;
        double v1y = e1.otherNode(n).y - n.y;
        double v2x = e2.otherNode(n).x - n.x;
        double v2y = e2.otherNode(n).y - n.y;
        double v1 = Math.sqrt(v1x * v1x + v1y * v1y);
        double v2 = Math.sqrt(v2x * v2x + v2y * v2y);
        v1x = v1x / v1;
        v1y = v1y / v1;
        v2x = v2x / v2;
        v2y = v2y / v2;
        return v1x * v2x + v1y * v2y;
    }

    private double computeConvexityMeasure() {
        if (hullArea == 0) // points in general position
            return 1;
        else {
            double t = (double) totalCount / 500;
            double correction = .7 + .3 / (1 + t * t);
            double convexity = alphaArea / hullArea;
            return correction * convexity;
        }
    }

    private double computeSkinnyMeasure() {
        if (alphaPerimeter > 0)
            return 1 - Math.sqrt(4 * Math.PI * alphaArea) / alphaPerimeter;
        else
            return 1;
    }

    private void computeAlphaArea() {
        double area = 0;
        Iterator tri = triangles.iterator();
        while (tri.hasNext()) {
            Triangle t = (Triangle) tri.next();
            if (t.onComplex) {
                Node p1 = t.anEdge.p1;
                Node p2 = t.anEdge.p2;
                Node p3 = t.anEdge.nextE.p2;
                area += Math.abs(p1.x * p2.y + p1.y * p3.x + p2.x * p3.y
                        - p3.x * p2.y - p3.y * p1.x - p1.y * p2.x);
            }
        }
        alphaArea = area / 2;
    }

    private void computeHullArea() {
        double area = 0.0;
        Iterator tri = triangles.iterator();
        while (tri.hasNext()) {
            Triangle t = (Triangle) tri.next();
            Node p1 = t.anEdge.p1;
            Node p2 = t.anEdge.p2;
            Node p3 = t.anEdge.nextE.p2;
            area += Math.abs(p1.x * p2.y + p1.y * p3.x + p2.x * p3.y
                    - p3.x * p2.y - p3.y * p1.x - p1.y * p2.x);
        }
        hullArea = area / 2.;
    }

    private void computeAlphaPerimeter() {
        double sum = 0;
        Iterator it = edges.iterator();
        while (it.hasNext()) {
            Edge e = (Edge) it.next();
            if (e.onShape) {
                sum += e.weight;
            }
        }
        alphaPerimeter = sum;
    }

    private void computeHullPerimeter() {
        double sum = 0;
        Edge e = hullStart;
        do {
            sum += e.p1.distToNode(e.p2.x, e.p2.y);
            e = e.nextH;
        } while (!e.isEqual(hullStart));
        hullPerimeter = sum;
    }

    private void setNeighbors() {
        Iterator it = edges.iterator();
        while (it.hasNext()) {
            Edge e = (Edge) it.next();
            if (e.isNewEdge(e.p1))
                e.p1.setNeighbor(e);
            if (e.isNewEdge(e.p2))
                e.p2.setNeighbor(e);
        }
    }

    private void insert(int px, int py, double count, int id) {
        int eid;
        Node nd = new Node(px, py, count, id);
        nodes.add(nd);
        if (nodes.size() < 3) return;
        if (nodes.size() == 3)    // create the first triangle
        {
            Node p1 = (Node) nodes.get(0);
            Node p2 = (Node) nodes.get(1);
            Node p3 = (Node) nodes.get(2);
            Edge e1 = new Edge(p1, p2);
            if (e1.onSide(p3) == 0) {
                nodes.remove(nd);
                return;
            }
            if (e1.onSide(p3) == -1)  // right side
            {
                p1 = (Node) nodes.get(1);
                p2 = (Node) nodes.get(0);
                e1.update(p1, p2);
            }
            Edge e2 = new Edge(p2, p3);
            Edge e3 = new Edge(p3, p1);
            e1.nextH = e2;
            e2.nextH = e3;
            e3.nextH = e1;
            hullStart = e1;
            triangles.add(new Triangle(edges, e1, e2, e3));
            return;
        }
        //when the size of nodes is bigger than 3
        actE = (Edge) edges.get(0);
        if (actE.onSide(nd) == -1) {
            if (actE.invE == null)
                eid = -1;
            else
                eid = searchEdge(actE.invE, nd);
        } else
            eid = searchEdge(actE, nd);
        if (eid == 0) {
            nodes.remove(nd);
            return;
        }
        if (eid > 0)
            expandTri(actE, nd, eid);   // nd is inside or on a triangle
        else
            expandHull(nd);                // nd is outside convex hull
    }

    private void expandTri(Edge e, Node nd, int type) {
        Edge e1 = e;
        Edge e2 = e1.nextE;
        Edge e3 = e2.nextE;
        Node p1 = e1.p1;
        Node p2 = e2.p1;
        Node p3 = e3.p1;
        if (type == 2) {   // nd is inside of the triangle
            Edge e10 = new Edge(p1, nd);
            Edge e20 = new Edge(p2, nd);
            Edge e30 = new Edge(p3, nd);
            e.inT.removeEdges(edges);
            triangles.remove(e.inT);     // remove old triangle
            Edge e100 = e10.makeSymm();
            Edge e200 = e20.makeSymm();
            Edge e300 = e30.makeSymm();
            triangles.add(new Triangle(edges, e1, e20, e100));
            triangles.add(new Triangle(edges, e2, e30, e200));
            triangles.add(new Triangle(edges, e3, e10, e300));
            swapTest(e1);   // swap test for the three new triangles
            swapTest(e2);
            swapTest(e3);
        } else {          // nd is on the edge e
            Edge e4 = e1.invE;
            if (e4 == null || e4.inT == null) {          // one triangle involved
                Edge e30 = new Edge(p3, nd);
                Edge e02 = new Edge(nd, p2);
                Edge e10 = new Edge(p1, nd);
                Edge e03 = e30.makeSymm();
//								shareEdges(e03,e30);
                e10.asIndex();
                e1.mostLeft().nextH = e10;
                e10.nextH = e02;
                e02.nextH = e1.nextH;
                hullStart = e02;
                triangles.remove(e1.inT);  // remove oldtriangle and add two new triangles
                edges.remove(e1);
                edges.add(e10);
                edges.add(e02);
                edges.add(e30);
                edges.add(e03);
                triangles.add(new Triangle(e2, e30, e02));
                triangles.add(new Triangle(e3, e10, e03));
                swapTest(e2);   // swap test for the two new triangles
                swapTest(e3);
                swapTest(e30);
            } else {        // two triangle involved
                Edge e5 = e4.nextE;
                Edge e6 = e5.nextE;
                Node p4 = e6.p1;
                Edge e10 = new Edge(p1, nd);
                Edge e20 = new Edge(p2, nd);
                Edge e30 = new Edge(p3, nd);
                Edge e40 = new Edge(p4, nd);
                triangles.remove(e.inT);                   // remove oldtriangle
                e.inT.removeEdges(edges);
                triangles.remove(e4.inT);               // remove old triangle
                e4.inT.removeEdges(edges);
                e5.asIndex();   // because e, e4 removed, reset edge sortOrder of node p1 and p2
                e2.asIndex();
                triangles.add(new Triangle(edges, e2, e30, e20.makeSymm()));
                triangles.add(new Triangle(edges, e3, e10, e30.makeSymm()));
                triangles.add(new Triangle(edges, e5, e40, e10.makeSymm()));
                triangles.add(new Triangle(edges, e6, e20, e40.makeSymm()));
                swapTest(e2);   // swap test for the three new triangles
                swapTest(e3);
                swapTest(e5);
                swapTest(e6);
                swapTest(e10);
                swapTest(e20);
                swapTest(e30);
                swapTest(e40);
            }
        }
    }

    private void expandHull(Node nd) {
        Edge e1, e2, e3 = null, enext;
        Edge e = hullStart;
        Edge comedge = null, lastbe = null;
        while (true) {
            enext = e.nextH;
            if (e.onSide(nd) == -1) {  // right side
                if (lastbe != null) {
                    e1 = e.makeSymm();
                    e2 = new Edge(e.p1, nd);
                    e3 = new Edge(nd, e.p2);
                    if (comedge == null) {
                        hullStart = lastbe;
                        lastbe.nextH = e2;
                        lastbe = e2;
                    } else
                        comedge.linkSymm(e2);


                    comedge = e3;
                    triangles.add(new Triangle(edges, e1, e2, e3));
                    swapTest(e);
                }
            } else {
                if (comedge != null) break;
                lastbe = e;
            }
            e = enext;
        }

        lastbe.nextH = e3;
        e3.nextH = e;
    }

    private int searchEdge(Edge e, Node nd) {
        int f2, f3;
        Edge e0 = null;
        if ((f2 = e.nextE.onSide(nd)) == -1) {
            if (e.nextE.invE != null)
                return searchEdge(e.nextE.invE, nd);
            else {
                actE = e;
                return -1;
            }
        }
        if (f2 == 0) e0 = e.nextE;
        Edge ee = e.nextE;
        if ((f3 = ee.nextE.onSide(nd)) == -1) {
            if (ee.nextE.invE != null)
                return searchEdge(ee.nextE.invE, nd);
            else {
                actE = ee.nextE;
                return -1;
            }
        }
        if (f3 == 0) e0 = ee.nextE;
        if (e.onSide(nd) == 0) e0 = e;
        if (e0 != null) {
            actE = e0;
            if (e0.nextE.onSide(nd) == 0) {
                actE = e0.nextE;
                return 0;
            }
            if (e0.nextE.nextE.onSide(nd) == 0) return 0;
            return 1;
        }
        actE = ee;
        return 2;
    }

    private void swapTest(Edge e11) {
        Edge e21 = e11.invE;
        if (e21 == null || e21.inT == null) return;
        Edge e12 = e11.nextE;
        Edge e13 = e12.nextE;
        Edge e22 = e21.nextE;
        Edge e23 = e22.nextE;
        if (e11.inT.inCircle(e22.p2) || e21.inT.inCircle(e12.p2)) {
            e11.update(e22.p2, e12.p2);
            e21.update(e12.p2, e22.p2);
            e11.linkSymm(e21);
            e13.inT.update(e13, e22, e11);
            e23.inT.update(e23, e12, e21);
            e12.asIndex();
            e22.asIndex();
            swapTest(e12);
            swapTest(e22);
            swapTest(e13);
            swapTest(e23);
        }
    }

    private void markHull() {
        Edge e = hullStart;
        if (e != null)
            do {
                e.onHull = true;
                e.p1.onHull = true;
                e.p2.onHull = true;
                e = e.nextH;
            } while (!e.isEqual(hullStart));
    }
}
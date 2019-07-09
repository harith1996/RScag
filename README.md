# Improving the Robustness of Scagnostics

Code and dataset of the paper "Improving the Robustness of Scagnostics" sumbitted to IEEE InfoVis 2019.  

## File Description

The sourcecode is modified based on the original Java implementation (i.e., code on the Scag-06 folder) provided by Leland Wilkinson from the R package "scagnostics" of CRAN.  

### RScag

File "RScag" is the code for robust scagnostics presented in our paper.  

### Scag-05

File "Scag-05" contains the code for Scag-05 proposed on the "Graph-Theoretic Scagnostics" paper on InfoVis 2005.  

Scag-05 defines a point as outlier when the point has degree 1 and the adjacent edge of this point is larger than a threshold.  

### Scag-06

File "Scag-06" contains the code of new version of scagnostics proposed on the "High-Dimensional Visual Analytics: Interactive Exploration Guided by Pairwise Views of Point Distributions" paper on TVCG 2006.  

Scag-06 defines a point as outlier then the adjacent edges of this point are larger than a threshold. It can delete interior outliers, but it becomes less robust than Scag-05.  

### Dataset

File "dataset" are two types of data we used. "dataset/real_data" contains more than 69K real scatterplots. "dataset/data_analysis" contains the R code we used to automatically generate synthetic scatterplots with different distributions, such as binormal, clustered, exponential, funnel etc.

## Contact

If you have any problem with the code or the dataset, please contact zywangx@gmail.com.

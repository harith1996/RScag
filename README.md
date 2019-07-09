# RScag
Code and dataset of the paper "Improving the Robustness of Scagnostics" sumbitted to IEEE InfoVis 2019.  
File "RScag" is the code for robust scagnostics.  
File "Scag-05" is the code for Scag-05, which defines a point as outlier when the point has degree 1 and the adjacent edge of this point is larger than a threshold.  
File "Scag-06" is the new version of scagnostics, which defines a point as outlier then the adjacent edges of this point are larger than a threshold. It can delete interior outliers, but it becomes less robust than Scag-05.  
File "dataset" are two types of data we used. "dataset/real_data" contains more than 60M real scatterplots. "dataset/data_analysis" contains the R code we used to automatically generate different scatterplots with different distributions, such as binormal, clustered, exponential, funnel etc.

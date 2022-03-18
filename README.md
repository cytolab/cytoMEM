# cytoMEM

Marker Enrichment Modeling (MEM) is a tool designed to calculate enrichment scores.  MEM generates human and machine readable labels that quantify the features enriched in a sample.  The classic use of MEM is to identify multiple populations of cells and to compare each population to all of the other remaining cells from the original sample.  MEM enrichment scores range from +10 (meaning greatly enriched) through 0 (meaning not enriched) to -10 (meaning greatly lacking).  MEM scores are built form two fundamental statistics, the median and interquartile range, and the output of MEM can be represented as a heatmap of the scores where the rows are each population and the columns are measured features.  This information can also be represented in a compact label where the most enriched features are listed first.

### Installing cytoMEM

if (!require("BiocManager", quietly = TRUE))  
        install.packages("BiocManager")  
        BiocManager::install("MEM")

### Citation

If you use the cytoMEM package, please use the following citation:

Diggins, K., Greenplate, A., Leelatian, N. et al. Characterizing cell subsets using marker enrichment modeling. Nat Methods 14, 275–278 (2017). https://doi.org/10.1038/nmeth.4149

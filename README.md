# Conversion of steppe to agriculture increases spatial heterogeneity of soil functional genes
This repository contains code usted to generate results for the manuscript "Conversion of steppe to agriculture increases spatial heterogeneity of soil functional genes". Citation: ...

The folder **data** contains two KO tables presented in this paper, "kegg_KO_percent.xls" is compositional profile normalized by sum of all KO reads number, and "kegg_KO_profile.xls" is read count profile which is used to calculate alpha diversity. The folder **scripts** contains universal functions.

## Explanations for the code
1. vpa.r is VPA(variation partitioning analysis) code.
2. ddr.r calculates similarity and distance used in DDRs(dissimilarity distance curve).
3. randomforest.r is randomforest model mentioned in this manuscript.
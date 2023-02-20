This repository contains all analyses and figure scripts used in ''Dissection of the role of a SH3 domain in the evolution of binding preference of paralogous proteins ''. 

PCA_analysisX.R scripts are used for the analysis of raw PCA data (out put of pyphe-quant). The scripts are used to treshold, normalize and save the data for subsequent visualization. There is one script PCA_analysis for each of the three PCA experiment described in the paper.

PCA_SH3_dependency.R is used to test if the interactions are different between the paralog variant. It uses Wilcoxon tests( Benjamin-Hochberg correction) to test if the extant paralogs show a different PPI scores than other paralog variants. The significant tests between the extant paralog and the SH3-depleted variant indicate the PPIs which are SH3-dependent.

FigureX.R and FigureSX.R are scripts to visualize the results and create the figures and supplementary figures of the paper. Some supplementary panels and figures are generated in FigureX.R or in the PCA_analysisX.R scripts.

# EDC_grouping (1_Physicopathway_grouping.R)
## Preprocessing of data
1. Retrieving compounds with EDC average score > 0.9
2. Using RDkit to calculate physico chemical descriptors
3. Removing compounds with constant descriptors
4. Removing descriptors with constant values
5. Removing compounds with NA desciptor values
6. Retrieving pathway NES scores across all 15 layers 
7. Removing constant pathways across all chemicals for each data layer
## Single view HCA 
1. Scaling the descriptors between 0 and 1
2. Calculation of euclidean distance matrix for descriptors matrix
3. Performing agnes with ward method on descriptors space
4. Calculation of euclidean distance for the pathway score in each data layer
5. Performing agnes with ward method on transcriptome space for each data layer
## Multi view clustering on toxicogenomics space
1. Performing multi view clustering using ward method and euclidean distance on toxicogenomics  data space
# Regularized canonical correlation analysis (3_rCCA_independent.R)
1. Preprocessing of Rdkit descriptors and removing pathways with constant values
2. Tuning lambda for regularized canonical correlation analysis using M-FOLD CV (folds=5) for each toxicogenomics data layer-physico chemical properties pair
3. Making final r-CCA using optimized paramters for each pair
4. Visualization of r-CCA based on Correlation distance matrix for each pair as heatmap plots
5. 


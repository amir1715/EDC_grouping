# A. EDC_grouping (scripts/1_Physicopathway_grouping.R)
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
3. Performing agglomerative nested HCA using ward method on physico-chemical descriptors (chemical) space
4. Calculation of euclidean distance for the pathway score in each data layer
5. Performing agnes with ward method on transcriptome space for each data layer
## Multi view clustering on toxicogenomics space
1. Performing multi view clustering using ward method and euclidean distance on toxicogenomics  data space
# B. Regularized canonical correlation analysis (r-CCA) (scripts/3_rCCA_independent.R)
1. Preprocessing of Rdkit descriptors and removing pathways with constant values
2. Tuning of lambda for each toxicogenomics data layer-physico chemical properties pair by means of grid  M-FOLD CV (folds=5)   
3. Making final r-CCA using optimized paramters for each pair from step 2
4. Heatmap visualization of r-CCA results with descriptors in terms of Correlation distance matrix for each data layer
5. Chor diagram of top 5% most correlated descriptors with pathways scores for each data layer


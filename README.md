# EDC_grouping
## Preprocessing
1. Retrieving compounds with EDC average score more than 0.9
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
4. 




# 1. Preprocessing --------------------------------------------------------



# Data prep
load('outputData/pathway_scores.RData')
load('outputData/properties.RData')


physico_data<-RDkit_properties[,!colnames(RDkit_properties) %in% c('cas','comp_name','mesh')]
rownames(physico_data)<-RDkit_properties$mesh
# data analysis
physico_data<-physico_data[,!apply(physico_data, 2, sd,na.rm=T) ==0] # removing descriptors with constant values
source('functions/general_function.R')
physico_data<-zero_one_scale(physico_data)
patways_data=lapply(pathway_scores, function(x)x[,!apply(x, 2, sd,na.rm=T) ==0]) # removing pathaways with zero SD across all compounds


# 2. Pearson-spearman pairwise Correlation -----------------------------------------

cor_pair<-function(pat_mat,physico_data,method='pearson'){
pearson_Cor<-matrix(NA, nrow = ncol(physico_data), ncol = ncol(pat_mat))
for (i in 1:nrow(pearson_Cor)){
    for (j in 1:ncol(pearson_Cor)){
 pearson_Cor[i,j]<- cor(physico_data[,i],pat_mat[,j],use="pairwise.complete.obs", method=method)  
    }
}
return(pearson_Cor)
}

pearson_correlations<-lapply(patways_data, function(x)cor_pair(x,physico_data,method='pearson'))
spearman_correlations<-lapply(patways_data, function(x)cor_pair(x,physico_data,method='spearman'))
sapply(pearson_correlations, max)
sapply(spearman_correlations, max)
save(spearman_correlations,pearson_Cor,file='outputData/spearman_pearson_correlations.RData')













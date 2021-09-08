
# 1. Data preparation and calculation of RDKit properties -----------------------------------------------
# Retrieving the cas ids for the grouping
all_Scores<-xlsx::read.xlsx('../EDCmet_new/outputData/excel_files/all_compounds_vam_moa_15networks.xlsx',sheetIndex = 1)
cas_ids<-all_Scores$cas[which(all_Scores$average_edc_score>0.9)]
cas_ids<-unique(cas_ids)[!unique(cas_ids) %in% '']
f<-file('scripts/python/cas.txt')
writeLines(cas_ids,f)
# CALCULATION OF DESCRIPTORS with RDkit
# using RDkit and the cas.txt file as input in Solu to calculate 200 descriptors for each chemical
 setwd('scripts/python')
 system('sh rdkit2pc.sh') # a pipeline with bash to activate RDkit retireving the smiles from CAS and caclulating the features
 setwd('../..')

# 2. Preprocessing of pathway scores and RDkit properties----------------------------


load('../EDCmet_new/outputData/glm/integrated_fgsea_allcompounds_results_final_moa.RData')
remove_pattern='TG_Gates_Rat_invivo_Repeated_15_day|TG_Gates_Rat_invivo_Repeated_8_day|TG_Gates_Rat_invivo_Repeated_29_day'
fgs<-fgs[names(fgs)[-grep(x=names(fgs),pattern=remove_pattern)]]
comp_mesh_fgs<-Reduce('intersect',lapply(fgs, function(x)rownames(x$x))) # the compounds in all layers

RDkit_properties<-read.csv('scripts/python/RD_kit_res.csv') # calculated descriptors by RDkit
RDkit_properties<-RDkit_properties[-which(apply(RDkit_properties[,
                                  !colnames(RDkit_properties) %in% 'cas'],
                                   1, sd,na.rm=T) ==0),] # removing compounds with constant descriptors

RDkit_properties<-RDkit_properties[,-which(apply(RDkit_properties[,
                                  !colnames(RDkit_properties) %in% 'cas'],
                                   2, sd,na.rm=T) ==0)] # removing descriptors with constant values
RDkit_properties<-na.omit(RDkit_properties) # removing compounds with NA descriptor values
source('functions/annotation_functions.R')
RDkit_properties$mesh<-cas2mesh(RDkit_properties$cas) # Cas to mesh conversion of the compounds
#keeping compounds with pathway scores for all layers
RDkit_properties<-RDkit_properties[RDkit_properties$mesh %in% intersect(RDkit_properties$mesh,comp_mesh_fgs),] 
RDkit_properties$comp_name<-mesh2name(RDkit_properties$mesh)

pathway_scores<-lapply(fgs, function(x){y=x$x[RDkit_properties$mesh,];
                                        y=y[,!apply(y, 2, sd)==0]
                                        y}) # removing constant pathways across all chemicals

save(RDkit_properties,file = 'outputData/properties.RData')   # properties
save(pathway_scores,file = 'outputData/pathway_scores.RData') # pathway scores in 15 layers



# 3. Single view Clustering of compounds based on RDkit physico-chemical properties ----------
library(cluster)
load('outputData/properties.RData')
non_desc_cols<-c('cas','mesh','comp_name')
df<-as.matrix(RDkit_properties[,!colnames(RDkit_properties) %in% non_desc_cols])
source('functions/general_function.R')
df<-zero_one_scale(df) # scaling of the descriptors between 0-1

# distance matrix
dx<-dist(df, method = "euclidean") # euclidean distance matrix

# clustering
pc_cluster<-agnes(dx,diss = T,method = 'ward') # HCA ward method

dend<-as.dendrogram(pc_cluster)
plot(dend,main = 'Clustering by Euclidean distance')

grp_pc <- cutree(pc_cluster, k = 5) # we select 5 groups based on the plot
names(grp_pc)<-RDkit_properties$comp_name

rownames(df)<-RDkit_properties$comp_name
library(factoextra)
fviz_cluster(list(data = df, cluster = grp_pc),main = 'Physico Chemical Descriptors')
names(grp_pc)<-RDkit_properties$mesh
save(grp_pc,file = 'outputData/physico_chemical_based_grouping.RData')

# 4. Single view Clustering of compounds based on pathway scores for each data layer ---------------------------
library(cluster)
load('outputData/pathway_scores.RData')
source('functions/general_function.R')
scaled_pat_scores<-lapply(pathway_scores, function(x)zero_one_scale(x))
dx_sclaed_pat_scores<-lapply(scaled_pat_scores, function(x)dist(x,method = 'euclidean'))
save(scaled_pat_scores,dx_sclaed_pat_scores,file = 'outputData/euclidean_dist_scaled_patway_scores.RData')
##preparation of clustering denrogram plot
library(cluster)
library(ggplotify)
load('outputData/euclidean_dist_scaled_patway_scores.RData')
pat_cluster<-lapply(dx_sclaed_pat_scores, function(x)agnes(x,diss = T,method = 'ward')) 
grp<-list()
dend_plt_list<-list()

for (i in 1:length(pat_cluster)){
  layer_name<-names(pat_cluster)[[i]]
  pat_dend<-as.dendrogram(pat_cluster[[i]])
  #plot(pat_dend,main = layer_name) # one by one plot
  dend_plt_list[[i]]<-as.ggplot(~plot(pat_dend,main=layer_name,legend = F,cex.lab=2))
  grp[[i]]<-cutree(pat_cluster[[i]], k = 5)
  names(grp[[i]])<-rownames(scaled_pat_scores[[i]])
}
names(grp)<-names(pat_cluster)
save(grp,file = 'outputData/patway_based_grouping.RData')
save(dend_plt_list,file = 'outputData/dendplot_single_cluster_pathways.RData')
# final dendogram plot
load('outputData/dendplot_single_cluster_pathways.RData')
library(ggpubr)
ggarrange(plotlist = dend_plt_list,nrow = 5,ncol = 3,common.legend = T)                  

## Fviz PLOT preparation
load('outputData/patway_based_grouping.RData')
load('outputData/euclidean_dist_scaled_patway_scores.RData')
source('functions/annotation_functions.R')
comp_names<-mesh2name(names(grp$Drug_Matrix_Rat_invitro_Single_Dose_1_day)) # rownames of all layers are same
library(ggplotify)
plt_list<-list()
library(factoextra)
 for (i in 1:length(scaled_pat_scores)){
  # rownames(scaled_pat_scores[[i]])<-names(grp[[i]])<-comp_names
   pl<-fviz_cluster(list(data = scaled_pat_scores[[i]],
   cluster = grp[[i]]),main=names(scaled_pat_scores)[[i]])
   plt_list[[i]]<-as.ggplot(~plot(pl,legend = F,cex.lab=2))
   }
 save(plt_list,file = 'outputData/plots_single_cluster_pathways.RData')

 ##final fviz plot for pathway scores
load('outputData/plots_single_cluster_pathways.RData')
library(ggpubr)
ggarrange(plotlist = plt_list,nrow = 5,ncol = 3,common.legend = T)                  


# 5. Multi  view Clusteing based on the pathway scores of all 15 toxicogenomics data layers ---------------------------
load('outputData/euclidean_dist_scaled_patway_scores.RData')
library(IntClust)
library(factoextra)

dx_sclaed_pat_scores<-lapply(dx_sclaed_pat_scores, function(x)as.matrix(x))
MVC<-WonM(List = dx_sclaed_pat_scores,
          type = 'dist',
          linkage = rep('ward',15),
          nrclusters = seq(5,10,1))  #multi view clustering
 dend<-as.dendrogram(MVC$Clust)
 plot(dend,main = 'Multi view Clustering by Euclidean distance')
 grp_mvc <- cutree(MVC$Clust, k = 5) # we select 5 groups based on the plot
source('functions/annotation_functions.R')

names(grp_mvc)<- rownames(scaled_pat_scores$Drug_Matrix_Rat_invitro_Single_Dose_1_day)
save(grp_mvc,file = 'outputData/MVC_based_grouping')


# 6. Selection of compounds by Similarities between different grouping approaches ----------
# we used both HCA with hamming distance and the jaccard similarity between
# the compounds in each cluster
library(e1071)
library(cluster)
library(dendextend)
load('outputData/physico_chemical_based_grouping.RData')
load('outputData/patway_based_grouping.RData')
load('outputData/MVC_based_grouping')

grp[[length(grp)+1]]<-grp_mvc
names(grp)[[length(grp)]]<-'Multi_view_clustering'

grp[[length(grp)+1]]<-grp_pc
names(grp)[[length(grp)]]<-'physico_chemical'


freq_table<-as.data.frame(sapply(grp, table)) # frequency of compounds in each cluster 
write.csv(freq_table,file = 'outputData/excel_files/frequency_clusters.csv')
df<-t(as.data.frame(grp))
dx<-as.dist(hamming.distance(df[!rownames(df) %in%  "Multi_view_clustering", ])) #hamming distance matrix on binary data

# clustering
cluster<-agnes(dx,diss = T,method = 'ward') # HCA ward method
dend<-as.dendrogram(cluster)
dend_labels <- labels(dend)
labels(dend) <- ""
plot(dend,main='Similarities between grouping methods by hamming distance')
text(x = 1:length(dend_labels), labels = dend_labels, cex = .5,srt = 30, adj = c(1,1), xpd = T)

# jaccard similarity
jac_mat<-list()
for (k in 1:length(grp)){
conf_mat<-matrix(NA,nrow = length(unique(unlist(grp))),ncol = length(unique(unlist(grp))))
rownames(conf_mat)<-paste(names(grp)[k],'_',1:nrow(conf_mat),sep = '')
colnames(conf_mat)<-paste('Physico_chemical',1:ncol(conf_mat),sep = '')
  for (i in 1:nrow(conf_mat)){
    clust_mvc<-i
     for (j in 1:ncol(conf_mat)){
     clust_pc<-j
     g1<-names(grp$Multi_view_clustering)[which(grp[[k]]==clust_mvc)]
     g2<-names(grp$physico_chemical)[which(grp$physico_chemical==clust_pc)]
     jac<-length(intersect(g1,g2))/length(union(g1,g2))
     conf_mat[i,j]<-jac  
     } # for j
   }     # for i
 jac_mat[[k]]<-conf_mat
 }          # for k
 names(jac_mat)<-names(grp)
 jac_mat$physico_chemical<-c()
  #final selection of compounds
 max_jaccard<-as.data.frame(lapply(jac_mat[!names(jac_mat) %in% "Multi_view_clustering" ], function(x)apply(x,2,max)))
# we select the top cluster of physico chemical property clustering
 pc_clusters_overlapped_with_tox_layers<-which(apply(max_jaccard, 1,
                                                                     mean)>=quantile(apply(max_jaccard, 1, 
                                                                     mean),probs = .75))
 
 # the compounds in the top selected clusters
 selected_by_pc<-names(grp$physico_chemical)[grp$physico_chemical %in% pc_clusters_overlapped_with_tox_layers]
 
 mvc_category<-sort(apply(jac_mat$Multi_view_clustering[,
                    pc_clusters_overlapped_with_tox_layers],
                    1, mean),decreasing = T)[1]
 library(magrittr)
 temp<-names(mvc_category) %>% strsplit('_')
 selected_MVC_category<-unlist(temp)[length(unlist(temp))]
 selected_by_mvc<-names(grp$Multi_view_clustering)[grp$Multi_view_clustering%in% selected_MVC_category]
 final_list_CCA<-intersect(selected_by_mvc,selected_by_pc)



output_excel_file='outputData/excel_files/Jaccard_similarity.xlsx'

# excel file 
for (i in 1:length(jac_mat)){
                        
  to_save<-jac_mat[[i]]
  xlsx::write.xlsx(to_save,file = output_excel_file,sheetName = paste('layer_',i,
                                                    sep = ''),append = T,
                                                   row.names = T)  # appending all in one excel file 
}
save(jac_mat,final_list_CCA,file = 'outputData/selected_for_CCA.RData')



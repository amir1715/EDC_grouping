
# 1. rCCA -----------------------------------------------------------------


library(mixOmics)
# Data prep
load('outputData/pathway_scores.RData')
load('outputData/properties.RData')
load('outputData/selected_for_CCA.RData')
patways_data<-lapply(pathway_scores, function(x){y=x[final_list_CCA,];y})
physico_data<-RDkit_properties[RDkit_properties$mesh  %in% final_list_CCA,]
rownames(physico_data)<-physico_data$mesh
physico_data<-physico_data[final_list_CCA,
                           !colnames(physico_data) %in% c('cas','comp_name','mesh')]
# data analysis
physico_data<-physico_data[,!apply(physico_data, 2, sd,na.rm=T) ==0] # removing descriptors with constant values
source('functions/general_function.R')
physico_data<-zero_one_scale(physico_data)
patways_data=lapply(patways_data, function(x)x[,!apply(x, 2, sd,na.rm=T) ==0]) # removing pathaways with zero SD across all compounds

grid1 <- grid2 <- seq(0.05, 0.2, length=5)

library(doParallel)               #For paralellization 
cl<-makeCluster(length(patways_data))            
registerDoParallel(cl)           
CV_resuls<-foreach(i=1:length(patways_data),
                   .export = 'tune.rcc') %dopar% tune.rcc(physico_data,
                                                          patways_data[[i]],
                                                          grid1 = grid1,
                                                          grid2 = grid2,
                                                          validation = 'Mfold',
                                                          folds=5)
stopCluster(cl)
names(CV_resuls)<-names(patways_data)

#cim(cor(x, y), cluster = "none")
cl<-makeCluster(length(patways_data))            
registerDoParallel(cl)           
final_resuls<-foreach(i=1:length(CV_resuls),
                   .export = 'rcc') %dopar% rcc(physico_data,
                                                patways_data[[i]],
                                                ncomp=3,
                                                lambda1 =CV_resuls[[i]]$opt.lambda1,
                                                lambda2 =CV_resuls[[i]]$opt.lambda2)
                                                
stopCluster(cl)
names(final_resuls)<-names(CV_resuls)
save(final_resuls,CV_resuls,file = 'outputData/cvCCA.RData')


# 2. Plot analysis --------------------------------------------------------



load('outputData/cvCCA.RData')

graphics.off()
par("mar")
par(mar=c(1,1,1,1))
gr_list<-cor_mat<-list()
library(ggplotify)
library(mixOmics)
source('functions/general_function.R')
for (i in 1:length(final_resuls)){
  gr_list[[i]]<-ggplotify::as.ggplot(~cim(final_resuls[[i]],
                                          xlab = "",
                                          ylab = "",
                                          title = names(final_resuls)[i],
                                          #cluster='none',
                                          dist.method = c('correlation','correlation'),
                                          margins = c(2, 2)))
  
  cor_mat[[i]]<-cim_cormat(final_resuls[[i]])$mat
  }

names(cor_mat)<-names(final_resuls)
save(gr_list,cor_mat,file = 'outputData/after_HCA_graphs_CCA.RData')
#final heatmap
library(ggpubr)
load('outputData/after_HCA_graphs_CCA.RData')
pdf('outputData/plots/CCA_subsetfromHCA.pdf',width = 50,height = 50)
plot(ggarrange(plotlist = gr_list,nrow = 5,ncol = 3,common.legend = T) )
dev.off()


# 3.Chor diagram ----------------------------------------------------------
library(ggpubr)
library(circlize)
library(ggplotify)
load('outputData/after_HCA_graphs_CCA.RData')
source('functions/general_function.R')

most_related_desc_mat<-lapply(cor_mat, function(x)most_related_desc(x))

circ_plot<-list()
for (i in 1:length(most_related_desc_mat)){
mat<-most_related_desc_mat[[i]]
 grid.col <- setNames(rainbow(length(unlist(dimnames(mat)))),
                      union(rownames(mat),
                            colnames(mat)))
#chordDiagram(mat, grid.col = grid.col) 
circ_plot[[i]]<-as.ggplot(function(){ 
chordDiagram(mat, annotationTrack = "grid",preAllocateTracks = 1, grid.col = grid.col)
 circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
   xlim = get.cell.meta.data("xlim");ylim = get.cell.meta.data("ylim")
   sector.name = get.cell.meta.data("sector.index")
   circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise",
               niceFacing = TRUE, adj = c(0, 0.5),cex = .6)
   circos.axis(h = "top", 
               labels.cex = .1,
               major.tick.percentage = 0.2, 
               sector.index = sector.name, track.index = 2)
 }, bg.border = NA)
 title(names(most_related_desc_mat)[i])})
}
names(circ_plot)<-names(most_related_desc_mat)
pdf('outputData/plots/CCA_subsetCircos.pdf',width = 30,height = 30)
plot(ggarrange(plotlist = circ_plot,nrow = 5,ncol = 3,common.legend = T) )
dev.off()



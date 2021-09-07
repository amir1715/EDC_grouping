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
x=physico_data

x<-physico_data[,-which(apply(physico_data, 2, sd,na.rm=T) ==0)] # removing descriptors with constant values

y=patways_data$Drug_Matrix_Rat_invitro_Single_Dose_1_day[,-which(apply(patways_data$Drug_Matrix_Rat_invitro_Single_Dose_1_day, 2, sd,na.rm=T) ==0)]


design = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)
data<-list(pc=x,pat=y)




grid1 <- grid2 <- seq(0.05, 0.2, length=5)
#cv <- tune.rcc(x, y, grid1 = grid1, grid2 = grid2, validation = "Mfold",folds = 5)
result <- rcc(x, y, ncomp = 3, lambda1 = 0.2, lambda2 = 0.125)
#cim(cor(x, y), cluster = "none")
cim(result, xlab = "genes", ylab = "lipids", margins = c(5, 6))
circosPlot(result, cutoff = 0.7, ncol.legend = 2, size.legend = 1.1)

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


load('outputData/cvCCA.RData')
library(mixOmics)
graphics.off()
par("mar")
par(mar=c(1,1,1,1))
gr_list<-list()
for (i in 1:length(final_resuls)){
gr_list[[i]]<-ggplotify::as.ggplot(~cim(final_resuls[[i]],
                  xlab = "",
                  ylab = "",
                  title = names(final_resuls)[i],
                  margins = c(2, 2)))}
save(gr_list,file = 'outputData/graphs_CCA.RData')
library(ggpubr)
ggarrange(plotlist = gr_list,nrow = 5,ncol = 3,common.legend = T)                  



nutrimouse.sgccda <- wrapper.sgccda(X=data,
                                    Y = Y,
                                    design = design,
                                    keepX = list(gene=c(10,10), lipid=c(15,15)),
                                    ncomp = 2,
                                    scheme = "horst")

circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2, size.legend = 1.1)

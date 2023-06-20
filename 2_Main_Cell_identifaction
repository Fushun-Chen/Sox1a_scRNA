### Analyze the Sox1a  quality filtered data
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggrepel)

the.result <- readRDS("~/the.results.rds")

# sox1aAll: 1dpf 2dpf 3dpf 5dpf

Sox1a_samples <- the.result$sox1aAll
Sox1a_samples <- ScaleData(Sox1a_samples, verbose = FALSE)
saveRDS(Sox1a_samples, "~/Sox1a_samples.rds")
##
# UMAP plot
DimPlot(Sox1a_samples, reduction = "umap", group.by = "orig.ident", label = F)
DimPlot(Sox1a_samples, reduction = "umap", group.by = "integrated_snn_res.1.5", label = T)

DimPlot(Sox1a_samples, reduction = "umap", group.by = "integrated_snn_res.0.5", label = T) + DimPlot(Sox1a_samples, reduction = "umap", group.by = "integrated_snn_res.1.5", label = T)
#
Idents(Sox1a_samples) = Sox1a_samples$integrated_snn_res.1.5 
##
Sox1a_samples$Cell_labels <- "Not Classified" # Not Classified # Null
Sox1a_samples$Cell_labels[Sox1a_samples$integrated_snn_res.1.5 %in% c(0,7, 28, 13, 10,2)] <- "Neurons"
Sox1a_samples$Cell_labels[Sox1a_samples$integrated_snn_res.1.5 %in% c(8, 30, 17, 32, 14, 9)] <- "Neural progenitors"
Sox1a_samples$Cell_labels[Sox1a_samples$integrated_snn_res.1.5 %in% c(25, 16, 1, 5, 24, 15, 20, 36, 23, 33)] <- "Lateral line"
Sox1a_samples$Cell_labels[Sox1a_samples$integrated_snn_res.1.5 %in% c(19, 11, 6, 22, 18)] <- "Mesoderm"
Sox1a_samples$Cell_labels[Sox1a_samples$integrated_snn_res.1.5 %in% c(3, 4, 29)] <- "Fin fold"
Sox1a_samples$Cell_labels[Sox1a_samples$integrated_snn_res.1.5 %in% c(38, 26, 12)] <- "Hematopoietic system" # "Blood vessel" 
Sox1a_samples$Cell_labels[Sox1a_samples$integrated_snn_res.1.5 %in% c(27)] <- "Cardiomyocytes"
Sox1a_samples$Cell_labels[Sox1a_samples$integrated_snn_res.1.5 %in% c(31)] <- "Fast skeletal muscle"
Sox1a_samples$Cell_labels[Sox1a_samples$integrated_snn_res.1.5 %in% c(34)] <- "Vascular"

Sox1a_samples$Cell_labels <- factor(Sox1a_samples$Cell_labels,  levels = c("Neural progenitors", "Neurons", "Lateral line", "Mesoderm", "Fin fold", "Hematopoietic system", "Cardiomyocytes", "Fast skeletal muscle", "Vascular", "Not Classified"))

Idents(Sox1a_samples) = Sox1a_samples$Cell_labels 

DimPlot(Sox1a_samples, reduction = "umap", group.by = "Cell_labels", label = T)

DimPlot(Sox1a_samples, reduction = "umap", group.by = "Cell_labels", label = T, split.by = "orig.ident")

#
library(viridis)

p1 <- DoHeatmap(Sox1a_samples, features = c("sox2", "fabp7a", "gfap", "elavl3", "sv2a", "cldnb", "meox1", "vcanb", "ccl25b", "pitx2", "mylpfa", "fli1a"), group.by = "Cell_labels") + scale_fill_viridis(option = "E") +  theme(text = element_text(size = 17)) + theme(axis.text.y = element_text(size = 17))

pdf("~/Neurons_Progenitors_Mesoderm_markers_Heatmap2.pdf", width = 28, height = 9)
print(p1)
dev.off()



## plot figure:
Sox1a_umap = Sox1a_samples@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  cbind(cell_type = Sox1a_samples@meta.data$Cell_labels) 


p <- ggplot(Sox1a_umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) + geom_point(size = 1 , alpha =1 ) + scale_color_npg() #+ scale_color_manual(values = c(brewer.pal(9, "Set1"), "grey"))

p2 <- p  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))
p2
p3 <- p2 +         
        theme(
          legend.title = element_blank(),  
          legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=15), 
        legend.key.size=unit(1,'cm') ) +  
  guides(color = guide_legend(override.aes = list(size=5))) 

p4 <- p3 + 
  geom_segment(aes(x = min(Sox1a_umap$UMAP_1) , y = min(Sox1a_umap$UMAP_2) ,
                   xend = min(Sox1a_umap$UMAP_1) +3, yend = min(Sox1a_umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(Sox1a_umap$UMAP_1)  , y = min(Sox1a_umap$UMAP_2)  ,
                   xend = min(Sox1a_umap$UMAP_1) , yend = min(Sox1a_umap$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(Sox1a_umap$UMAP_1) +1.5, y = min(Sox1a_umap$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(Sox1a_umap$UMAP_1) -1, y = min(Sox1a_umap$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p4

pdf("~/All_Cell_UMAP_CellLabels2.pdf", width = 10)
print(p4)
dev.off()
# DoHeatmap:
library(ggplot2)
library(viridis)

group.colors = c("Neural progenitors" = "#E64B35B2", "Neurons" = "#4DBBD5B2", "Lateral line" = "#00A087B2", "Mesoderm" = "#3C5488B2", "Fin fold" = "#F39B7FB2", "Hematopoietic system" = "#8491B4B2", "Cardiomyocytes" = "#91D1C2B2", "Fast skeletal muscle" = "#DC0000B2", "Vascular" = "#7E6148B2", "Not Classified" = "#B09C85B2")
Sox1a_Doheatmap_plot <- DoHeatmap(Sox1a_samples, features = c("sox2", "fabp7a", "gfap", "elavl3", "sv2a", "cldnb", "meox1", "vcanb", "ccl25b", "pitx2", "mylpfa", "fli1a"), group.by = "Cell_labels", group.colors = c("Neural progenitors" = "#E64B35B2", "Neurons" = "#4DBBD5B2", "Lateral line" = "#00A087B2", "Mesoderm" = "#3C5488B2", "Fin fold" = "#F39B7FB2", "Hematopoietic system" = "#8491B4B2", "Cardiomyocytes" = "#91D1C2B2", "Fast skeletal muscle" = "#DC0000B2", "Vascular" = "#7E6148B2", "Not Classified" = "#B09C85B2")
) + scale_fill_viridis(option = "E") +  theme(text = element_text(size = 20)) + theme(axis.text.y = element_text(size = 25)) + scale_colour_manual(values=c("Neural progenitors" = "#E64B35B2", "Neurons" = "#4DBBD5B2", "Lateral line" = "#00A087B2", "Mesoderm" = "#3C5488B2", "Fin fold" = "#F39B7FB2", "Hematopoietic system" = "#8491B4B2", "Cardiomyocytes" = "#91D1C2B2", "Fast skeletal muscle" = "#DC0000B2", "Vascular" = "#7E6148B2", "Not Classified" = "#B09C85B2"), 
                      labels=c("Neural progenitors", "Neurons", "Lateral line", "Mesoderm", "Fin fold", "Hematopoietic system", "Cardiomyocytes", "Fast skeletal muscle", "Vascular", "Not Classified")) 

pdf("~/All_Cell_DoHeatmap_CellLabels.pdf", width = 30, height = 9)
print(Sox1a_Doheatmap_plot)
dev.off()
# raster = FALSE, If true, plot with geom_raster, else use geom_tile.
Sox1a_Doheatmap_plot <- DoHeatmap(Sox1a_samples, raster = FALSE, features = c("sox2", "fabp7a", "gfap", "elavl3", "sv2a", "cldnb", "meox1", "vcanb", "ccl25b", "pitx2", "mylpfa", "fli1a"), group.by = "Cell_labels", group.colors = c("Neural progenitors" = "#E64B35B2", "Neurons" = "#4DBBD5B2", "Lateral line" = "#00A087B2", "Mesoderm" = "#3C5488B2", "Fin fold" = "#F39B7FB2", "Hematopoietic system" = "#8491B4B2", "Cardiomyocytes" = "#91D1C2B2", "Fast skeletal muscle" = "#DC0000B2", "Vascular" = "#7E6148B2", "Not Classified" = "#B09C85B2")
) + scale_fill_viridis(option = "E") +  theme(text = element_text(size = 20)) + theme(axis.text.y = element_text(size = 25)) + scale_colour_manual(values=c("Neural progenitors" = "#E64B35B2", "Neurons" = "#4DBBD5B2", "Lateral line" = "#00A087B2", "Mesoderm" = "#3C5488B2", "Fin fold" = "#F39B7FB2", "Hematopoietic system" = "#8491B4B2", "Cardiomyocytes" = "#91D1C2B2", "Fast skeletal muscle" = "#DC0000B2", "Vascular" = "#7E6148B2", "Not Classified" = "#B09C85B2"), 
                      labels=c("Neural progenitors", "Neurons", "Lateral line", "Mesoderm", "Fin fold", "Hematopoietic system", "Cardiomyocytes", "Fast skeletal muscle", "Vascular", "Not Classified")) 

pdf("~/All_Cell_DoHeatmap_CellLabels_geom_tile.pdf", width = 30, height = 9)
print(Sox1a_Doheatmap_plot)
dev.off()




# 
library(gridExtra)
library(ggeasy)

GFP_plot <- FeaturePlot(Sox1a_samples, "GFP", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
Sox1a_plot  <- FeaturePlot(Sox1a_samples, "sox1a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes() 
Sox2_plot <- FeaturePlot(Sox1a_samples, "sox2", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
Elavl3_plot  <- FeaturePlot(Sox1a_samples, "elavl3", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
Cldnb_plot  <- FeaturePlot(Sox1a_samples, "cldnb", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
Sv2a_plot  <- FeaturePlot(Sox1a_samples, "sv2a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()


pdf("~/NP_GFP_Sox1a_FeaturePlot.pdf", width = 60, height = 30)
#print(GFP_plot + Sox1a_plot)
grid.arrange(Sox2_plot, Elavl3_plot, Cldnb_plot, Sv2a_plot, GFP_plot,Sox1a_plot, ncol=3, nrow=2)

dev.off()


png("~/NP_GFP_Sox1a_FeaturePlot.png", width = 480 * 9, height = 480 * 4.5)
grid.arrange(Sox2_plot, Elavl3_plot, Cldnb_plot, Sv2a_plot, GFP_plot,Sox1a_plot, ncol=3, nrow=2)

dev.off()
#
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
##
## DimPlot by Stage

pdf("~/All_Cells_DimPlot_Stage.pdf", width = 15, height = 10)
print(DimPlot(Sox1a_samples, label = FALSE, group.by = "orig.ident", cols = c("#999999", "#333333", "#ff0000", "#00cdcd")) + ggtitle(""))
dev.off()

# res = 0.5
pdf("~/All_Cells_DimPlot_res0.5.pdf", width = 15, height = 10)
print(DimPlot(Sox1a_samples, label = TRUE, group.by = "integrated_snn_res.0.5") + ggtitle(""))
dev.off()

# res = 1.5
pdf("~/All_Cells_DimPlot_res1.5.pdf", width = 15, height = 10)
print(DimPlot(Sox1a_samples, label = TRUE, group.by = "integrated_snn_res.1.5") + ggtitle(""))
dev.off()
#  pt.size = 3
library(gridExtra)
library(ggeasy)

pdf("~/All_Cells_DimPlot_Stage_pt3.pdf", width = 15, height = 10)
print(DimPlot(Sox1a_samples, label = FALSE, pt.size = 3, group.by = "orig.ident", cols = c("#999999", "#333333", "#ff0000", "#00cdcd"))  + easy_remove_axes() + theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=20)) +  ggtitle(""))
dev.off()

# res = 0.5
pdf("~/All_Cells_DimPlot_res0.5_pt3.pdf", width = 15, height = 10)
print(DimPlot(Sox1a_samples, label = TRUE, pt.size = 3, label.size = 10, group.by = "integrated_snn_res.0.5") +  easy_remove_axes() + theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=20)) +  ggtitle(""))
dev.off()

# res = 1.5
pdf("~/All_Cells_DimPlot_res1.5_pt3.pdf", width = 15, height = 10)
print(DimPlot(Sox1a_samples, label = TRUE, pt.size = 3, label.size = 10,group.by = "integrated_snn_res.1.5") +  easy_remove_axes() + theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=20)) +  ggtitle(""))
dev.off()
# the bar plot for cell numbers at every stage:
library(ggplot2)

df_cell_num <- data.frame(table(Sox1a_samples$orig.ident))
colnames(df_cell_num) <- c("Stages", "Num")

p_bar_cell_num <- ggplot(df_cell_num, aes(x=Stages, y=Num, fill=Stages))+ geom_bar(stat="identity")+ scale_fill_manual(values=c("#999999", "#333333", "#ff0000", "#00cdcd")) + geom_text(aes(label=Num), vjust=-0.25, position = position_dodge(0.9), size=5) + theme_minimal() + xlab("") + ylab("Cell Numbers") + theme(legend.position = "none")
pdf("~/All_Cells_Number_bar.pdf")
print(p_bar_cell_num)
dev.off()
## cell number and cell types
df_cell_number_cell_types <- data.frame(table( Sox1a_samples$orig.ident, Sox1a_samples$Cell_labels))
colnames(df_cell_number_cell_types) <- c("Stages", "Cell_labels", "Num")
df_cell_number_cell_types$Cell_labels <- factor(df_cell_number_cell_types$Cell_labels,  levels = c("Neural progenitors", "Neurons", "Lateral line", "Mesoderm", "Fin fold", "Hematopoietic system", "Cardiomyocytes", "Fast skeletal muscle", "Vascular", "Not Classified"))

p_bar_cellNum_cellTypes <- ggplot(data=df_cell_number_cell_types, aes(x=Stages, y=Num, fill=Cell_labels)) +
  geom_bar(stat="identity", position=position_dodge())+ scale_fill_manual(values = c("Neural progenitors" = "#E64B35B2", "Neurons" = "#4DBBD5B2", "Lateral line" = "#00A087B2", "Mesoderm" = "#3C5488B2", "Fin fold" = "#F39B7FB2", "Hematopoietic system" = "#8491B4B2", "Cardiomyocytes" = "#91D1C2B2", "Fast skeletal muscle" = "#DC0000B2", "Vascular" = "#7E6148B2", "Not Classified" = "#B09C85B2")) + 
  geom_text(aes(label=Num), vjust=-0.25,
            position = position_dodge(0.9), size=3.5)+ theme_minimal() + xlab("") + ylab("Cell Numbers")

pdf("~/All_Cells_Number_CellTypes_bar.pdf", width = 15)
print(p_bar_cellNum_cellTypes)
dev.off()

p_bar_cellPercent_cellTypes <- ggplot(data=df_cell_number_cell_types, aes(x=Stages, y=Num, fill=Cell_labels)) +
  geom_bar(stat="identity", position="fill")+ scale_fill_manual(values = c("Neural progenitors" = "#E64B35B2", "Neurons" = "#4DBBD5B2", "Lateral line" = "#00A087B2", "Mesoderm" = "#3C5488B2", "Fin fold" = "#F39B7FB2", "Hematopoietic system" = "#8491B4B2", "Cardiomyocytes" = "#91D1C2B2", "Fast skeletal muscle" = "#DC0000B2", "Vascular" = "#7E6148B2", "Not Classified" = "#B09C85B2")) + theme_minimal() + xlab("") + ylab("Cell Percentage")
pdf("~/All_CellsPercent_CellTypes_bar.pdf", width = 10)
print(p_bar_cellPercent_cellTypes)
dev.off()
# contribution of each stage to the cell types identified
p_bar_cellPercent_cellTypes2 <- ggplot(data=df_cell_number_cell_types, aes(x=Cell_labels, y=Num, fill=Stages)) +
  geom_bar(stat="identity", position="fill")+ scale_fill_manual(values=c("#999999", "#333333", "#ff0000", "#00cdcd")) + theme_minimal() + xlab("") + ylab("Cell Percentage")
pdf("~/All_CellsStagesPercent_CellTypes_bar.pdf", width = 15)
print(p_bar_cellPercent_cellTypes2)
dev.off()

## Marker gene table at res 1.5
## FindMarkers 
Idents(Sox1a_samples) <- Sox1a_samples$integrated_snn_res.1.5
All_cell_res1.5_obj.all.mrkrs = FindAllMarkers(Sox1a_samples, only.pos = T)
#  save markers
All_cell_res1.5_obj.all.mrkrs.list <- list()
for(i in 0:39){
 
    All_cell_res1.5_obj.all.mrkrs.list[[paste("cluster", i, sep = "")]] <- All_cell_res1.5_obj.all.mrkrs[All_cell_res1.5_obj.all.mrkrs$cluster == i, ]

}

writexl::write_xlsx(All_cell_res1.5_obj.all.mrkrs.list, "~/All_Cells_res1.5_Cluster_markers.xlsx")

### Cell Label to find markers
Idents(Sox1a_samples) <- Sox1a_samples$Cell_labels

All_cell_CellType_obj.all.mrkrs = FindAllMarkers(Sox1a_samples, only.pos = T)

All_Cell_CellType.mrkrs.list <- list()
#for(i in sort(as.character(unique(Cell_NP_obj.all.mrkrs$cluster)))){
for(i in as.character(unique(sort(Sox1a_samples$Cell_labels)))){
 
    All_Cell_CellType.mrkrs.list[[i]] <- All_cell_CellType_obj.all.mrkrs[All_cell_CellType_obj.all.mrkrs$cluster == i, ]

}



writexl::write_xlsx(All_Cell_CellType.mrkrs.list, "~/All_Cells_CellType_Cluster_markers.xlsx")


# add Ensembl
features <- read.table("~/raw_feature_bc_matrix/features.tsv.gz", header = FALSE, sep = "\t")

getUniqs = function(xx){
  idx2 = duplicated(xx)
  for(i in seq(1, length(idx2),1)) if(idx2[i]) xx[i] = paste0(xx[i],"#",i)
  return(xx)
}

features$gene <- getUniqs(features$V2)
colnames(features) <- c("Ensembl", "V2", "V3", "gene")
#
All_Cell_CellType.mrkrs.list <- list()

for(i in as.character(unique(sort(Sox1a_samples$Cell_labels)))){

    One_cell_CellType_obj.all.mrkrs  <- All_cell_CellType_obj.all.mrkrs[All_cell_CellType_obj.all.mrkrs$cluster == i, ]
    One_cell_CellType_obj.all.mrkrs.Ensembl <- merge(One_cell_CellType_obj.all.mrkrs, features[, c("Ensembl", "gene")], by = "gene", sort = FALSE)
    All_Cell_CellType.mrkrs.list[[i]] <- One_cell_CellType_obj.all.mrkrs.Ensembl

}

writexl::write_xlsx(All_Cell_CellType.mrkrs.list, "~/All_Cells_CellType_Cluster_markers.xlsx")

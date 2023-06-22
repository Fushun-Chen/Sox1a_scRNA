
## Neuron and Progenitors

library(Seurat)

the.results = readRDS(paste0("~/New_integration/Step1/Step1multiRDAfiles/Progenitor_Neuronal_population_integration.rds"))

NP_obj = the.results$all_samples
DimPlot(NP_obj, label =T)
#
FeaturePlot(NP_obj, features = c("GFP", "sox1a"))
VlnPlot(NP_obj, features = c("GFP", "sox1a"))
##
# the bar plot for cell numbers at every stage:
library(ggplot2)

df_cell_num <- data.frame(table(NP_obj$orig.ident))
colnames(df_cell_num) <- c("Stages", "Num")

p_bar_cell_num <- ggplot(df_cell_num, aes(x=Stages, y=Num, fill=Stages))+ geom_bar(stat="identity")+ scale_fill_manual(values=c("#999999", "#333333", "#ff0000", "#00cdcd")) + geom_text(aes(label=Num), vjust=-0.25, position = position_dodge(0.9), size=5) + theme_minimal() + xlab("") + ylab("Cell Numbers") + theme(legend.position = "none")
pdf("~/NP_Cells_Number_bar.pdf")
print(p_bar_cell_num)
dev.off()

#
FeaturePlot(NP_obj, c("tph2", "fev", "nkx2.2a", "sst1.1", "urp1","zic1", "olig2", "ccna2", "slc6a5"))

NP_obj$Cell_labels <- as.character(NP_obj$integrated_snn_res.2)
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("0", "1", "2", "4", "10", "24")] <- "KA'"
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("5", "8", "12", "22")] <- "KA\""
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("16", "25", "30", "29", "23")] <- "V2s"
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("21")] <- "OPC"
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("18", "19")] <- "Zic+"
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("9")] <- "LFP"
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("14", "15", "32")] <- "ISN"
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("13")] <- "ERG" 
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("11")] <- "NPC" # neuronal precursor cells
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("17")] <- "GLU IN" # glutamatergic interneurons
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("26")] <- "V2-pre" # V2 interneuron precursor
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("20")] <- "NPro4" # neural progenitor in proliferation
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("3")] <- "NPro1" # neural progenitor in proliferation
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("7")] <- "NPro3" # neural progenitor in proliferation
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("6")] <- "NPro2" # neural progenitor in proliferation
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("27")] <- "MN" # motor neuron
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("28")] <- "MC" # mesenchymal cells
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("31")] <- "KC" # Keratynocytes


##
NP_obj <- ScaleData(NP_obj, verbose = FALSE)


library(viridis)
library(ggplot2)
NP_obj$Cell_labels <- factor(NP_obj$Cell_labels, levels = c("KA'", "KA\"", "V2s", "OPC", "Zic+", "LFP",  "ISN", "ERG", "NPC", "GLU IN", "V2-pre", "NPro1", "NPro2", "NPro3", "NPro4", "MN", "MC","KC"))
Idents(NP_obj) = NP_obj$Cell_labels 

p1 <- DoHeatmap(NP_obj, features = c("sst1.1", "urp1", "slc6a5", "nkx1.2lb", "olig2", "sox10", "zic1", "zic2a", "nkx2.9", "sulf1", "fev","lmx1bb", "tph2",  "gfap", "fabp7a","foxj1a", "neurod4", "neurog1", "slc17a6b", "slc17a6a", "foxn4", "vsx1", "mki67", "sox2", "sox19a", "fabp7a", "her2", "sox2", "isl1", "isl2a", "twist1b", "epcam"), group.by = "Cell_labels") + scale_fill_viridis(option = "E") +  theme(text = element_text(size = 17)) + theme(axis.text.y = element_text(size = 17))

pdf("~/Neurons_Progenitors_Mesoderm_markers_Heatmap.pdf", width = 29, height = 10)
print(p1)
dev.off()


saveRDS(NP_obj, "~/NP_obj_label.rds")
#

DimPlot(NP_obj, label =T) + DimPlot(NP_obj, label =T, group.by = "Cell_labels")

#00cdcd

pdf("~/Neurons_Progenitors_DimPlot_Stage.pdf", width = 15, height = 10)
print(DimPlot(NP_obj, label =T, group.by = "orig.ident", cols = c("#999999", "#333333", "#ff0000", "#00cdcd") ))
dev.off()


pdf("~/Neurons_Progenitors_DimPlot.pdf", width = 15, height = 10)
print(DimPlot(NP_obj, label =T, group.by = "Cell_labels"))
dev.off()


pdf("~/Neurons_Progenitors_DimPlot_split.pdf", width = 15 * 1.5, height = 10)
print(DimPlot(NP_obj, label =F, group.by = "Cell_labels", split.by = "orig.ident"))
dev.off()

Idents(NP_obj) = NP_obj$Cell_labels 

##

##

DefaultAssay(NP_obj) <- "RNA"


## FindMarkers 
NP_obj.all.mrkrs = FindAllMarkers(NP_obj, only.pos = T)
saveRDS(NP_obj.all.mrkrs, file = "~/NP_obj.all.mrkrs.rds")

#  save markers
NP_obj.all.mrkrs.list <- list()
for(i in 0:32){
 
    NP_obj.all.mrkrs.list[[paste("cluster", i, sep = "")]] <- NP_obj.all.mrkrs[NP_obj.all.mrkrs$cluster == i, ]

}

writexl::write_xlsx(NP_obj.all.mrkrs.list, "~/Sox1a_data_Cluster_markers.xlsx")

### Cell Label to find markers
Cell_NP_obj.all.mrkrs = FindAllMarkers(NP_obj, only.pos = T)
saveRDS(Cell_NP_obj.all.mrkrs, file = "~/NP_Cell_type_obj.all.mrkrs_new.rds")

Cell_NP_obj.all.mrkrs.list <- list()
#for(i in sort(as.character(unique(Cell_NP_obj.all.mrkrs$cluster)))){
for(i in as.character(unique(sort(NP_obj$Cell_labels)))){
 
    Cell_NP_obj.all.mrkrs.list[[i]] <- Cell_NP_obj.all.mrkrs[Cell_NP_obj.all.mrkrs$cluster == i, ]

}

Cell_NP_obj.all.mrkrs.list_save <- Cell_NP_obj.all.mrkrs.list
names(Cell_NP_obj.all.mrkrs.list_save) <- c("KA1", "KA2", names(Cell_NP_obj.all.mrkrs.list)[3:18])


writexl::write_xlsx(Cell_NP_obj.all.mrkrs.list_save, "~/Sox1a_data_Cell_label_Cluster_markers_new.xlsx")
# add Ensembl
features <- read.table("~/features.tsv.gz", header = FALSE, sep = "\t")

getUniqs = function(xx){
  idx2 = duplicated(xx)
  for(i in seq(1, length(idx2),1)) if(idx2[i]) xx[i] = paste0(xx[i],"#",i)
  return(xx)
}

features$gene <- getUniqs(features$V2)
colnames(features) <- c("Ensembl", "V2", "V3", "gene")
#
Cell_NP_obj.all.mrkrs.list <- list()

for(i in as.character(unique(sort(Cell_NP_obj.all.mrkrs$cluster)))){

    One_Cell_NP_obj.all.mrkrs <- Cell_NP_obj.all.mrkrs[Cell_NP_obj.all.mrkrs$cluster == i, ]
    One_Cell_NP_obj.all.mrkrs.Ensembl <- merge(One_Cell_NP_obj.all.mrkrs, features[, c("Ensembl", "gene")], by = "gene", sort = FALSE)
    Cell_NP_obj.all.mrkrs.list[[i]] <- One_Cell_NP_obj.all.mrkrs.Ensembl

}

names(Cell_NP_obj.all.mrkrs.list) <- c("KA1", "KA2", names(Cell_NP_obj.all.mrkrs.list)[3:18])

writexl::write_xlsx(Cell_NP_obj.all.mrkrs.list, "~/Sox1a_data_Cell_label_Cluster_markers_new.xlsx")


## top 5 markers 
sox1a_top5_markers <- c()
#for(i in sort(as.character(unique(Cell_NP_obj.all.mrkrs$cluster)))){
for(i in as.character(unique(sort(NP_obj$Cell_labels)))){

    sox1a_top5_markers <- c(sox1a_top5_markers, Cell_NP_obj.all.mrkrs.list[[i]]$gene[1:5])
}

pdf("~/Sox1a_data_Cell_Type_Top5_markers3.pdf", width = 22, height = 5)

print(DotPlot(NP_obj, features = unique(sox1a_top5_markers), group.by = "Cell_labels") + scale_color_viridis_c(direction = 1) + cowplot::theme_cowplot() + theme(axis.line  = element_blank())  + ylab('') + theme(axis.ticks = element_blank()) + theme(axis.text.x = element_text(angle = 45,  hjust=1)))
dev.off()

pdf("~/Sox1a_data_Cell_Type_markers_split2.pdf", width = 25, height = 15)

print(DotPlot(NP_obj, features = unique(sox1a_top5_markers), cols = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#80B1D3"), split.by = "orig.ident", group.by = "Cell_labels")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)))

dev.off()

### transcription factors
library(Seurat)
library(ggplot2)
## read TF data
TF_DB <- read.table("~/AnimalTFDB/Danio_rerio_TF.txt", header = T, sep = "\t")
##
obj.list <- NP_obj
all.mrkrs <- Cell_NP_obj.all.mrkrs

all.mrkrs_order2 <- all.mrkrs[order(all.mrkrs$cluster, all.mrkrs$avg_log2FC, decreasing = T),]

## TF
all.mrkrs_order2_TF <- all.mrkrs_order2[all.mrkrs_order2$gene %in% TF_DB$Symbol, ]

# save TF

NP_obj.all.mrkrs.TF.list <- list()
#for(i in sort(as.character(unique(all.mrkrs_order2_TF$cluster)))){
for(i in as.character(unique(sort(NP_obj$Cell_labels)))){
 
    NP_obj.all.mrkrs.TF.list[[i]] <- all.mrkrs_order2_TF[all.mrkrs_order2_TF$cluster == i, ]

}

names(NP_obj.all.mrkrs.TF.list) <- c("KA1", "KA2", names(NP_obj.all.mrkrs.TF.list)[3:18])

writexl::write_xlsx(NP_obj.all.mrkrs.TF.list, "~/TFs/Sox1a_data_Cell_type_markers_TF.xlsx")
# 
# filter p_val_adj < 0.05
all.mrkrs_order2_TF_sel <- all.mrkrs_order2_TF[all.mrkrs_order2_TF$p_val_adj < 0.05, ]
## TF
# Top 5

TFs_top5_list <- list()

for(i in as.character(unique(sort(NP_obj$Cell_labels)))){

    TFs_top5_list[[as.character(i)]] <- na.omit( all.mrkrs_order2_TF_sel[all.mrkrs_order2_TF_sel$cluster == i, ][1:5,])

}

TFs_top5 <- do.call(rbind, TFs_top5_list)
TFs_top5_sel <- TFs_top5[!is.na(TFs_top5$gene), ]
TFs_Top5_print <- DotPlot(NP_obj, features = unique(TFs_top5_sel$gene), group.by = "Cell_labels" ) + scale_color_viridis_c(direction = -1) + cowplot::theme_cowplot() + theme(axis.line  = element_blank())  + ylab('') + theme(axis.ticks = element_blank()) + coord_flip()


pdf("~/TFs_DotPlot_Cluster_p0.05_top5_2.pdf", width = 16, height = 17)

print(TFs_Top5_print)
dev.off()


TFs_Top5_print_v <- DotPlot(NP_obj, features = unique(TFs_top5_sel$gene), group.by = "Cell_labels") + scale_color_viridis_c(direction = 1) + cowplot::theme_cowplot() + theme(axis.line  = element_blank())  + ylab('') + theme(axis.ticks = element_blank()) + theme(axis.text.x = element_text(angle = 45,  hjust=1))

pdf("~/TFs_DotPlot_Cluster_p0.05_top5_V3.pdf", width = 22, height = 5)
print(TFs_Top5_print_v)
dev.off()

TFs_Top5_print_v_apart <- DotPlot(NP_obj, features = c("foxa", "nkx2.2a", "nkx2.2b", "fev","lmx1bb", "lmx1ba"), group.by = "Cell_labels") + scale_color_viridis_c(direction = -1) + cowplot::theme_cowplot() + theme(axis.line  = element_blank())  + ylab('') + theme(axis.ticks = element_blank()) + theme(axis.text.x = element_text(angle = 0,  hjust=1))



## plot  figure: 
library(ggplot2)
library(tidyverse)
library(ggsci)
my36colors <- c('#58A4C3', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#E5D2DD', '#C5DEBA', '#5F3D69',   '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')

my36colors2 <- c('#58A4C3', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#E95C59', '#476D87', '#57C3F3', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#E5D2DD', '#C5DEBA', '#5F3D69',   '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')


Sox1a_umap = NP_obj@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  cbind(cell_type = NP_obj@meta.data$Cell_labels) 


p <- ggplot(Sox1a_umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) + geom_point(size = 1 , alpha =1 ) + scale_color_aaas() + scale_color_manual(values = my36colors) #+ scale_color_manual(values = c(brewer.pal(9, "Set1"), "grey"))

#ggplot(Sox1a_umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) + geom_point(size = 1 , alpha =1 ) + scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"))
p2 <- p  +
  theme(panel.grid.major = element_blank(), #
        panel.grid.minor = element_blank(), #
        panel.border = element_blank(), #
        axis.title = element_blank(),  #
        axis.text = element_blank(), # 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))
p2
p3 <- p2 +         
        theme(
          legend.title = element_blank(), #
          legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=15), #
        legend.key.size=unit(1,'cm') ) +  # 
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

pdf("~/NP_UMAP_CellLabels2.pdf", width = 10)
print(p4)
dev.off()
png("~/NP_UMAP_CellLabels2.png", width = 480 * 2, height = 480 * 1.5 )
print(p4)
dev.off()

pdf("~/NP_UMAP_CellLabels3.pdf", width = 10)
print(p3)
dev.off()


# split by stage
pdf("~/Neurons_Progenitors_DimPlot_split_Stage2.pdf", width = 15 * 2.5, height = 12)
print(DimPlot(NP_obj, label =F, split.by = "orig.ident",pt.size = 2 ) + scale_color_manual(values = my36colors) + theme(text = element_text(size = 30), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=30))) 
dev.off()

png("~/Neurons_Progenitors_DimPlot_split_Stage.png", width = 480 * 4.5, height = 480 * 1.5)
print(DimPlot(NP_obj, label =F, split.by = "orig.ident",pt.size = 2 ) + scale_color_manual(values = my36colors) + theme(text = element_text(size = 30), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=30))) 
dev.off()
# split fev and tph2 expression by stage
pdf("~/Fev_Tph2_FeaturePlot_split_Stage.pdf", width = 15 * 2.5, height = 12)
FeaturePlot(NP_obj, c("fev", "tph2" ), split.by = "orig.ident", pt.size = 2) 
dev.off()


## stage

pdf("~/Neurons_Progenitors_DimPlot_Stage.pdf", width = 15, height = 10)
print(DimPlot(NP_obj, label =F, group.by = "orig.ident", pt.size = 2, cols = c("#999999", "#333333", "#ff0000", "#00cdcd") )+ easy_remove_axes() + theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=20)) +  ggtitle(""))
dev.off()

png("~/Neurons_Progenitors_DimPlot_Stage.png", width = 480 * 3, height = 480 * 2)
print(DimPlot(NP_obj, label =F, group.by = "orig.ident", pt.size = 2, cols = c("#999999", "#333333", "#ff0000", "#00cdcd") )+ easy_remove_axes() + theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=20)) +  ggtitle(""))
dev.off()

# DoHeatmap:
library(viridis)
library(ggplot2)

group.colors = c("KA'" = "#58A4C3", "KA\"" = '#53A85F', "V2s" = '#F1BB72', "OPC" = '#F3B1A0', "Zic+" = '#D6E7A3', "LFP" = '#57C3F3', "ISN" = '#476D87', "ERG" = '#E95C59', "NPC" = '#E59CC4', "GLU IN" = '#AB3282', "V2-pre" = '#23452F', "NPro1" = '#BD956A', "NPro2" = '#8C549C', "NPro3" = '#585658', "NPro4" = '#9FA3A8', "MN" = '#E0D4CA', "MC" = '#E5D2DD',"KC" = '#C5DEBA')

NP_Doheatmap_plot <- DoHeatmap(NP_obj, size = 10,features = c("sst1.1", "urp1", "slc6a5", "nkx1.2lb", "olig2", "sox10", "zic1", "zic2a", "nkx2.9", "sulf1", "fev","lmx1bb", "tph2",  "gfap", "fabp7a","foxj1a", "neurod4", "neurog1", "slc17a6b", "slc17a6a", "foxn4", "vsx1", "mki67", "sox2", "sox19a", "fabp7a", "her2", "sox2", "isl1", "isl2a", "twist1b", "epcam"), group.by = "Cell_labels", group.colors = c("KA'" = "#58A4C3", "KA\"" = '#53A85F', "V2s" = '#F1BB72', "OPC" = '#F3B1A0', "Zic+" = '#D6E7A3', "LFP" = '#57C3F3', "ISN" = '#476D87', "ERG" = '#E95C59', "NPC" = '#E59CC4', "GLU IN" = '#AB3282', "V2-pre" = '#23452F', "NPro1" = '#BD956A', "NPro2" = '#8C549C', "NPro3" = '#585658', "NPro4" = '#9FA3A8', "MN" = '#E0D4CA', "MC" = '#E5D2DD',"KC" = '#C5DEBA')
) + scale_fill_viridis(option = "E") +  theme(text = element_text(size = 40), legend.key.height = unit(2, 'cm')) + theme(axis.text.y = element_text(size = 40)) + scale_colour_manual(values=c("KA'" = "#58A4C3", "KA\"" = '#53A85F', "V2s" = '#F1BB72', "OPC" = '#F3B1A0', "Zic+" = '#D6E7A3', "LFP" = '#57C3F3', "ISN" = '#476D87', "ERG" = '#E95C59', "NPC" = '#E59CC4', "GLU IN" = '#AB3282', "V2-pre" = '#23452F', "NPro1" = '#BD956A', "NPro2" = '#8C549C', "NPro3" = '#585658', "NPro4" = '#9FA3A8', "MN" = '#E0D4CA', "MC" = '#E5D2DD',"KC" = '#C5DEBA'), labels=c("KA'", "KA\"", "V2s", "OPC", "Zic+", "LFP",  "ISN", "ERG", "NPC", "GLU IN", "V2-pre", "NPro1", "NPro2", "NPro3", "NPro4", "MN", "MC","KC")) + guides(color="none")

pdf("~/Neurons_Progenitors_DoHeatmap_CellLabels2.pdf", width = 55, height = 20)
print(NP_Doheatmap_plot)
dev.off()
png("~/Neurons_Progenitors_DoHeatmap_CellLabels.png", width = 480 * 8, height = 480 * 3)
print(NP_Doheatmap_plot)
dev.off()
## raster = FALSE, If true, plot with geom_raster, else use geom_tile.


NP_Doheatmap_plot <- DoHeatmap(NP_obj, size = 10, raster = FALSE ,features = c("sst1.1", "urp1", "slc6a5", "nkx1.2lb", "olig2", "sox10", "zic1", "zic2a", "nkx2.9", "sulf1", "fev","lmx1bb", "tph2",  "gfap", "fabp7a","foxj1a", "neurod4", "neurog1", "slc17a6b", "slc17a6a", "foxn4", "vsx1", "mki67", "sox2", "sox19a", "fabp7a", "her2", "sox2", "isl1", "isl2a", "twist1b", "epcam"), group.by = "Cell_labels", group.colors = c("KA'" = "#58A4C3", "KA\"" = '#53A85F', "V2s" = '#F1BB72', "OPC" = '#F3B1A0', "Zic+" = '#D6E7A3', "LFP" = '#57C3F3', "ISN" = '#476D87', "ERG" = '#E95C59', "NPC" = '#E59CC4', "GLU IN" = '#AB3282', "V2-pre" = '#23452F', "NPro1" = '#BD956A', "NPro2" = '#8C549C', "NPro3" = '#585658', "NPro4" = '#9FA3A8', "MN" = '#E0D4CA', "MC" = '#E5D2DD',"KC" = '#C5DEBA')
) + scale_fill_viridis(option = "E") +  theme(text = element_text(size = 40), legend.key.height = unit(2, 'cm')) + theme(axis.text.y = element_text(size = 40)) + scale_colour_manual(values=c("KA'" = "#58A4C3", "KA\"" = '#53A85F', "V2s" = '#F1BB72', "OPC" = '#F3B1A0', "Zic+" = '#D6E7A3', "LFP" = '#57C3F3', "ISN" = '#476D87', "ERG" = '#E95C59', "NPC" = '#E59CC4', "GLU IN" = '#AB3282', "V2-pre" = '#23452F', "NPro1" = '#BD956A', "NPro2" = '#8C549C', "NPro3" = '#585658', "NPro4" = '#9FA3A8', "MN" = '#E0D4CA', "MC" = '#E5D2DD',"KC" = '#C5DEBA'), labels=c("KA'", "KA\"", "V2s", "OPC", "Zic+", "LFP",  "ISN", "ERG", "NPC", "GLU IN", "V2-pre", "NPro1", "NPro2", "NPro3", "NPro4", "MN", "MC","KC")) + guides(color="none")


pdf("~/Neurons_Progenitors_DoHeatmap_CellLabels_geom_tile2.pdf", width = 55, height = 20)
print(NP_Doheatmap_plot)
dev.off()
# 
library(gridExtra)
library(ggeasy)

Sox2_plot <- FeaturePlot(NP_obj, "sox2", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
Sv2a_plot  <- FeaturePlot(NP_obj, "sv2a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
Sox1a_plot  <- FeaturePlot(NP_obj, "sox1a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()


pdf("~/NP_Sox2_Sv2a_FeaturePlot.pdf", width = 20* 2, height = 15)
#print(GFP_plot + Sox1a_plot)
grid.arrange(Sox2_plot, Sv2a_plot, ncol=2, nrow=1)

dev.off()


png("~/NP_Sox2_Sv2a_FeaturePlot.png", width = 800 * 2.5, height = 800)
grid.arrange(Sox2_plot, Sv2a_plot, ncol=2, nrow=1)

dev.off()
#
pdf("~/NP_Sox2_Sv2a_Sox1a_FeaturePlot2.pdf", width = 20* 3, height = 20)
#print(GFP_plot + Sox1a_plot)
grid.arrange(Sox2_plot, Sv2a_plot, Sox1a_plot, ncol=3, nrow=1)

dev.off()

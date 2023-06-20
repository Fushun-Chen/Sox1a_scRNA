# using filter LFP and ISN data to try
library(Seurat)
library(monocle)

the.results = readRDS(paste0("~/Progenitor_Neuronal_population_integration.rds"))

NP_obj = the.results$all_samples
DimPlot(NP_obj, label =T)
#

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
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("26")] <- "V2_pre" # V2 interneuron precursor
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("20")] <- "NP_pro4" # neural progenitor in proliferation
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("3")] <- "NP_pro1" # neural progenitor in proliferation
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("7")] <- "NP_pro3" # neural progenitor in proliferation
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("6")] <- "NP_pro2" # neural progenitor in proliferation
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("27")] <- "MN" # motor neuron
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("28")] <- "MC" # mesenchymal cells
NP_obj$Cell_labels[NP_obj$integrated_snn_res.2 %in% c("31")] <- "KC" # Keratynocytes
##
######################## publication ################################################
# ISN markers:
library(Seurat)
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
       p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
               xlab("") + ylab(feature) + ggtitle("") +
               theme(legend.position = "none",
               axis.text.x = element_blank(),
               #axis.text.x = element_text(angle = 45),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(),
               axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
               plot.margin = margin(0, 0, 0, 0, "cm"),
               plot.title= element_blank(),
               axis.title.x = element_blank(), text = element_text(size = 20))
               #plot.margin = plot.margin )
       return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
       plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
            theme(axis.text.x=element_text(angle = 45, hjust = 0.4, size = 20), axis.ticks.x = element_line())
       p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
       return(p)
}

#StackedVlnPlot(sdata, c('Retnlg', 'Pygl', 'Anxa1', 'Igf1r', 'Stfa2l1'), pt.size=0, cols=my36colors)
p1 <- StackedVlnPlot(NP_obj, c("fev", "lmx1bb", "gata3", "gata2a", "sox1a", "sox1b", "tph2", "slc6a4a", "htr1d"),  pt.size=0, cols=my36colors)
#
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175')
pdf("~/ISN_markers_StackedVlPlot2.pdf", width = 12, height = 6)
print(p1)
dev.off()

### ISN subcluster
ISN_obj <- subset(NP_obj, Cell_labels %in% c("ISN"))

## filtering other neuron cells mixed into ISN at 1dpf
dpf1_obj <- readRDS("~/dpf1_obj.rds")

q1 <- dpf1_obj
q1_apart <- subset(q1, RNA_snn_res.1 %in% c(4, 11, 17, 19))

##
ISN_obj$other_1dpf_lable <- "No"
ISN_obj$other_1dpf_lable[colnames(ISN_obj) %in% colnames(q1_apart)] <- "Yes"
ISN_obj_sel <- subset(ISN_obj, other_1dpf_lable %in% "No")
#
DefaultAssay(ISN_obj_sel) <- "integrated"

#

#
#NP_LFP_ISNpre_ISN_obj_sel2 <- subset(NP_LFP_ISNpre_ISN_obj_sel, Cell_label_Stage %in% names(table(NP_LFP_ISNpre_ISN_obj_sel$Cell_label_Stage))[!c(names(table(NP_LFP_ISNpre_ISN_obj_sel$Cell_label_Stage)) %in% c("LFP_2dpf","LFP_3dpf", "LFP_5dpf"))])

# ISN subcluter 
DefaultAssay(ISN_obj_sel) <- "RNA"
NP_ISNpre_ISN_obj_sel2_list <- SplitObject(ISN_obj_sel, split.by = "orig.ident")

NP_ISNpre_ISN_obj_sel2_list <- lapply(X = NP_ISNpre_ISN_obj_sel2_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = NP_ISNpre_ISN_obj_sel2_list)
anchors <- FindIntegrationAnchors(object.list = NP_ISNpre_ISN_obj_sel2_list, anchor.features = features, k.filter = 5, k.score = 5, dims = 1:10)

# set k.weight so that it is not larger than the number of cells in the smallest dataset
# https://github.com/satijalab/seurat/issues/3936
ISN_combined1 <- IntegrateData(anchorset = anchors, k.weight = 10)
DefaultAssay(ISN_combined1) <- "integrated"

ISN_combined1 <- ScaleData(ISN_combined1,  vars.to.regress = "nCount_RNA", verbose = FALSE)
ISN_combined1 <- RunPCA(ISN_combined1, npcs = 30, verbose = FALSE)
ISN_combined1 <- RunUMAP(ISN_combined1, reduction = "pca", dims = 1:30)
ISN_combined1 <- FindNeighbors(ISN_combined1, reduction = "pca", dims = 1:30)
ISN_combined1 <- FindClusters(ISN_combined1, resolution = 0.3)

DefaultAssay(ISN_combined1) <- "RNA"
DimPlot(ISN_combined1) + DimPlot(ISN_combined1, group.by = "orig.ident") + FeaturePlot(ISN_combined1, c("tph2", "fev", "nkx2.2a"))

FeaturePlot(ISN_combined1, c("fev", "nkx2.2a","lmx1bb","tph2"), split.by = "integrated_snn_res.0.3")

FeaturePlot(ISN_combined1, c("fev",  "gata3", "gata2a", "nkx2.2a", "nkx2.2b", "lmx1bb", "lmx1ba","tph2", "insm1a", "insm1b"))
## plot figure(publication ):
ISN_combined1_umap = ISN_combined1@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  cbind(cell_type = ISN_combined1@meta.data$integrated_snn_res.0.3) 


p <- ggplot(ISN_combined1_umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) + geom_point(size = 1 , alpha =1 ) + scale_color_aaas() #+ scale_color_manual(values = my36colors) #+ scale_color_manual(values = c(brewer.pal(9, "Set1"), "grey"))

#ggplot(Sox1a_umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) + geom_point(size = 1 , alpha =1 ) + scale_color_manual(values = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"))
p2 <- p  +
  theme(panel.grid.major = element_blank(), #
        panel.grid.minor = element_blank(), #
        panel.border = element_blank(), #
        axis.title = element_blank(),  #
        axis.text = element_blank(), # 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #
        plot.background=element_rect(fill="white"))
p2
p3 <- p2 +         
        theme(
          legend.title = element_blank(), # 
          legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=15), #
        legend.key.size=unit(1,'cm') ) +  #
  guides(color = guide_legend(override.aes = list(size=5))) #

p4 <- p3 + geom_segment(aes(x = min(ISN_combined1_umap$UMAP_1)  , y = min(ISN_combined1_umap$UMAP_2) -0.5, xend = min(ISN_combined1_umap$UMAP_1) +1, yend = min(ISN_combined1_umap$UMAP_2) -0.5 ), colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))  + geom_segment(aes(x = min(ISN_combined1_umap$UMAP_1)  , y = min(ISN_combined1_umap$UMAP_2) - 0.5, xend = min(ISN_combined1_umap$UMAP_1) , yend = min(ISN_combined1_umap$UMAP_2) + 0.5), colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(ISN_combined1_umap$UMAP_1) +0.5, y = min(ISN_combined1_umap$UMAP_2) -0.8, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(ISN_combined1_umap$UMAP_1) -0.4, y = min(ISN_combined1_umap$UMAP_2) + 0.03, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)

p4

pdf("~/ISN_combined1_SubClusters.pdf", width = 10)
print(p4)
dev.off()
png("~/ISN_combined1_SubClusters.png", width = 480 * 2, height = 480 * 1.5 )
print(p4)
dev.off()
##
library(ggeasy)
pdf("~/ISN_combined1_SubClusters2.pdf", width = 15, height = 10)
print(DimPlot(ISN_combined1, label =F, group.by = "integrated_snn_res.0.3", pt.size = 6 )+ easy_remove_axes() + theme(text = element_text(size =25), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=30)) +  ggtitle("") + scale_color_aaas() + guides(colour = guide_legend(override.aes = list(size=5))))
dev.off()
#

##
# group by stage
## stage
library(ggeasy)
pdf("~/ISN_combined1_SubClusters_DimPlot_Stage2.pdf", width = 15, height = 10)
print(DimPlot(ISN_combined1, label =F, group.by = "orig.ident", pt.size = 6, cols = c("#999999", "#333333", "#ff0000", "#00cdcd") )+ easy_remove_axes() + theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=25)) +  ggtitle("") + guides(colour = guide_legend(override.aes = list(size=5))))
dev.off()

png("~/ISN_combined1_SubClusters_DimPlot_Stage.png", width = 480 * 3, height = 480 * 2)
print(DimPlot(ISN_combined1, label =F, group.by = "orig.ident", pt.size = 3, cols = c("#999999", "#333333", "#ff0000", "#00cdcd") )+ easy_remove_axes() + theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=20)) +  ggtitle(""))
dev.off()
# label ISN_pre and ISN
ISN_combined1$ISN_labels <- "ISN"
ISN_combined1$ISN_labels[ISN_combined1$integrated_snn_res.0.3 == 1] <- "ISN_pre"

ISN_subcluster_plot <- DimPlot(ISN_combined1, label =F, group.by = "ISN_labels", pt.size = 10 )+ easy_remove_axes() + scale_color_manual(values = c( "ISN_pre" = "#91D0BE", "ISN" = '#476D87'))  + theme(text = element_text(size = 50), legend.text = element_text(size=50)) +  ggtitle("") + guides(colour = guide_legend(override.aes = list(size=10)))

##
pdf("~/ISNsubcluster_ISN_pre_DimPlot.pdf", width = 15, height = 10)
#print(GFP_plot + Sox1a_plot)
DimPlot(ISN_combined1, label =F, group.by = "ISN_labels", pt.size = 6 )+ easy_remove_axes() + scale_color_manual(values = c( "ISN_pre" = "#91D0BE", "ISN" = '#476D87'))  + theme(text = element_text(size = 20), legend.text = element_text(size=25)) +  ggtitle("") + guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()
###


pdf("~/ISNsubcluster_FeaturePlot2.pdf", width = 20* 3, height = 15*1)
#print(GFP_plot + Sox1a_plot)
grid.arrange(ISN_subcluster_plot, nkx2.2a_plot, nkx2.2b_plot, fev_plot, lmx1bb_plot, gata2a_plot, gata3_plot, tph2_plot, ncol=4, nrow=2)

dev.off()

pdf("~/ISNsubcluster_FeaturePlot4.pdf", width = 20* 3, height = 15*1.5)
#print(GFP_plot + Sox1a_plot)
grid.arrange(nkx2.2a_plot, fev_plot, gata2a_plot, tph2_plot, nkx2.2b_plot, lmx1bb_plot, gata3_plot,  ddc_plot, ncol=4, nrow=2)

dev.off()


# Feature
library(gridExtra)
library(ggeasy)

c("nkx2.2a", "nkx2.2b", "fev", "lmx1ba", "lmx1bb", "gata2a", "gata3", "tph2")

nkx2.2a_plot <- FeaturePlot(ISN_combined1, "nkx2.2a", pt.size = 10) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
nkx2.2b_plot <- FeaturePlot(ISN_combined1, "nkx2.2b", pt.size = 10) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
fev_plot <- FeaturePlot(ISN_combined1, "fev", pt.size = 10) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
lmx1ba_plot <- FeaturePlot(ISN_combined1, "lmx1ba", pt.size = 10) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
lmx1bb_plot <- FeaturePlot(ISN_combined1, "lmx1bb", pt.size = 10) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
gata2a_plot <- FeaturePlot(ISN_combined1, "gata2a", pt.size = 10) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
gata3_plot <- FeaturePlot(ISN_combined1, "gata3", pt.size = 10) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
tph2_plot <- FeaturePlot(ISN_combined1, "tph2", pt.size = 10) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

ddc_plot <- FeaturePlot(ISN_combined1, "ddc", pt.size = 10) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

pdf("~/ISNsubcluster_FeaturePlot.pdf", width = 20* 3, height = 15*1)
#print(GFP_plot + Sox1a_plot)
grid.arrange(nkx2.2a_plot, nkx2.2b_plot, fev_plot, lmx1ba_plot, lmx1bb_plot, gata2a_plot, gata3_plot, tph2_plot, ncol=4, nrow=2)

dev.off()


png("~/ISNsubcluster_FeaturePlot.png", width = 800 * 3.5, height = 800*1.5)
grid.arrange(nkx2.2a_plot, nkx2.2b_plot, fev_plot, lmx1ba_plot, lmx1bb_plot, gata2a_plot, gata3_plot, tph2_plot, ncol=4, nrow=2)

dev.off()
##
pdf("~/ISNsubcluster_FeaturePlot_Test.pdf", width = 15*2.5 , height = 20* 2)
#print(GFP_plot + Sox1a_plot)
grid.arrange(nkx2.2a_plot, nkx2.2b_plot, fev_plot, lmx1bb_plot, gata2a_plot, gata3_plot, tph2_plot, ddc_plot,ncol=2, nrow=4)

dev.off()



## find LFP, ISN, ISN_pre markers and TF: c(this is LFP compare with NP all other cells)

NP_obj$ISN_ISN_Pre_labels <- "Other"
NP_obj$ISN_ISN_Pre_labels[colnames(subset(ISN_combined1, integrated_snn_res.0.3 %in% 1))] <- "ISN_pre"
NP_obj$ISN_ISN_Pre_labels[colnames(subset(ISN_combined1, integrated_snn_res.0.3 %in% c(0, 2)))] <- "ISN"
NP_obj$ISN_ISN_Pre_labels[colnames(subset(NP_obj, Cell_labels %in% c("LFP")))] <- "LFP"
saveRDS(NP_obj, file = "~/NP_obj_with_ISNpre_labels.rds")
# LFP markers
LFP_markers <- FindMarkers(NP_obj, ident.1 = "LFP", group.by = "ISN_ISN_Pre_labels", only.pos = TRUE)
LFP_markers$gene <- rownames(LFP_markers)
# ISN_pre markers
ISN_pre_markers <- FindMarkers(NP_obj, ident.1 = "ISN_pre", group.by = "ISN_ISN_Pre_labels", only.pos = TRUE)
ISN_pre_markers$gene <- rownames(ISN_pre_markers)
# ISN markers
ISN_markers <- FindMarkers(NP_obj, ident.1 = "ISN", group.by = "ISN_ISN_Pre_labels", only.pos = TRUE)
ISN_markers$gene <- rownames(ISN_markers)

LFP_ISN_ISN_pre_markers_list <- list()
LFP_ISN_ISN_pre_markers_list[["LFP"]] <- LFP_markers
LFP_ISN_ISN_pre_markers_list[["ISN_pre"]] <- ISN_pre_markers
LFP_ISN_ISN_pre_markers_list[["ISN"]] <- ISN_markers

writexl::write_xlsx(LFP_ISN_ISN_pre_markers_list, "~/LFP_ISN_ISN_pre_obj.all.mrkrs.xlsx")
##
# add Ensembl
features <- read.table("~/features.tsv.gz", header = FALSE, sep = "\t")

getUniqs = function(xx){
  idx2 = duplicated(xx)
  for(i in seq(1, length(idx2),1)) if(idx2[i]) xx[i] = paste0(xx[i],"#",i)
  return(xx)
}

features$gene <- getUniqs(features$V2)
colnames(features) <- c("Ensembl", "V2", "V3", "gene")

LFP_markers_ensembl <- merge(LFP_markers, features[, c("Ensembl", "gene")], by = "gene", sort = FALSE)
ISN_pre_markers_ensembl <- merge(ISN_pre_markers, features[, c("Ensembl", "gene")], by = "gene", sort = FALSE)
ISN_markers_ensembl <- merge(ISN_markers, features[, c("Ensembl", "gene")], by = "gene", sort = FALSE)

LFP_ISN_ISN_pre_markers_list <- list()
LFP_ISN_ISN_pre_markers_list[["LFP"]] <- LFP_markers_ensembl
LFP_ISN_ISN_pre_markers_list[["ISN_pre"]] <- ISN_pre_markers_ensembl
LFP_ISN_ISN_pre_markers_list[["ISN"]] <- ISN_markers_ensembl

writexl::write_xlsx(LFP_ISN_ISN_pre_markers_list, "~/LFP_ISN_ISN_pre_obj.all.mrkrs.xlsx")
#


# select by p_val_adj < 0.05
LFP_markers_sel <- LFP_markers[LFP_markers$p_val_adj < 0.05, ]
ISN_pre_markers_sel <- ISN_pre_markers[ISN_pre_markers$p_val_adj < 0.05, ]
ISN_markers_sel <- ISN_markers[ISN_markers$p_val_adj < 0.05, ]

LFP_ISN_ISN_pre_markers_sel_list <- list()
LFP_ISN_ISN_pre_markers_sel_list[["LFP"]] <- LFP_markers_sel
LFP_ISN_ISN_pre_markers_sel_list[["ISN_pre"]] <- ISN_pre_markers_sel
LFP_ISN_ISN_pre_markers_sel_list[["ISN"]] <- ISN_markers_sel

writexl::write_xlsx(LFP_ISN_ISN_pre_markers_sel_list, "~/LFP_ISN_ISN_pre_obj.pval0.05.mrkrs.xlsx")

LFP_ISN_Pre_overlap_markers <- LFP_markers_sel[LFP_markers_sel$gene %in% ISN_pre_markers_sel$gene, ]$gene


##select TF in LFP, ISN_pre and ISN
TF_DB <- read.table("~/AnimalTFDB/Danio_rerio_TF.txt", header = T, sep = "\t")
LFP_TF_markers <- LFP_markers[LFP_markers$gene %in% TF_DB$Symbol, ]
ISN_pre_TF_markers <- ISN_pre_markers[ISN_pre_markers$gene %in% TF_DB$Symbol, ]
ISN_TF_markers <- ISN_markers[ISN_markers$gene %in% TF_DB$Symbol, ]
# select
LFP_TF_markers_sel <- LFP_markers_sel[LFP_markers_sel$gene %in% TF_DB$Symbol, ]
ISN_pre_TF_markers_sel <- ISN_pre_markers_sel[ISN_pre_markers_sel$gene %in% TF_DB$Symbol, ]
ISN_TF_markers_sel <- ISN_markers_sel[ISN_markers_sel$gene %in% TF_DB$Symbol, ]

LFP_ISN_ISN_pre_TF_markers_sel_list <- list()
LFP_ISN_ISN_pre_TF_markers_sel_list[["LFP"]] <- LFP_TF_markers_sel
LFP_ISN_ISN_pre_TF_markers_sel_list[["ISN_pre"]] <- ISN_pre_TF_markers_sel
LFP_ISN_ISN_pre_TF_markers_sel_list[["ISN"]] <- ISN_TF_markers_sel

writexl::write_xlsx(LFP_ISN_ISN_pre_TF_markers_sel_list, "~/LFP_ISN_ISN_pre_obj.pval0.05.mrkrs_TF.xlsx")


## Plot markers(gene or TF) 
pdf("~/LFP_ISN_ISN_pre_feature_dotplot.pdf", width = 15, height = 4)
DotPlot(subset(NP_obj, ISN_ISN_Pre_labels %in% c("LFP", "ISN_pre", "ISN")), group.by = "ISN_ISN_Pre_labels", features= c("nkx2.9", "nkx6.2", "sulf1","foxa","nkx2.2a", "nkx2.2b", "sox3", "sox11b", "sox5","fev", "lmx1ba", "lmx1bb", "shox2", "gata2a","tph2", "gata3", "slc7a8a"), scale = F) + scale_color_viridis_c(direction = 1) + cowplot::theme_cowplot() + theme(axis.line  = element_blank())  + ylab('') + theme(axis.ticks = element_blank()) + theme(axis.text.x = element_text(angle = 45,  hjust=1))+  theme(legend.key.height = unit(0.4, 'cm'))
dev.off()

png("~/LFP_ISN_ISN_pre_feature_dotplot.png", width = 480 * 2, height = 250)
DotPlot(subset(NP_obj, ISN_ISN_Pre_labels %in% c("LFP", "ISN_pre", "ISN")), group.by = "ISN_ISN_Pre_labels", features= c("nkx2.9", "nkx6.2", "sulf1","foxa","nkx2.2a", "nkx2.2b", "sox3", "sox11b", "sox5","fev", "lmx1ba", "lmx1bb", "shox2", "gata2a","tph2", "gata3", "slc7a8a"), scale = F) + scale_color_viridis_c(direction = 1) + cowplot::theme_cowplot() + theme(axis.line  = element_blank())  + ylab('') + theme(axis.ticks = element_blank()) + theme(axis.text.x = element_text(angle = 45,  hjust=1))+  theme(legend.key.height = unit(0.4, 'cm'))
dev.off()
# ISN subcluster findMarkers:

### ISN to find markers
Idents(ISN_combined1) = ISN_combined1$integrated_snn_res.0.3 

ISN_combined1_subcluster.all.mrkrs = FindAllMarkers(ISN_combined1, only.pos = T)
saveRDS(ISN_combined1_subcluster.all.mrkrs, file = "~/ISN_combined1_subcluster.all.mrkrs_new.rds")

ISN_combined1_subcluster_obj.all.mrkrs.list <- list()
#for(i in sort(as.character(unique(ISN_combined1_subcluster_obj.all.mrkrs$cluster)))){
for(i in as.character(unique(sort(ISN_combined1$integrated_snn_res.0.3)))){
 
    ISN_combined1_subcluster_obj.all.mrkrs.list[[i]] <- ISN_combined1_subcluster.all.mrkrs[ISN_combined1_subcluster.all.mrkrs$cluster == i, ]

}


writexl::write_xlsx(ISN_combined1_subcluster_obj.all.mrkrs.list, "~/ISN_combined1_subcluster_obj.all.mrkrs.xlsx")
#
# add Ensembl
ISN_combined1_subcluster.all.mrkrs <- readRDS("~/ISN_combined1_subcluster.all.mrkrs_new.rds")
features <- read.table("~/features.tsv.gz", header = FALSE, sep = "\t")

getUniqs = function(xx){
  idx2 = duplicated(xx)
  for(i in seq(1, length(idx2),1)) if(idx2[i]) xx[i] = paste0(xx[i],"#",i)
  return(xx)
}

features$gene <- getUniqs(features$V2)
colnames(features) <- c("Ensembl", "V2", "V3", "gene")
#
ISN_combined1_subcluster_obj.all.mrkrs.list <- list()

for(i in as.character(unique(sort(ISN_combined1_subcluster.all.mrkrs$cluster)))){

    One_Cell_ISNsubcluster_obj.all.mrkrs <- ISN_combined1_subcluster.all.mrkrs[ISN_combined1_subcluster.all.mrkrs$cluster == i, ]
    One_Cell_ISNsubcluster_obj.all.mrkrs.Ensembl <- merge(One_Cell_ISNsubcluster_obj.all.mrkrs, features[, c("Ensembl", "gene")], by = "gene", sort = FALSE)
    ISN_combined1_subcluster_obj.all.mrkrs.list[[i]] <- One_Cell_ISNsubcluster_obj.all.mrkrs.Ensembl
    
}

writexl::write_xlsx(ISN_combined1_subcluster_obj.all.mrkrs.list, "~/ISN_combined1_subcluster_obj.all.mrkrs.xlsx")


###
### ISN to compare markers for ISN_pre and ISN
ISN_combined1$ISN_labels <- "ISN"
ISN_combined1$ISN_labels[ISN_combined1$integrated_snn_res.0.3 == 1] <- "ISN_pre"

Idents(ISN_combined1) = ISN_combined1$ISN_labels 

saveRDS(ISN_combined1, file = "~/ISN_combined1_subcluster.rds")


ISN_combined1_subcluster.all.mrkrs = FindAllMarkers(ISN_combined1, only.pos = T)

ISN_combined1_subcluster_obj.all.mrkrs.list <- list()
#for(i in sort(as.character(unique(ISN_combined1_subcluster_obj.all.mrkrs$cluster)))){
for(i in as.character(unique(sort(ISN_combined1$ISN_labels)))){
 
    ISN_combined1_subcluster_obj.all.mrkrs.list[[i]] <- ISN_combined1_subcluster.all.mrkrs[ISN_combined1_subcluster.all.mrkrs$cluster == i, ]

}


writexl::write_xlsx(ISN_combined1_subcluster_obj.all.mrkrs.list, "~/ISN_combined1_subcluster_CellLabels_obj.all.mrkrs.xlsx")

##################################################
## LFP, ISN
ISN_obj_sel$ISN_labels <- ISN_combined1$ISN_labels
DimPlot(ISN_combined1, group.by = "ISN_labels")
DimPlot(ISN_obj_sel, group.by = "ISN_labels")



##
NP_LFP_ISNpre_ISN_obj <- subset(NP_obj, Cell_labels %in% c("LFP", "ISN"))
NP_LFP_ISNpre_ISN_obj$Cell_label_Stage <- paste(NP_LFP_ISNpre_ISN_obj$Cell_labels, NP_LFP_ISNpre_ISN_obj$orig.ident, sep = "_")

##
dpf1_obj <- readRDS("~/dpf1_obj.rds")

q1 <- dpf1_obj
q1_apart <- subset(q1, RNA_snn_res.1 %in% c(4, 11, 17, 19))

##
NP_LFP_ISNpre_ISN_obj$other_1dpf_lable <- "No"
NP_LFP_ISNpre_ISN_obj$other_1dpf_lable[colnames(NP_LFP_ISNpre_ISN_obj) %in% colnames(q1_apart)] <- "Yes"
NP_LFP_ISNpre_ISN_obj_sel <- subset(NP_LFP_ISNpre_ISN_obj, other_1dpf_lable %in% "No")

#
NP_LFP_ISNpre_ISN_obj_sel2 <- subset(NP_LFP_ISNpre_ISN_obj_sel, Cell_label_Stage %in% names(table(NP_LFP_ISNpre_ISN_obj_sel$Cell_label_Stage))[!c(names(table(NP_LFP_ISNpre_ISN_obj_sel$Cell_label_Stage)) %in% c("LFP_2dpf","LFP_3dpf", "LFP_5dpf"))])

#
## LFP , filtering ISN_Pre, ISN
NP_LFP_ISNpre_ISN_obj <- readRDS("~/NP_LFP_ISN_pre_ISN_filterlabel_LFP_ISN.rds")
NP_LFP_ISNpre_ISN_obj_sel <- subset(NP_LFP_ISNpre_ISN_obj, other_1dpf_lable %in% "No")

obj.combined2 <- readRDS(file = "~/ISN_combined1_subcluster.rds")

#  names 
NP_LFP_ISNpre_ISN_obj_sel$Cell_labels_New <- "1"
NP_LFP_ISNpre_ISN_obj_sel$Cell_labels_New[colnames(NP_LFP_ISNpre_ISN_obj_sel) %in% colnames(subset(obj.combined2, ISN_labels %in% "ISN"))] <- "ISN"
NP_LFP_ISNpre_ISN_obj_sel$Cell_labels_New[colnames(NP_LFP_ISNpre_ISN_obj_sel) %in% colnames(subset(obj.combined2, ISN_labels %in% "ISN_pre"))] <- "ISN_pre"
NP_LFP_ISNpre_ISN_obj_sel$Cell_labels_New[colnames(NP_LFP_ISNpre_ISN_obj_sel) %in% colnames(subset(NP_LFP_ISNpre_ISN_obj_sel, Cell_labels %in% "LFP"))] <- "LFP"

# Plot LFP, ISN_pre, ISN
pdf("~/LFP_ISN_DimPlot_split_Stage.pdf", width = 15 * 2.5, height = 12)
DimPlot(NP_LFP_ISNpre_ISN_obj_sel, label =F, group.by = "Cell_labels_New", pt.size = 2 , split.by = "orig.ident") + scale_color_manual(values = c("LFP" = '#57C3F3', "ISN_pre" = "#91D0BE", "ISN" = '#476D87')) + theme(text = element_text(size = 30), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=30)) + ggtitle("")

dev.off()


###
nkx2.9_plot <- FeaturePlot(NP_LFP_ISNpre_ISN_obj_sel, c("nkx2.9"), pt.size = 2 , split.by = "orig.ident") 
nkx2.2a_plot <- FeaturePlot(NP_LFP_ISNpre_ISN_obj_sel, c("nkx2.2a"), pt.size = 2 , split.by = "orig.ident") 
fev_plot <- FeaturePlot(NP_LFP_ISNpre_ISN_obj_sel, c("fev"), pt.size = 2 , split.by = "orig.ident") 
lmx1bb_plot <- FeaturePlot(NP_LFP_ISNpre_ISN_obj_sel, c("lmx1bb"), pt.size = 2 , split.by = "orig.ident") 
tph2_plot <- FeaturePlot(NP_LFP_ISNpre_ISN_obj_sel, c("tph2"), pt.size = 2 , split.by = "orig.ident") 

#
ISN_FeaturePlot <- FeaturePlot(NP_LFP_ISNpre_ISN_obj_sel, c("nkx2.9", "nkx2.2a", "fev", "lmx1bb", "tph2"), pt.size = 2 , split.by = "orig.ident") 

pdf("~/ISN_nkx2.9_2.2a_fev_lmx1bb_tph2_FeaturePlot.pdf", width = 20* 2, height = 15*2)

print(ISN_FeaturePlot)
dev.off()
##
### Cell Label to find markers: Comparing markers between LFP. ISN-pre and ISN
Idents(NP_LFP_ISNpre_ISN_obj_sel) = NP_LFP_ISNpre_ISN_obj_sel$Cell_labels_New 

LFP_ISN_obj.all.mrkrs = FindAllMarkers(NP_LFP_ISNpre_ISN_obj_sel, only.pos = T)
LFP_ISN_obj.all.mrkrs_sel <- LFP_ISN_obj.all.mrkrs[LFP_ISN_obj.all.mrkrs$p_val_adj < 0.05, ]
LFP_ISN_obj.all.mrkrs.list <- list()
for(i in as.character(unique(LFP_ISN_obj.all.mrkrs$cluster))){
 
    LFP_ISN_obj.all.mrkrs.list[[i]] <- LFP_ISN_obj.all.mrkrs_sel[LFP_ISN_obj.all.mrkrs_sel$cluster == i, ]

}


writexl::write_xlsx(LFP_ISN_obj.all.mrkrs.list, "~/LFP_ISN_Cell_label_Cluster_markers_new.xlsx")
## LFP, ISN_pre vs ISN
LFP_vs_ISN_markers <- FindMarkers(object = NP_LFP_ISNpre_ISN_obj_sel, ident.1 = "LFP", ident.2 = "ISN", only.pos = T)
LFP_vs_ISN_markers_sel <- LFP_vs_ISN_markers[LFP_vs_ISN_markers$p_val_adj < 0.05, ]
LFP_vs_ISN_markers_sel2 <- LFP_vs_ISN_markers_sel[!c(substr(rownames(LFP_vs_ISN_markers_sel), 1, 2) %in% c("hb", "rp")), ]

ISN_pre_vs_ISN_markers <- FindMarkers(object = NP_LFP_ISNpre_ISN_obj_sel, ident.1 = "ISN_pre", ident.2 = "ISN", only.pos = T)
ISN_pre_vs_ISN_markers_sel <- ISN_pre_vs_ISN_markers[ISN_pre_vs_ISN_markers$p_val_adj < 0.05, ]
ISN_pre_vs_ISN_markers_sel2 <- ISN_pre_vs_ISN_markers_sel[!c(substr(rownames(ISN_pre_vs_ISN_markers_sel), 1, 2) %in% c("hb", "rp")), ]
length(rownames(LFP_vs_ISN_markers_sel2[rownames(LFP_vs_ISN_markers_sel2) %in% rownames(ISN_pre_vs_ISN_markers_sel2), ]))
# ISN ISN_pre vs LFP
ISN_vs_LFP_markers <- FindMarkers(object = NP_LFP_ISNpre_ISN_obj_sel, ident.1 = "ISN", ident.2 = "LFP", only.pos = T)
ISN_vs_LFP_markers_sel <- ISN_vs_LFP_markers[ISN_vs_LFP_markers$p_val_adj < 0.05, ]
ISN_vs_LFP_markers_sel2 <- ISN_vs_LFP_markers_sel[!c(substr(rownames(ISN_vs_LFP_markers_sel), 1, 2) %in% c("hb", "rp")), ]

ISN_pre_vs_LFP_markers <- FindMarkers(object = NP_LFP_ISNpre_ISN_obj_sel, ident.1 = "ISN_pre", ident.2 = "LFP", only.pos = T)
ISN_pre_vs_LFP_markers_sel <- ISN_pre_vs_LFP_markers[ISN_pre_vs_LFP_markers$p_val_adj < 0.05, ]
ISN_pre_vs_LFP_markers_sel2 <- ISN_pre_vs_LFP_markers_sel[!c(substr(rownames(ISN_pre_vs_LFP_markers_sel), 1, 2) %in% c("hb", "rp")), ]

length(rownames(ISN_vs_LFP_markers_sel2[rownames(ISN_vs_LFP_markers_sel2) %in% rownames(ISN_pre_vs_LFP_markers_sel2), ]))

####################################################################################################################
######## Trajectory analysis by dataset(filtering LFP, and ISN_pre ISN)

NP_LFP_ISNpre_ISN_obj_sel2 <- readRDS("~/NP_LFP_ISN_pre_ISN_filter_LFP_ISN.rds")

NP_obj_LFP_ISN <- NP_LFP_ISNpre_ISN_obj_sel2
# read New ISN subcluster
obj.combined2 <- readRDS(file = "~/ISN_combined1_subcluster.rds")

#
NP_obj_LFP_ISN$Cell_labels_New <- "1"
NP_obj_LFP_ISN$Cell_labels_New[colnames(NP_obj_LFP_ISN) %in% colnames(subset(obj.combined2, ISN_labels %in% "ISN"))] <- "ISN"
NP_obj_LFP_ISN$Cell_labels_New[colnames(NP_obj_LFP_ISN) %in% colnames(subset(obj.combined2, ISN_labels %in% "ISN_pre"))] <- "ISN_pre"
NP_obj_LFP_ISN$Cell_labels_New[colnames(NP_obj_LFP_ISN) %in% colnames(subset(NP_obj_LFP_ISN, Cell_labels %in% "LFP"))] <- "LFP"

NP_obj_LFP_ISN$Cell_labels <- NULL
saveRDS(NP_obj_LFP_ISN, "~/NP_filtered_LFP_ISNpre_ISN.rds")
## cell number in these 
# the bar plot for cell numbers at every stage:
library(ggplot2)

df_cell_num <- data.frame(table(NP_obj_LFP_ISN$orig.ident))
colnames(df_cell_num) <- c("Stages", "Num")

p_bar_cell_num <- ggplot(df_cell_num, aes(x=Stages, y=Num, fill=Stages))+ geom_bar(stat="identity")+ scale_fill_manual(values=c("#999999", "#333333", "#ff0000", "#00cdcd")) + geom_text(aes(label=Num), vjust=-0.25, position = position_dodge(0.9), size=5) + theme_minimal() + xlab("") + ylab("Cell Numbers") + theme(legend.position = "none")
pdf("~/LFP_ISNpre_ISN_Cells_Number_bar.pdf")
print(p_bar_cell_num)
dev.off()
##
## cell number and cell types
df_cell_number_cell_types <- data.frame(table( NP_obj_LFP_ISN$orig.ident, NP_obj_LFP_ISN$Cell_labels_New))
colnames(df_cell_number_cell_types) <- c("Stages", "Cell_labels", "Num")
df_cell_number_cell_types$Cell_labels <- factor(df_cell_number_cell_types$Cell_labels,  levels = c( "LFP",  "ISN_pre", "ISN"))

# contribution of each stage to the cell types identified
p_bar_cellPercent_cellTypes2 <- ggplot(data=df_cell_number_cell_types, aes(x=Cell_labels, y=Num, fill=Stages)) +
  geom_bar(stat="identity", position="fill")+ scale_fill_manual(values=c("#999999", "#333333", "#ff0000", "#00cdcd")) + theme_minimal() + xlab("") + ylab("Cell Percentage")
pdf("~/LFP_ISNpre_ISN_CellsStagesPercent_CellTypes_bar.pdf")
print(p_bar_cellPercent_cellTypes2)
dev.off()
#
df_cell_number_cell_types2 <- data.frame(table( NP_obj_LFP_ISN$Cell_labels_New))
colnames(df_cell_number_cell_types2) <- c( "Cell_labels", "Num")
df_cell_number_cell_types2$Cell_labels <- factor(df_cell_number_cell_types2$Cell_labels,  levels = c( "LFP",  "ISN_pre", "ISN"))

p_bar_cellNum_cellTypes <- ggplot(data=df_cell_number_cell_types2, aes(x=Cell_labels, y=Num, fill=Cell_labels)) +
  geom_bar(stat="identity")+ scale_fill_manual(values = c("LFP" = '#57C3F3', "ISN_pre" = "#91D0BE", "ISN" = '#476D87')) + geom_text(aes(label=Num), vjust=-0.25, position = position_dodge(0.9), size=5) + theme_minimal() + xlab("") + ylab("Cell Numbers") + theme(legend.position = "none")
pdf("~/LFP_ISNpre_ISN_CellsNum_CellTypes_bar.pdf")
print(p_bar_cellNum_cellTypes)
dev.off()

## Monocle2
LFP_ISN_data <- as(as.matrix(NP_obj_LFP_ISN@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = NP_obj_LFP_ISN@meta.data)
fData <- data.frame(gene_short_name = row.names(LFP_ISN_data), row.names = row.names(LFP_ISN_data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(LFP_ISN_data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

#monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

print(head(fData(monocle_cds))) # 22804

expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= 10)) # 10855 
print(head(pData(monocle_cds)))

## clustering cells without marker genes
## ##monocle
#disp_table <- dispersionTable(monocle_cds)
#unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
#unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit )

#monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
#plot_ordering_genes(monocle_cds)
## bulid trajectory 
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes, ], fullModelFormulaStr = "~Cell_labels_New")

#diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes, ], fullModelFormulaStr = "~Cell_labels_New + orig.ident", reducedModelFormulaStr = " ~orig.ident", relative_expr=TRUE)

df_ordering_genes <- subset(diff_test_res, qval < 0.01) # qval < 0.1, qval < 0.01
df_ordering_genes_order <- df_ordering_genes[order(df_ordering_genes$qval, decreasing = F), ]
writexl::write_xlsx(df_ordering_genes_order, "~/LFP_ISN_monocle2_trajectory_DEGs_qval0.01.xlsx")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01)) # qval < 0.1, qval < 0.01
saveRDS(ordering_genes, "~/LFP_ISN_monocle2_trajectory_DEGs_qval0.01_ordering_genes.rds" )
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)

plot_ordering_genes(monocle_cds)

# reduce data dimensionality

monocle_cds <- reduceDimension(monocle_cds, max_components = 2, method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds, reverse = F)

monocle_cds$Cell_labels_New <- factor(monocle_cds$Cell_labels_New, levels = c("LFP", "ISN_pre", "ISN"))

saveRDS(monocle_cds, "~/monocle_cds.rds")

pdf("~/LFP_ISN_Trajectory.pdf", width = 10)
plot_cell_trajectory(monocle_cds, color_by = "Cell_labels_New")
dev.off()

pdf("~/LFP_ISN_Trajectory_stage.pdf", width = 10)
plot_cell_trajectory(monocle_cds, color_by = "orig.ident")
dev.off()
#
#monocle_cds <- reduceDimension(monocle_cds, max_components = 2, reduction_method = 'DDRTree', sigma = 0.25)
pdf("~/LFP_ISN_Pseudotime.pdf", width = 10)

plot(plot_cell_trajectory(monocle_cds, show_cell_names = F, color_by = "Pseudotime") + scale_color_viridis_c())
dev.off()
plot_genes_in_pseudotime(monocle_cds[c("nkx2.9", "foxa", "nkx2.2a", "fev", "tph2")], color_by = "Cell_labels_New", label_by_short_name = c("nkx2.9", "foxa", "nkx2.2a", "fev", "tph2"))
#plot_genes_in_pseudotime(monocle_cds[c("dld", "dlb", "foxa")], color_by = "Cell_labels_New")

plot_genes_jitter(monocle_cds[c("nkx2.9", "nkx2.2a", "fev", "tph2")], grouping = "Cell_labels_New",color_by = "Cell_labels_New")
plot_genes_violin(monocle_cds[c("nkx2.9", "nkx2.2a", "fev", "tph2")], grouping = "Cell_labels_New",color_by = "Cell_labels_New")
## plot gene in pseudotime
colnames(pData(monocle_cds))

pData(monocle_cds)$fev <- log2(exprs(monocle_cds)['fev', ] + 1)
#plot_cell_trajectory(monocle_cds, color_by="fev")
pData(monocle_cds)$tph2 <- log2(exprs(monocle_cds)['tph2', ] + 1)
#plot_cell_trajectory(monocle_cds, color_by="tph2")+ scale_color_gradient2()
pData(monocle_cds)$nkx2.2a <- log2(exprs(monocle_cds)['nkx2.2a', ] + 1)
pData(monocle_cds)$nkx2.9 <- log2(exprs(monocle_cds)['nkx2.9', ] + 1)

#plot_cell_trajectory(monocle_cds, color_by="nkx2.9")+ scale_color_gradient2(name = NULL, high = "#FF0000", mid = "#FFFFFF", low = "#0000FF")
p1 <- plot_cell_trajectory(monocle_cds, color_by="fev")+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF")
p2 <- plot_cell_trajectory(monocle_cds, color_by="tph2")+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF")
p3 <- plot_cell_trajectory(monocle_cds, color_by="nkx2.2a")+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF")
p4 <- plot_cell_trajectory(monocle_cds, color_by="nkx2.9")+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF")

p4 + p3 + p1 + p2
#library(ggsci)

#plot_cell_trajectory(monocle_cds, color_by="nkx2.9")+ scale_color_gsea()

## Find Genes that Change as a Function of Pseudotime
Time_diff <- differentialGeneTest(monocle_cds[ordering_genes, ], cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_diff <- Time_diff[, c(5, 2, 3, 4, 1, 6, 7)]

writexl::write_xlsx(Time_diff, "~/LFP_ISN_monocle2_trajectory_DEGs_qval0.01_todo_Time_diff.xlsx")

library(tidyverse)

Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

pdf("~/Time_genes_heatmap.pdf", height = 8)

Time_heatmap_plot <- plot_pseudotime_heatmap(monocle_cds[Time_genes, ],  num_clusters = 4, show_rownames = F, return_heatmap = T)

dev.off()
#
df_Time_sig_gene <- subset(Time_diff, qval < 0.00001)

Time_sig_gene_names <- subset(Time_diff, qval < 0.00001)$gene_short_name

pdf("~/Time_sig_genes_heatmap.pdf", height = 8)

p_heatmap <- plot_pseudotime_heatmap(monocle_cds[Time_sig_gene_names, ],   show_rownames = F,  num_clusters = 4, return_heatmap = T)
dev.off()
plot_pseudotime_heatmap(monocle_cds[c("nkx2.9", "nkx2.2a", "foxa", "fev",  "lmx1bb", "tph2"), ],   show_rownames = TRUE,  num_clusters = 4, return_heatmap = T, cluster_rows = F)
# getting the heatmap cluster genes
clusters <- cutree(p_heatmap$tree_row, k = 4)
df_clusters <- data.frame(clusters)
df_clusters$gene_names <- rownames(df_clusters)

write.csv(df_clusters,  "~/LFP_ISN_monocle2_trajectory_DEGs_qval0.00001_heatmap_clusters.csv", row.names = F)
Time_sig_gene_names_clusters <- merge(df_Time_sig_gene, df_clusters, by.x = "gene_short_name", by.y = "gene_names", sort = F)
Time_sig_gene_names_clusters_order <- Time_sig_gene_names_clusters[order(Time_sig_gene_names_clusters$qval), ]

Time_sig_gene_names_clusters_order_list <- list()
for(i in 1:4){
    Time_sig_gene_names_clusters_order_list[[paste("Cluster", i, sep = "")]]<- Time_sig_gene_names_clusters_order[Time_sig_gene_names_clusters_order$clusters == i, ]

}

writexl::write_xlsx(Time_sig_gene_names_clusters_order_list, "~/LFP_ISN_monocle2_trajectory_Sig_gene_Heatmap_Gene_clusters.xlsx")

# cluster3
pdf("~/Time_sig_genes_heatmap_clusterN3.pdf", height = 8)

p_heatmap2 <- plot_pseudotime_heatmap(monocle_cds[Time_sig_gene_names, ],   show_rownames = F,  num_clusters = 3, return_heatmap = T)
dev.off()
###
#Time_sig_gene_names <- row.names(subset(Time_diff, qval < 0.00000001))
#p_heatmap3 <- plot_pseudotime_heatmap(monocle_cds[Time_sig_gene_names, ],   show_rownames = F,  num_clusters = 4, return_heatmap = T)


plot_pseudotime_heatmap(monocle_cds[c("nkx2.9", "nkx2.2a", "foxa", "fev",  "lmx1bb", "tph2"), ],   show_rownames = TRUE,  num_clusters = 4, return_heatmap = T, cluster_rows = F, norm_method = c( "vstExprs"))


#pData(monocle_cds)$dcn <- log2(exprs(monocle_cds)['dcn', ] + 1)
#plot_cell_trajectory(monocle_cds, color_by="dcn")+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF")
### transcription factors
library(Seurat)
library(ggplot2)
## read TF data
TF_DB <- read.table("~/Danio_rerio_TF.txt", header = T, sep = "\t")

Time_sig_gene <- subset(Time_diff, qval < 0.00001)

Time_sig_gene_TF <- Time_sig_gene[rownames(Time_sig_gene) %in%  TF_DB$Symbol, ]

pdf("~/Time_sig_gene_names_TF_heatmap.pdf", height = 15)
plot_pseudotime_heatmap(monocle_cds[rownames(Time_sig_gene_TF), ],   show_rownames = TRUE,  num_clusters = 4, return_heatmap = T, cluster_rows = T)
dev.off()

pdf("~/Time_sig_gene_names_TF_reName_heatmap.pdf", height = 8)
plot_pseudotime_heatmap(monocle_cds[rownames(Time_sig_gene_TF), ],   show_rownames = F,  num_clusters = 4, return_heatmap = T, cluster_rows = T)
dev.off()


Time_sig_gene_TF_order <- Time_sig_gene_TF[order(Time_sig_gene_TF$qval), ]

pdf("~/Time_sig_gene_names_TF_heatmap_top100.pdf", height = 8)
plot_pseudotime_heatmap(monocle_cds[rownames(Time_sig_gene_TF_order)[1:100], ],   show_rownames = TRUE,  num_clusters = 4, return_heatmap = T, cluster_rows = T)
dev.off()

#
# getting the heatmap cluster genes 2 clusters
pdf("~/Time_sig_genes_heatmapClusterN2.pdf", height = 8)

p_heatmap <- plot_pseudotime_heatmap(monocle_cds[Time_sig_gene_names, ],   show_rownames = F,  num_clusters = 2, return_heatmap = T)
dev.off()
p_heatmap <- plot_pseudotime_heatmap(monocle_cds[Time_sig_gene_names, ],   show_rownames = F,  num_clusters = 2, return_heatmap = T)
clusters <- cutree(p_heatmap$tree_row, k = 2)
df_clusters <- data.frame(clusters)
df_clusters$gene_names <- rownames(df_clusters)

Time_sig_gene_names_clusters <- merge(df_Time_sig_gene, df_clusters, by.x = "gene_short_name", by.y = "gene_names", sort = F)
Time_sig_gene_names_clusters_order <- Time_sig_gene_names_clusters[order(Time_sig_gene_names_clusters$qval), ]

Time_sig_gene_names_clusters_order_list <- list()
for(i in 1:2){
    Time_sig_gene_names_clusters_order_list[[paste("Cluster", i, sep = "")]]<- Time_sig_gene_names_clusters_order[Time_sig_gene_names_clusters_order$clusters == i, ]

}

writexl::write_xlsx(Time_sig_gene_names_clusters_order_list, "~/LFP_ISN_monocle2_trajectory_Sig_gene_Heatmap_Gene_clusters2.xlsx")
# 
# add Ensembl
features <- read.table("~/features.tsv.gz", header = FALSE, sep = "\t")

getUniqs = function(xx){
  idx2 = duplicated(xx)
  for(i in seq(1, length(idx2),1)) if(idx2[i]) xx[i] = paste0(xx[i],"#",i)
  return(xx)
}

features$gene <- getUniqs(features$V2)
colnames(features) <- c("Ensembl", "V2", "V3", "gene_short_name")
#
Time_sig_gene_names_cluster1 <- readxl::read_xlsx("~/LFP_ISN_monocle2_trajectory_Sig_gene_Heatmap_Gene_clusters2.xlsx")
Time_sig_gene_names_cluster2 <- readxl::read_xlsx("~/LFP_ISN_monocle2_trajectory_Sig_gene_Heatmap_Gene_clusters2.xlsx", sheet = 2)
df_Time_sig_gene_names_cluster1 <- data.frame(Time_sig_gene_names_cluster1)
df_Time_sig_gene_names_cluster2 <- data.frame(Time_sig_gene_names_cluster2)
df_Time_sig_gene_names_cluster1_ensembl <- merge(df_Time_sig_gene_names_cluster1, features[, c("Ensembl", "gene_short_name")], by = "gene_short_name", sort = FALSE)
df_Time_sig_gene_names_cluster2_ensembl <- merge(df_Time_sig_gene_names_cluster2, features[, c("Ensembl", "gene_short_name")], by = "gene_short_name", sort = FALSE)

#
Time_sig_gene_names_clusters_order_list <- list()
Time_sig_gene_names_clusters_order_list[["Cluster1"]] <- df_Time_sig_gene_names_cluster1_ensembl
Time_sig_gene_names_clusters_order_list[["Cluster2"]] <- df_Time_sig_gene_names_cluster2_ensembl

writexl::write_xlsx(Time_sig_gene_names_clusters_order_list, "~/LFP_ISN_monocle2_trajectory_Sig_gene_Heatmap_Gene_clusters2.xlsx")



### sigfincate gene (qval < 0.00001), heatmap 2 cluster to do GO KEGG 
###
## cluster to do kegg
library(stringr)
df_clusters_split <- str_split_fixed(df_clusters$gene_names, "#", 2)
df_clusters$Gene_names <- str_split_fixed(df_clusters$gene_names, "#", 2)[, 1]
#
library("org.Dr.eg.db")
library(clusterProfiler)


gene_ids <- bitr(df_clusters$Gene_names,'SYMBOL','ENTREZID','org.Dr.eg.db')

df_clusters_ENTREZID <- merge(df_clusters, gene_ids, by.x='Gene_names', by.y='SYMBOL', sort = F)
df_clusters_ENTREZID_uniq <- unique(df_clusters_ENTREZID[, c("ENTREZID", "clusters", "Gene_names")])
gcSample <- split(df_clusters_ENTREZID_uniq$ENTREZID, df_clusters_ENTREZID_uniq$clusters)
gcSample # entrez id , compareCluster 
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="dre", pvalueCutoff= 1, qvalueCutoff = 1 )
xx <- setReadable(xx, OrgDb = org.Dr.eg.db, keyType="ENTREZID")
df_xx <- data.frame(xx)

#
#  save markers
Time_sig_genes_heatmap.KEGG.list <- list()
for(i in sort(unique(as.character(df_xx$Cluster)))){
 
    Time_sig_genes_heatmap.KEGG.list[[i]] <- df_xx[df_xx$Cluster == i, ]

}

writexl::write_xlsx(Time_sig_genes_heatmap.KEGG.list,"~/Time_sig_genes_heatmap_AllCluster_KEGG.xlsx")
## filter 
xx_filter <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="dre", pvalueCutoff= 0.05, qvalueCutoff = 0.05 )
xx_filter <- setReadable(xx_filter, OrgDb = org.Dr.eg.db, keyType="ENTREZID")
df_xx_filter <- data.frame(xx_filter)
#
p <- dotplot(xx_filter, showCategory = length(unique(data.frame(df_xx_filter)$Description))) 
p_KEGG <- p + theme(axis.text.x = element_text(angle = 45, 
                                  vjust = 0.5, hjust=0.5))
pdf("~/Time_sig_genes_heatmap_KEGG.pdf", width = 10, height = 10)
print(p_KEGG)

dev.off()
#
p <- dotplot(xx_filter, showCategory = 10) 
p_KEGG <- p + theme(axis.text.x = element_text(angle = 45, 
                                  vjust = 0.5, hjust=0.5))
pdf("~/Time_sig_genes_heatmap_KEGG.pdf", width = 10, height = 10)
print(p_KEGG)

dev.off()

###
xx_Go3_BP <- compareCluster(gcSample, fun="enrichGO", OrgDb = org.Dr.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'ENTREZID' , readable = TRUE )

df_xx_Go3_BP <- data.frame(xx_Go3_BP)
Time_sig_genes_heatmap.Go_BP.list <- list()
for(i in sort(unique(as.character(df_xx_Go3_BP$Cluster)))){
 
    Time_sig_genes_heatmap.Go_BP.list[[i]] <- df_xx_Go3_BP[df_xx_Go3_BP$Cluster == i, ]

}
writexl::write_xlsx(Time_sig_genes_heatmap.Go_BP.list, "~/Time_sig_genes_heatmap_AllCluster_Go_BP.xlsx")


#
p_Go_BP <- dotplot(xx_Go3_BP, showCategory = 20) 

p_Go_BP_1 <- p_Go_BP + theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))

pdf("~/Time_sig_genes_heatmap_GO_BP_top20.pdf", width = 10, height = 13)
print(p_Go_BP_1)

dev.off()
#
############# publication figure
library(monocle)
monocle_cds <- readRDS("~/monocle_cds.rds") 

pdf("~/Publication_LFP_ISN_Trajectory.pdf", width = 10)
plot_cell_trajectory(monocle_cds, color_by = "Cell_labels_New", cell_size = 5)+ scale_color_manual(values = c("LFP" = '#57C3F3', "ISN_pre" = "#91D0BE", "ISN" = '#476D87')) + theme(text = element_text(size = 30), legend.text = element_text(size=30)) + guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

# stage:
pdf("~/Publication_LFP_ISN_Trajectory_stage.pdf", width = 10)
plot_cell_trajectory(monocle_cds, color_by = "orig.ident", cell_size = 5) + scale_color_manual(values = c("#999999", "#333333", "#ff0000", "#00cdcd")) + theme(text = element_text(size = 30), legend.text = element_text(size=30)) + guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

#monocle_cds <- reduceDimension(monocle_cds, max_components = 2, reduction_method = 'DDRTree', sigma = 0.25)
p1 <- plot(plot_cell_trajectory(monocle_cds, show_cell_names = F, color_by = "Pseudotime", cell_size = 5) + scale_color_viridis_c()) + theme(text = element_text(size = 30), legend.text = element_text(size=30),  legend.key.size = unit(1.5, 'cm')) + easy_remove_axes() + easy_remove_x_axis() + easy_remove_y_axis()

pdf("~/Publication_LFP_ISN_Pseudotime.pdf", width = 10)
print(p1)
dev.off()

##
pdf("~/Publication_Genes_In_Pseudotime2.pdf", width = 15)

plot_genes_in_pseudotime(monocle_cds[c("nkx2.9", "lmx1bb", "nkx2.2a", "fev",  "foxa","tph2")], cell_size = 2, color_by = "Cell_labels_New", label_by_short_name =F, ncol = 2) + scale_color_manual(values = c("LFP" = '#57C3F3', "ISN_pre" = "#91D0BE", "ISN" = '#476D87'), name = "Cell_labels") + theme(text = element_text(size = 30), legend.text = element_text(size=30)) + guides(colour = guide_legend(override.aes = list(size=5)))+ geom_line(aes(x = Pseudotime, y = expectation), size = 1.5) + theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1), axis.ticks = element_line(size = 1)) + xlab("Pseudotime")
dev.off()


## plot gene in pseudotime

pData(monocle_cds)$fev <- log2(exprs(monocle_cds)['fev', ] + 1)
pData(monocle_cds)$tph2 <- log2(exprs(monocle_cds)['tph2', ] + 1)
pData(monocle_cds)$nkx2.2a <- log2(exprs(monocle_cds)['nkx2.2a', ] + 1)
pData(monocle_cds)$nkx2.9 <- log2(exprs(monocle_cds)['nkx2.9', ] + 1)
pData(monocle_cds)$foxa <- log2(exprs(monocle_cds)['foxa', ] + 1)
pData(monocle_cds)$lmx1bb <- log2(exprs(monocle_cds)['lmx1bb', ] + 1)

#plot_cell_trajectory(monocle_cds, color_by="nkx2.9")+ scale_color_gradient2(name = NULL, high = "#FF0000", mid = "#FFFFFF", low = "#0000FF")
library(gridExtra)
library(ggeasy)

p1 <- plot_cell_trajectory(monocle_cds, color_by="nkx2.9", cell_size = 5)+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF") + theme(text = element_text(size = 30), legend.text = element_text(size=30),  legend.key.size = unit(1, 'cm'))+ easy_remove_axes() + easy_remove_x_axis() + easy_remove_y_axis()
p2 <- plot_cell_trajectory(monocle_cds, color_by="foxa", cell_size = 5)+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF")  + theme(text = element_text(size = 30), legend.text = element_text(size=30),  legend.key.size = unit(1, 'cm'))+ easy_remove_axes() + easy_remove_x_axis() + easy_remove_y_axis()
p3 <- plot_cell_trajectory(monocle_cds, color_by="nkx2.2a", cell_size = 5)+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF")  + theme(text = element_text(size = 30), legend.text = element_text(size=30),  legend.key.size = unit(1, 'cm'))+ easy_remove_axes() + easy_remove_x_axis() + easy_remove_y_axis()
p4 <- plot_cell_trajectory(monocle_cds, color_by="fev", cell_size = 5)+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF") + theme(text = element_text(size = 30), legend.text = element_text(size=30),  legend.key.size = unit(1, 'cm')) + easy_remove_axes() + easy_remove_x_axis() + easy_remove_y_axis()
p5 <- plot_cell_trajectory(monocle_cds, color_by="lmx1bb", cell_size = 5)+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF") + theme(text = element_text(size = 30), legend.text = element_text(size=30),  legend.key.size = unit(1, 'cm'))+ easy_remove_axes() + easy_remove_x_axis() + easy_remove_y_axis()
p6 <- plot_cell_trajectory(monocle_cds, color_by="tph2", cell_size = 5)+ scale_color_gradient2(high = "#FF0000", mid = "grey", low = "#0000FF") + theme(text = element_text(size = 30), legend.text = element_text(size=30),  legend.key.size = unit(1, 'cm'))+ easy_remove_axes() + easy_remove_x_axis() + easy_remove_y_axis()


pdf("~/Publication_Genes_In_Trajectory.pdf", width = 16, height = 10)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2)

dev.off()
###############################

# Plot LFP, ISN_pre, ISN
pdf("~/LFP_ISN_DimPlot_split_Stage.pdf", width = 15 * 2.5, height = 12)
DimPlot(NP_LFP_ISNpre_ISN_obj_sel, label =F, group.by = "Cell_labels_New", pt.size = 2 , split.by = "orig.ident") + scale_color_manual(values = c("LFP" = '#57C3F3', "ISN_pre" = "#91D0BE", "ISN" = '#476D87')) + theme(text = element_text(size = 30), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=30)) + ggtitle("")

dev.off()


# Heatmap of Significal genes 
library(tidyverse)
library(readxl)

ordering_genes <- readRDS("~/LFP_ISN_monocle2_trajectory_DEGs_qval0.01_ordering_genes.rds" )
Time_diff <- read_xlsx("~/LFP_ISN_monocle2_trajectory_DEGs_qval0.01_todo_Time_diff.xlsx")
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
df_Time_sig_gene <- subset(Time_diff, qval < 0.00001)
Time_sig_gene_names <- subset(Time_diff, qval < 0.00001)$gene_short_name
plot_pseudotime_heatmap(monocle_cds[Time_sig_gene_names, ],   show_rownames = F,  num_clusters = 2) # fix it in the inkscape

## 

pdf("~/Publication_TopClustersGenes_In_Pseudotime.pdf", width = 25)

plot_genes_in_pseudotime(monocle_cds[c("her12", "her2", "notch3", "her6", "syt1a", "scg2b","insm1a", "ddc", "oprl1",  "gch1", "gata3", "gata2a")], cell_size = 2, color_by = "Cell_labels_New", label_by_short_name = F, ncol = 4) + scale_color_manual(values = c("LFP" = '#57C3F3', "ISN_pre" = "#91D0BE", "ISN" = '#476D87'), name = "Cell_labels") + theme(text = element_text(size = 30), legend.text = element_text(size=30)) + guides(colour = guide_legend(override.aes = list(size=5)))+ geom_line(aes(x = Pseudotime, y = expectation), size = 1.5) + theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1), axis.ticks = element_line(size = 1)) + xlab("Pseudotime")
dev.off()

##
saveRDS(gcSample, "~/gcSample.rds")

# KEGG
xx_filter <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="dre", pvalueCutoff= 0.05, qvalueCutoff = 0.05 )
xx_filter <- setReadable(xx_filter, OrgDb = org.Dr.eg.db, keyType="ENTREZID")
df_xx_filter <- data.frame(xx_filter)
# show all KEGG
p <- dotplot(xx_filter, showCategory = length(unique(data.frame(df_xx_filter)$Description))) 
p_KEGG <- p + theme(axis.text.x = element_text(angle = 45, 
                                  vjust = 0.5, hjust=0.5))
pdf("~/Time_sig_genes_heatmap_KEGG.pdf", width = 10, height = 10)
print(p_KEGG)

dev.off()
# show top10 KEGG
p <- dotplot(xx_filter, showCategory = 10) 
p_KEGG <- p + theme(axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30), axis.title=element_text(size=30,face="bold")) + theme(text = element_text(size = 20), legend.text = element_text(size=20)) 
pdf("~/Time_sig_genes_heatmap_Top10_KEGG.pdf", width = 15, height = 12)
print(p_KEGG)

dev.off()
## show top10 GO

xx_Go3_BP <- compareCluster(gcSample, fun="enrichGO", OrgDb = org.Dr.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05,keyType = 'ENTREZID' , readable = TRUE )
#
p_Go_BP <- dotplot(xx_Go3_BP, showCategory = 10) 

p_Go_BP_1 <- p_Go_BP + theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
p_Go_BP_1 <- p_Go_BP + theme(axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30), axis.title=element_text(size=30,face="bold")) + theme(text = element_text(size = 20), legend.text = element_text(size=20)) 

pdf("~/Time_sig_genes_heatmap_GO_BP_top10.pdf", width = 17, height = 12)
print(p_Go_BP_1)

dev.off()


pdf("~/Time_sig_genes_heatmap_GO_BP_top10_scaleY.pdf", width = 15, height = 12)
print(p_Go_BP_1 + scale_y_discrete(labels=function(x) str_wrap(x, width=40)))

dev.off()

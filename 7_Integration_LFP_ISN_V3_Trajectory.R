###### building the LFP -- V3 and -- ISN trajectory by using the Monocle2 and Monocle3, and comparing the similarities and differences.
##

monocle_func <- function(NP_obj_LFP_ISN, qvalue = 0.0001, reverse_log = TRUE){
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
	## bulid trajectory 
	#diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes, ], fullModelFormulaStr = "~Cell_labels")

	diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes, ], fullModelFormulaStr = "~Cell_labels + orig.ident", reducedModelFormulaStr = " ~orig.ident", relative_expr=TRUE)
	#diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes, ], fullModelFormulaStr = "~Cell_labels + data_source", reducedModelFormulaStr = " ~data_source", relative_expr=TRUE)
	ordering_genes <- row.names(subset(diff_test_res, qval < qvalue)) # qval < 0.1, qval < 0.01

	monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)

	#plot_ordering_genes(monocle_cds)

	# reduce data dimensionality

	monocle_cds <- reduceDimension(monocle_cds, max_components = 2, method = 'DDRTree')
	#monocle_cds <- reduceDimension(monocle_cds, max_components = 3, reduction_method = 'DDRTree', residualModelFormulaStr = "~orig.ident")

	monocle_cds <- orderCells(monocle_cds, reverse = reverse_log)

	#plot_cell_trajectory(monocle_cds, color_by = "Cell_labels")
	return(monocle_cds)
}


getIntegratedResults <- function(the.obj.list = NULL, the_samples = NULL, num_dims = 30, vars_regress = c("nCount_RNA"), resoLs = c(0.5, 1.0, 1.5, 2.0), anchor_features = 2000){
  print("Find Anchors ...")
  obj.combined <- FindIntegrationAnchors(object.list = the.obj.list[the_samples], dims = 1:num_dims, verbose = FALSE, anchor.features = anchor_features)
  print("Integrating ...")
  obj.combined <- IntegrateData(anchorset = obj.combined, dims = 1:num_dims, verbose = FALSE)
  DefaultAssay(obj.combined) <- "integrated"
  # Run the standard workflow for visualization and clustering
  all.genes <- rownames(obj.combined)
  print("Scaling ...")
  obj.combined <- ScaleData(obj.combined, features = all.genes, vars.to.regress = vars_regress, verbose = FALSE)
  print("Running PCA ...")
  obj.combined <- RunPCA(obj.combined, npcs = num_dims, verbose = FALSE)
  # t-SNE and Clustering
  obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:num_dims, verbose = FALSE)
  obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:num_dims, verbose = FALSE)
  for(resoL in resoLs) obj.combined <- FindClusters(obj.combined, resolution = resoL, verbose = FALSE)
  DefaultAssay(obj.combined) <- "RNA"
  return(obj.combined)
}
############## Integration of LFP V3 and ISN cells by using Olig2 and Sox1a datasets ##############################################
NP_LFP_ISNpre_ISN_obj_sel <- readRDS(file = "~/NP_LFP_ISN_pre_ISN_filterlabel_LFP_ISN_new.rds")
Olig2_V3_inteagrated_obj <- readRDS(file = "~/Olig2_V3_inteagrated_obj.rds")

Sox1_Olig2_V3_ISN_obj_list <- list()
Sox1_Olig2_V3_ISN_obj_list$Sox1a <- NP_LFP_ISNpre_ISN_obj_sel
Sox1_Olig2_V3_ISN_obj_list$Olig2 <- Olig2_V3_inteagrated_obj

Sox1_Olig2_V3_ISN_obj_Integrate6 <- getIntegratedResults(Sox1_Olig2_V3_ISN_obj_list, c("Sox1a", "Olig2" ), num_dims = 25, anchor_features = 1500, resoLs = c(0.5, 1.0, 1.5, 2.0, 2.5, 3))

DefaultAssay(Sox1_Olig2_V3_ISN_obj_Integrate6) <- "integrated"
Sox1_Olig2_V3_ISN_obj_Integrate6 <- FindClusters(Sox1_Olig2_V3_ISN_obj_Integrate6, resolution = 2.5, verbose = FALSE)
Sox1_Olig2_V3_ISN_obj_Integrate6 <- FindClusters(Sox1_Olig2_V3_ISN_obj_Integrate6, resolution = 3, verbose = FALSE)

DefaultAssay(Sox1_Olig2_V3_ISN_obj_Integrate6) <- "RNA"


DimPlot(Sox1_Olig2_V3_ISN_obj_Integrate6, label = T) + DimPlot(Sox1_Olig2_V3_ISN_obj_Integrate6, group.by = "orig.ident") +  FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, c("nkx2.2a", "nkx2.9", "dlb", "fev", "tph2", "sim1a"))

## define cell classes
Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New <- "LFP_1"
Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New[Sox1_Olig2_V3_ISN_obj_Integrate6$integrated_snn_res.3 %in% c(5, 8)] <- "LFP_2"
Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New[Sox1_Olig2_V3_ISN_obj_Integrate6$integrated_snn_res.3 %in% c(16)] <- "MN"
Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New[Sox1_Olig2_V3_ISN_obj_Integrate6$integrated_snn_res.3 %in% c(9)] <- "NPCs" # neuroanl precursor cells
Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New[Sox1_Olig2_V3_ISN_obj_Integrate6$integrated_snn_res.3 %in% c(7)] <- "ISN_pre"
Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New[Sox1_Olig2_V3_ISN_obj_Integrate6$integrated_snn_res.3 %in% c(15, 0, 17)] <- "ISN"
Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New[Sox1_Olig2_V3_ISN_obj_Integrate6$integrated_snn_res.3 %in% c(3, 4)] <- "V3_pre"
Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New[Sox1_Olig2_V3_ISN_obj_Integrate6$integrated_snn_res.3 %in% c(13)] <- "V3"


saveRDS(Sox1_Olig2_V3_ISN_obj_Integrate6, file = "~/Sox1_Olig2_V3_ISN_obj_25_1500_New.rds")
Sox1_Olig2_V3_ISN_obj_Integrate6_sel <- subset(Sox1_Olig2_V3_ISN_obj_Integrate6, Cell_labels_New %in% c("LFP_1", "NPCs", "ISN_pre", "ISN", "V3_pre", "V3"))

#
library(ggeasy)
library(gridExtra)

nkx2.9_plot <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "nkx2.9", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

nkx2.2a_plot <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "nkx2.2a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

nkx2.2b_plot  <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "nkx2.2b", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()


pdf("~/Sox1a_Olig2_LFP_V3_ISN_integration_FeaturePlot_nkx2.9_nkx2.2a_nkx2.2b.pdf", width = 20* 3, height = 15)
grid.arrange(nkx2.9_plot, nkx2.2a_plot, nkx2.2b_plot, ncol=3, nrow=1)

dev.off()

##
pcna_plot <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "pcna", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

mki67_plot <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "mki67", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
lig1_plot  <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "lig1", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
sulf1_plot  <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "sulf1", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

gfap_plot  <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "gfap", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
fabp7a_plot  <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "fabp7a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
shha_plot  <- FeaturePlot(Sox1_Olig2_V3_ISN_obj_Integrate6, "shha", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()


pdf("~/Sox1a_Olig2_ISN_V3_integration_FeaturePlot_pcna_mki67_lig1_gfap_fabp7a_shha.pdf", width = 20* 3, height = 15*2)
grid.arrange(pcna_plot, mki67_plot, lig1_plot, gfap_plot, fabp7a_plot, shha_plot, ncol=3, nrow=2)

dev.off()


pdf("~/Sox1a_Olig2_ISN_V3_integration_FeaturePlot_pcna_mki67_sulf1_gfap_nkx2.2a_shha.pdf", width = 20* 3, height = 15*2)
grid.arrange(pcna_plot, mki67_plot, sulf1_plot, gfap_plot, nkx2.2a_plot, shha_plot, ncol=3, nrow=2)

dev.off()

pdf("~/Sox1a_Olig2_ISN_V3_integration_FeaturePlot_nkx2.2a_nkx2.9_pcna_mki67_gfap_shha.pdf", width = 20* 3, height = 15*2)
grid.arrange(nkx2.2a_plot, nkx2.9_plot, pcna_plot, mki67_plot, gfap_plot, shha_plot, ncol=3, nrow=2)

dev.off()

# monocle analyze the ISN  and V3
Sox1_Olig2_V3_ISN_obj_Integrate6_sel$Cell_labels <- Sox1_Olig2_V3_ISN_obj_Integrate6_sel$Cell_labels_New
Sox1_Olig2_V3_ISN_obj_Integrate6_sel$data_source <- "Olig2"
Sox1_Olig2_V3_ISN_obj_Integrate6_sel$data_source[Sox1_Olig2_V3_ISN_obj_Integrate6_sel$orig.ident %in% c("1dpf","2dpf", "3dpf", "5dpf")] <- "Sox1a"

NP_obj_LFP_ISN <- Sox1_Olig2_V3_ISN_obj_Integrate6_sel

#monocle_cds_0.01 <- monocle_func(NP_obj_LFP_ISN, qvalue = 0.01, reverse_log = FALSE) # good
#monocle_cds_0.05 <- monocle_func(NP_obj_LFP_ISN, qvalue = 0.05, reverse_log = FALSE)
#monocle_cds_0.009 <- monocle_func(NP_obj_LFP_ISN, qvalue = 0.009, reverse_log = TRUE)

#plot_cell_trajectory(monocle_cds_0.009, color_by = "Cell_labels") + plot_cell_trajectory(monocle_cds_0.009, color_by = "orig.ident")
#plot_genes_branched_pseudotime(monocle_cds_0.009[c("nkx2.9", "nkx2.2a", "fev", "foxa", "sim1a", "tph2")], color_by = "Cell_labels", branch_point = 1, branch_labels = c("ISN", "V3"))
##
##
LFP_ISN_data <- as(as.matrix(NP_obj_LFP_ISN@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = NP_obj_LFP_ISN@meta.data)
fData <- data.frame(gene_short_name = row.names(LFP_ISN_data), row.names = row.names(LFP_ISN_data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(LFP_ISN_data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

#monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

print(head(fData(monocle_cds))) # 23191

expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= 10)) # 12696 
print(head(pData(monocle_cds)))

# using Seurat cluster markers to order cells
##
Idents(Sox1_Olig2_V3_ISN_obj_Integrate6_sel) <- Sox1_Olig2_V3_ISN_obj_Integrate6_sel$Cell_labels_New
V3_ISN_deg.cluster <- FindAllMarkers(Sox1_Olig2_V3_ISN_obj_Integrate6_sel)
diff.genes <- subset(V3_ISN_deg.cluster, p_val_adj < 0.001)
library(dplyr)

##  order diff genes then use top 500 genes ########################
diff.genes_order <- diff.genes[order(diff.genes$p_val_adj), ]
diff.genes_order_top500 <- diff.genes_order$gene[1:500]
monocle_cds <- setOrderingFilter(monocle_cds, diff.genes_order_top500)
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, reduction_method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds, reverse = T)

saveRDS(monocle_cds, "~/LFP_ISN_V3_trajectory_monocle_cds.rds")
saveRDS(diff.genes, "~/LFP_ISN_V3_trajectory_diff.genes.rds")

plot_cell_trajectory(monocle_cds, color_by = "Cell_labels_New") + plot_cell_trajectory(monocle_cds, color_by = "orig.ident")

plot_genes_branched_pseudotime(monocle_cds[c("nkx2.9", "nkx2.2a", "fev", "foxa", "sim1a", "tph2")], color_by = "Cell_labels_New", branch_point = 1, branch_labels = c( "V3", "ISN"))


## BEAM
BEAM_res <- BEAM(monocle_cds, branch_point = 1, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]

plot_genes_branched_heatmap(monocle_cds[rownames(subset(BEAM_res, qval < 0.00001)), ], branch_point = 1, num_clusters = 3, cores = 1, use_gene_short_name = T, show_rownames = T)

saveRDS(BEAM_res, "~/LFP_ISN_V3_trajectory_BEAM_res.rds")
#
## filter some bad genes(special in data source. dpf24h from olig2 vs 1dpf from sox1a)
Sox1_Olig2_V3_ISN_obj_Integrate6_sel$Cell_label_Stage_New <- paste(Sox1_Olig2_V3_ISN_obj_Integrate6_sel$Cell_labels_New, Sox1_Olig2_V3_ISN_obj_Integrate6_sel$orig.ident, sep ="_")
Idents(Sox1_Olig2_V3_ISN_obj_Integrate6_sel) <- Sox1_Olig2_V3_ISN_obj_Integrate6_sel$Cell_label_Stage_New
#  LFP_1_1dpf vs LFP_1_dpf24h diff
LFP_1dpf_24hpf_markers <- FindMarkers(Sox1_Olig2_V3_ISN_obj_Integrate6_sel, ident.1 = "LFP_1_1dpf", ident.2 = "LFP_1_dpf24h", verbose = FALSE)
LFP_1dpf_24hpf_markers_sel <- LFP_1dpf_24hpf_markers[LFP_1dpf_24hpf_markers$p_val_adj < 0.0001,]

#
BEAM_res_sel <- BEAM_res[!c(BEAM_res$gene_short_name %in% c(rownames(LFP_1dpf_24hpf_markers_sel))), ]
BEAM_res_sel2 <- BEAM_res_sel[BEAM_res_sel$qval < 0.01, ]
BEAM_res_sel3 <- BEAM_res_sel2[!c(substr(BEAM_res_sel2$gene_short_name, 1, 2) %in% c("rp", "NC")), ]
plot_genes_branched_heatmap(monocle_cds[rownames(BEAM_res_sel3), ], branch_point = 1, branch_labels = c("V3 Cell fate", "ISN Cell fate"),num_clusters = 4, cores = 1, use_gene_short_name = T, show_rownames = F)
BEAM_res_sel3_top100 <- BEAM_res_sel3[1:100, ]
pdf("~/V3_ISN_BEAM_res_sel3_top100.pdf", height = 10)
plot_genes_branched_heatmap(monocle_cds[rownames(BEAM_res_sel3_top100), ], branch_point = 1, branch_labels = c("V3 Cell fate", "ISN Cell fate"),num_clusters = 4, cores = 1, use_gene_short_name = T, show_rownames = T)

#print(p1)
dev.off()

### transcription factors
library(Seurat)
library(ggplot2)
## read TF data
TF_DB <- read.table("~/Danio_rerio_TF.txt", header = T, sep = "\t")

BEAM_res_sel3_TF <- BEAM_res_sel3[rownames(BEAM_res_sel3) %in%  TF_DB$Symbol, ]

pdf("~/V3_ISN_BEAM_res_sel3_TF_heatmap.pdf", height = 10)
plot_genes_branched_heatmap(monocle_cds[rownames(BEAM_res_sel3_TF), ], branch_point = 1, branch_labels = c("V3 Cell fate", "ISN Cell fate"),num_clusters = 4, cores = 1, use_gene_short_name = T, show_rownames = T)
dev.off()

################# Publication plot 

Sox1_Olig2_V3_ISN_obj_Integrate6 <- readRDS(file = "~/Sox1_Olig2_V3_ISN_obj_25_1500_New.rds")

### Cell Label to find markers
Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New <- factor(Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New, levels = c("LFP_1", "LFP_2","NPCs", "ISN_pre", "ISN", "V3_pre", "V3", "MN"))

Idents(Sox1_Olig2_V3_ISN_obj_Integrate6) <- Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New
##
## cell number in these 
# the bar plot for cell numbers at every stage:
library(ggplot2)

## cell number and cell types
df_cell_number_cell_types <- data.frame(table( Sox1_Olig2_V3_ISN_obj_Integrate6$orig.ident, Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New))
colnames(df_cell_number_cell_types) <- c("Stages", "Cell_labels", "Num")
df_cell_number_cell_types$Cell_labels <- factor(df_cell_number_cell_types$Cell_labels,  levels = c("LFP_1", "LFP_2","NPCs", "ISN_pre", "ISN", "V3_pre", "V3", "MN"))

# contribution of each stage to the cell types identified
p_bar_cellPercent_cellTypes2 <- ggplot(data=df_cell_number_cell_types, aes(x=Cell_labels, y=Num, fill=Stages)) +
  geom_bar(stat="identity", position="fill")+ scale_fill_manual(values = c( "1dpf" = "#999999", "2dpf" = "#333333", "3dpf" = "#ff0000", "5dpf"  = "#00cdcd", "dpf24h" = "#BD956A", "dpf36h" = '#8C549C' , "dpf48h"= '#E95C59')) + theme_minimal() + xlab("") + ylab("Cell Percentage")
pdf("~/V3_ISN_CellsStagesPercent_CellTypes_bar.pdf")
print(p_bar_cellPercent_cellTypes2)
dev.off()


df_cell_num <- data.frame(table(Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New))
colnames(df_cell_num) <- c("Cell_labels", "Num")
df_cell_num$Cell_labels <- factor(df_cell_num$Cell_labels, levels = c("LFP_1", "LFP_2","NPCs", "ISN_pre", "ISN", "V3_pre", "V3", "MN"))

p_bar_cell_num <- ggplot(df_cell_num, aes(x=Cell_labels, y=Num, fill=Cell_labels))+ geom_bar(stat="identity")+ scale_fill_manual(values = c("LFP_1" = '#57C3F3', "NPCs" = "#53A85F","ISN_pre" = "#91D0BE", "ISN" = '#476D87', "V3_pre" = "#E59CC4", "V3" = "#E95C59")) + geom_text(aes(label=Num), vjust=-0.25, position = position_dodge(0.9), size=5) + theme_minimal() + xlab("") + ylab("Cell Numbers") + theme(legend.position = "none")
pdf("~/V3_ISN_Cells_Number_bar.pdf")
print(p_bar_cell_num)
dev.off()
##
Sox1_Olig2_V3_ISN_obj_Integrate6.all.mrkrs = FindAllMarkers(Sox1_Olig2_V3_ISN_obj_Integrate6, only.pos = T)

Sox1_Olig2_V3_ISN_obj_Integrate6_CellType.mrkrs.list <- list()
#for(i in sort(as.character(unique(Cell_NP_obj.all.mrkrs$cluster)))){
for(i in as.character(unique(sort(Sox1_Olig2_V3_ISN_obj_Integrate6$Cell_labels_New)))){
 
    Sox1_Olig2_V3_ISN_obj_Integrate6_CellType.mrkrs.list[[i]] <- Sox1_Olig2_V3_ISN_obj_Integrate6.all.mrkrs[Sox1_Olig2_V3_ISN_obj_Integrate6.all.mrkrs$cluster == i, ]

}



writexl::write_xlsx(Sox1_Olig2_V3_ISN_obj_Integrate6_CellType.mrkrs.list, "~/Sox1a_Olig2_LFP_V3_ISN_Integrated_CellType_markers.xlsx")

####
pdf("~/Sox1_Olig2_V3_ISN_integration_DimPlot.pdf", height = 10, width = 14)

DimPlot(Sox1_Olig2_V3_ISN_obj_Integrate6, group.by = "Cell_labels_New", label = T, pt.size = 3, label.size = 10) + scale_color_manual(values = c("LFP_1" = '#57C3F3', "LFP_2" = "#58A4C3", "NPCs" = "#53A85F","ISN_pre" = "#91D0BE", "ISN" = '#476D87', "MN" = '#E0D4CA', "V3_pre" = "#E59CC4", "V3" = "#E95C59")) + theme(text = element_text(size = 30), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=30)) + ggtitle("")  + NoLegend()


dev.off()
## stage

pdf("~/Sox1_Olig2_V3_ISN_integration_DimPlot_Stage.pdf", height = 10, width = 14)

DimPlot(Sox1_Olig2_V3_ISN_obj_Integrate6, group.by = "orig.ident", label = F, pt.size = 3, label.size = 10) + scale_color_manual(values = c( "1dpf" = "#999999", "2dpf" = "#333333", "3dpf" = "#ff0000", "5dpf"  = "#00cdcd", "dpf24h" = "#BD956A", "dpf36h" = '#8C549C' , "dpf48h"= '#E95C59')) + theme(text = element_text(size = 30), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=30)) + ggtitle("")  

dev.off()
## trajectory:
monocle_cds <- readRDS( "~/LFP_ISN_V3_trajectory_monocle_cds.rds")
monocle_cds$Cell_labels <- factor(monocle_cds$Cell_labels, levels = c("LFP_1", "NPCs", "ISN_pre", "ISN", "V3_pre", "V3"))

plot_cell_trajectory(monocle_cds, color_by = "Cell_labels_New") + plot_cell_trajectory(monocle_cds, color_by = "orig.ident")

plot_genes_branched_pseudotime(monocle_cds[c("nkx2.9", "nkx2.2a", "fev", "foxa", "sim1a", "tph2")], color_by = "Cell_labels_New", branch_point = 1, branch_labels = c( "V3", "ISN"))


#############

pdf("~/Olig2_Sox1a_LFP_ISN_V3_Trajectory.pdf", width = 10, height = 10)
plot_cell_trajectory(monocle_cds, color_by = "Cell_labels", cell_size = 3)+ scale_color_manual(values = c("LFP_1" = '#57C3F3', "NPCs" = "#53A85F","ISN_pre" = "#91D0BE", "ISN" = '#476D87', "V3_pre" = "#E59CC4", "V3" = "#E95C59")) + theme(text = element_text(size = 30), legend.text = element_text(size=30)) + guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

# stage:
pdf("~/Olig2_Sox1a_LFP_ISN_V3_Trajectory_stage.pdf", width = 10, height = 10)
plot_cell_trajectory(monocle_cds, color_by = "orig.ident", cell_size = 5)+ scale_color_manual(values = c( "1dpf" = "#999999", "2dpf" = "#333333", "3dpf" = "#ff0000", "5dpf"  = "#00cdcd", "dpf24h" = "#BD956A", "dpf36h" = '#8C549C' , "dpf48h"= '#E95C59')) + theme(text = element_text(size = 30), legend.text = element_text(size=30)) + guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

#
p1 <- plot(plot_cell_trajectory(monocle_cds, show_cell_names = F, color_by = "Pseudotime", cell_size = 5) + scale_color_viridis_c()) + theme(text = element_text(size = 30), legend.text = element_text(size=30),  legend.key.size = unit(1.5, 'cm')) + easy_remove_axes() + easy_remove_x_axis() + easy_remove_y_axis()

pdf("~/Olig2_Sox1a_LFP_ISN_V3_Pseudotime.pdf", width = 10, height = 10)
print(p1)
dev.off()
#

pdf("~/Genes_In_Olig2_Sox1a_LFP_ISN_V3_Pseudotime.pdf", width = 15)

plot_genes_branched_pseudotime(monocle_cds[c("nkx2.9", "foxa", "nkx2.2a", "fev", "tph2", "sim1a")], cell_size = 2, color_by = "Cell_labels", branch_point = 1, branch_labels = c( "V3 Cell", "ISN Cell"), label_by_short_name =F, ncol = 2) + scale_color_manual(values = c("LFP_1" = '#57C3F3', "NPCs" = "#53A85F","ISN_pre" = "#91D0BE", "ISN" = '#476D87', "V3_pre" = "#E59CC4", "V3" = "#E95C59"), name = "Cell_labels") + theme(text = element_text(size = 30), legend.text = element_text(size=30)) + guides(colour = guide_legend(override.aes = list(size=5))) + theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1), axis.ticks = element_line(size = 1)) + xlab("Pseudotime")
dev.off()

##
BEAM_res <- readRDS("~/LFP_ISN_V3_trajectory_BEAM_res.rds")

## filter some bad genes(special in data source. dpf24h from olig2 vs 1dpf from sox1a)
Sox1_Olig2_V3_ISN_obj_Integrate6_sel$Cell_label_Stage_New <- paste(Sox1_Olig2_V3_ISN_obj_Integrate6_sel$Cell_labels_New, Sox1_Olig2_V3_ISN_obj_Integrate6_sel$orig.ident, sep ="_")
Idents(Sox1_Olig2_V3_ISN_obj_Integrate6_sel) <- Sox1_Olig2_V3_ISN_obj_Integrate6_sel$Cell_label_Stage_New
#  LFP_1_1dpf vs LFP_1_dpf24h diff
LFP_1dpf_24hpf_markers <- FindMarkers(Sox1_Olig2_V3_ISN_obj_Integrate6_sel, ident.1 = "LFP_1_1dpf", ident.2 = "LFP_1_dpf24h", verbose = FALSE)
LFP_1dpf_24hpf_markers_sel <- LFP_1dpf_24hpf_markers[LFP_1dpf_24hpf_markers$p_val_adj < 0.0001,]

#
BEAM_res_sel <- BEAM_res[!c(BEAM_res$gene_short_name %in% c(rownames(LFP_1dpf_24hpf_markers_sel))), ]
BEAM_res_sel2 <- BEAM_res_sel[BEAM_res_sel$qval < 0.01, ]
BEAM_res_sel3 <- BEAM_res_sel2[!c(substr(BEAM_res_sel2$gene_short_name, 1, 2) %in% c("rp", "NC")), ]

plot_genes_branched_heatmap(monocle_cds[rownames(BEAM_res_sel3), ], branch_point = 1, branch_labels = c("V3 Cell fate", "ISN Cell fate"),num_clusters = 3, cores = 1, use_gene_short_name = T, show_rownames = F)

plot_genes_branched_heatmap(monocle_cds[rownames(BEAM_res_sel3[BEAM_res_sel3$qval < 1e-3, ]), ], branch_point = 1, branch_labels = c("V3 Cell fate", "ISN Cell fate"),num_clusters = 4, cores = 1, use_gene_short_name = T, show_rownames = F)

###########
### qval < 0.0001

BEAM_res_sel5 <- BEAM_res_sel3[BEAM_res_sel3$qval < 1e-4, ]
pdf("~/LFP_ISN_V3_trajectory_Time_sig_genes_qval0.0001_heatmapClusterN4.pdf", height = 8)

plot_genes_branched_heatmap(monocle_cds[rownames(BEAM_res_sel5), ], branch_point = 1, branch_labels = c("V3 Cell fate", "ISN Cell fate"),num_clusters = 4, cores = 1, use_gene_short_name = T, show_rownames = F)
dev.off()
p_heatmap <- plot_genes_branched_heatmap(monocle_cds[rownames(BEAM_res_sel5), ], branch_point = 1, branch_labels = c("V3 Cell fate", "ISN Cell fate"),num_clusters = 4, cores = 1, use_gene_short_name = T, show_rownames = F, return_heatmap = TRUE)

clusters <- p_heatmap$annotation_row
df_clusters <- data.frame(clusters)
df_clusters$gene_names <- rownames(df_clusters)

BEAM_res_sel5_clusters <- merge(BEAM_res_sel5, df_clusters, by.x = "gene_short_name", by.y = "gene_names", sort = F)
BEAM_res_sel5_clusters_order <- BEAM_res_sel5_clusters[order(BEAM_res_sel5_clusters$qval), ]

BEAM_res_sel5_clusters_order_list <- list()
for(i in 1:4){
    BEAM_res_sel5_clusters_order_list[[paste("Cluster", i, sep = "")]]<- BEAM_res_sel5_clusters_order[BEAM_res_sel5_clusters_order$Cluster == i, ]

}

writexl::write_xlsx(BEAM_res_sel5_clusters_order_list, "~/LFP_ISN_V3_trajectory_BEAM_res_sel5_Heatmap_Gene_ClusterN4.xlsx")


pdf("~/Genes_BEAM_In_Olig2_Sox1a_LFP_ISN_V3_Pseudotime.pdf", width = 15, height = 10)

plot_genes_branched_pseudotime(monocle_cds[c("fev", "gata3", "sim1a", "nos1", "foxp4", "fosab", "nkx2.2a", "pcdh19")], cell_size = 2, color_by = "Cell_labels", branch_point = 1, branch_labels = c( "V3 Cell", "ISN Cell"), label_by_short_name =F, ncol = 2) + scale_color_manual(values = c("LFP_1" = '#57C3F3', "NPCs" = "#53A85F","ISN_pre" = "#91D0BE", "ISN" = '#476D87', "V3_pre" = "#E59CC4", "V3" = "#E95C59"), name = "Cell_labels") + theme(text = element_text(size = 30), legend.text = element_text(size=30)) + guides(colour = guide_legend(override.aes = list(size=5))) + theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1), axis.ticks = element_line(size = 1)) + xlab("Pseudotime")
dev.off()



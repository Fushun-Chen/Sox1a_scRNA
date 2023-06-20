library(Seurat)

GSE173350_data_obj = readRDS("~/GSE173350/GSE173350_appel.int.RDS")

DefaultAssay(GSE173350_data_obj) <- "RNA"

######
pdf("~/Olig2_DimPlot.pdf", height = 12, width = 12)
DimPlot(GSE173350_data_obj, label = T, label.size = 5) + NoLegend()
dev.off()


library(gridExtra)
library(ggeasy)

nkx2.2a_plot <- FeaturePlot(GSE173350_data_obj, "nkx2.2a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

fev_plot <- FeaturePlot(GSE173350_data_obj, "fev", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
lmx1bb_plot  <- FeaturePlot(GSE173350_data_obj, "lmx1bb", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
tph2_plot  <- FeaturePlot(GSE173350_data_obj, "tph2", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

gata2a_plot  <- FeaturePlot(GSE173350_data_obj, "gata2a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()


pdf("~/Olig2_FeaturePlot_nkx2.2a_fev_lmx1bb_tph2.pdf", width = 20* 2, height = 15*2)
grid.arrange(nkx2.2a_plot, fev_plot, lmx1bb_plot, tph2_plot, ncol=2, nrow=2)

dev.off()


pdf("~/Olig2_FeaturePlot_gata2a_fev_lmx1bb_tph2.pdf", width = 20* 2, height = 15*2)
grid.arrange(gata2a_plot, fev_plot, lmx1bb_plot, tph2_plot, ncol=2, nrow=2)

dev.off()
#####

DefaultAssay(GSE173350_data_obj) <- "RNA"
Idents(GSE173350_data_obj) <- GSE173350_data_obj$integrated_snn_res.1.2
##
GSE173350_data_obj <- FindSubCluster(GSE173350_data_obj, cluster = "6", graph.name = "integrated_snn", subcluster.name = "sub.cluster_v3_1", resolution = 0.6, algorithm = 1)

GSE173350_data_obj_apart <- subset(GSE173350_data_obj, integrated_snn_res.1.2 %in% "6")

DimPlot(GSE173350_data_obj_apart, label = T, group.by = "sub.cluster_v3_1") + DimPlot(GSE173350_data_obj_apart, label = T, group.by = "sub.cluster_v3_1", split.by = "orig.ident")

##
Idents(GSE173350_data_obj) <- GSE173350_data_obj$intCluster
GSE173350_data_obj <- FindSubCluster(GSE173350_data_obj, cluster = "V3 IN", graph.name = "integrated_snn", subcluster.name = "V3_SubClusters", resolution = 0.6, algorithm = 1)

GSE173350_data_obj_apart <- subset(GSE173350_data_obj, intCluster %in% "V3 IN")

DimPlot(GSE173350_data_obj_apart, label = T, group.by = "V3_SubClusters", label.size = 5) + NoLegend()

V3_sim1a <- FeaturePlot(GSE173350_data_obj_apart, c("sim1a"), pt.size = 3) + theme(text = element_text(size = 30), legend.key.height= unit(1, 'cm')) 
V3_nkx2.2a <- FeaturePlot(GSE173350_data_obj_apart, c("nkx2.2a"), pt.size = 3) + theme(text = element_text(size = 30), legend.key.height= unit(1, 'cm'))
V3_lhx5 <- FeaturePlot(GSE173350_data_obj_apart, c("lhx5"), pt.size = 3) + theme(text = element_text(size = 30), legend.key.height= unit(1, 'cm'))
V3_DimPlot <- DimPlot(GSE173350_data_obj_apart, label = F, group.by = "V3_SubClusters",  pt.size = 3) + theme(text = element_text(size = 30), legend.key.height= unit(1, 'cm')) 


pdf("~/V3_SubClusters_FeaturePlot_Sim1a.pdf", width = 20, height = 10)
grid.arrange(V3_sim1a, V3_DimPlot, ncol=2, nrow=1)

dev.off()
##


pdf("~/V3_SubClusters_FeaturePlot_Sim1a_nkx2.2a_lhx5_2.pdf", width = 40, height = 10)
grid.arrange(V3_DimPlot, V3_nkx2.2a, V3_sim1a, V3_lhx5,  ncol = 4, nrow=1)

dev.off()
##


##

Gene_names_change_fun <- function(tmp.data, rwns){
    theGenes = readRDS("~/theGenes_all_ensVersions.rds")
    theGenes = lapply(theGenes, function(x){ x = x$ens2sym$dre; colnames(x) = c("dreENS","dreSYM"); return(x)})
    #
    # now arrange rownames
    #table(rwns %in% theGenes$Version_92$dreSYM)
    new.rwns = data.frame("dreENS" = "", "dreSYM" = rwns, "symb" = rwns)
    idx = match(new.rwns$dreSYM,theGenes$Version_92$dreSYM)
    new.rwns$dreENS = theGenes$Version_92$dreENS[idx]
    rem.genes = new.rwns[is.na(new.rwns$dreENS),]
    
    new.rwns = new.rwns[!is.na(new.rwns$dreENS),]
    
    rem.theGenes = lapply(theGenes, function(x) x[!(x$dreENS %in% new.rwns$dreENS),])
    
    # check if any of the genes is in different ensemble versions
    aa = lapply(rem.theGenes, function(x) x[x$dreSYM %in% rem.genes$symb,])
    tt = aa[[1]]
    for(x in names(aa)[2:length(aa)]) tt = rbind(tt, aa[[x]])
    tt = tt[!duplicated(tt$dreENS),]
    tt$symb = tt$dreSYM
    
    new.rwns = rbind(new.rwns, tt)
    rem.genes = rem.genes[!(rem.genes$dreSYM %in% new.rwns$dreSYM),]
    
    # now split genes by space and find any genes that might be in different ensembl versions
    symbs = lapply(rem.genes$dreSYM, function(x) unlist(strsplit(x, split = " "))[1])
    symbs = unlist(symbs)
    rem.genes$symb = symbs
    
    aa = lapply(rem.theGenes, function(x) x[x$dreSYM %in% rem.genes$symb,])
    tt = aa[[1]]
    for(x in names(aa)[2:length(aa)]) tt = rbind(tt, aa[[x]])
    tt = tt[!duplicated(tt$dreENS),]
    tt$symb = tt$dreSYM
    colnames(tt) = colnames(rem.genes)
    
    idx = match(tt$symb, rem.genes$symb)
    tt$dreSYM = rem.genes$dreSYM[idx]
    new.rwns = rbind(new.rwns, tt)
    rem.genes = rem.genes[!(rem.genes$dreSYM %in% new.rwns$dreSYM),]
    
    # now split genes by space and find any genes that might be in different ensembl versions
    symbs = lapply(rem.genes$dreSYM, function(x) unlist(strsplit(x, split = "\\."))[1])
    symbs = unlist(symbs)
    rem.genes$symb = symbs
    
    aa = lapply(rem.theGenes, function(x) x[x$dreSYM %in% rem.genes$symb,])
    tt = aa[[1]]
    for(x in names(aa)[2:length(aa)]) tt = rbind(tt, aa[[x]])
    tt = tt[!duplicated(tt$dreENS),]
    tt$symb = tt$dreSYM
    colnames(tt) = colnames(rem.genes)
    
    idx = match(tt$symb, rem.genes$symb)
    tt$dreSYM = rem.genes$dreSYM[idx]
    new.rwns = rbind(new.rwns, tt)
    rem.genes = rem.genes[!(rem.genes$dreSYM %in% new.rwns$dreSYM),]
    #rem.genes = rem.genes[!(rem.genes$symb %in% new.rwns$symb),]
    
    new.rwns = new.rwns[!duplicated(new.rwns$dreENS),]
    #
    new.ensids = read.delim("~/features.tsv", header = F)
    colnames(new.ensids) = c("dreENS","dreSYM","noo")
    new.ensids = new.ensids[,c(1:2)]
    new.ensids$not.dups = getUniqs(new.ensids$dreSYM)

    #####
    new.rwns = new.rwns[new.rwns$dreENS %in% new.ensids$dreENS,]
    idx = match(new.rwns$dreENS, new.ensids$dreENS)
    new.rwns$newSYM = new.ensids$dreSYM[idx]
    new.rwns$symb = NULL
    ##
    new.rwns_rmduplicate <- new.rwns[!c(duplicated(new.rwns$dreSYM)), ]
    All_rem <- data.frame("dreENS" = "", "dreSYM" = rwns[!c(rwns %in% new.rwns_rmduplicate$dreSYM)], "newSYM" = rwns[!c(rwns %in% new.rwns_rmduplicate$dreSYM)])
    ALL_names <- rbind(new.rwns_rmduplicate, All_rem)
    #idx = match(ALL_names$dreSYM, rownames(tmp_counts))
    idx = match( rownames(tmp.data[[1]]), ALL_names$dreSYM)
    for(trt in names(tmp.data)) rownames(tmp.data[[trt]]) = ALL_names$newSYM[idx]

    return(tmp.data)
    ##
    #idx = match(new.rwns$dreSYM, rownames(tmp_counts))
    #new.counts = lapply(tmp.data, function(x) x[idx,])
    
    #new.ensids = new.ensids[,c(1:2)]
    #
    #new.ensids$not.dups = getUniqs(new.ensids$dreSYM)
    #idx = match(new.rwns$dreENS,new.ensids$dreENS)
    #for(trt in names(new.counts)) rownames(new.counts[[trt]]) = new.ensids$not.dups[idx]
}

######################################################

## get LFP, V3_pre, V3 from olig2 
GSE173350_data_obj_LFP_V3 <- subset(GSE173350_data_obj, sub.cluster_v3_1 %in% c("6_0", "6_1",  "21", "1"))

DimPlot(GSE173350_data_obj_LFP_V3, label = T, group.by = "sub.cluster_v3_1")


GSE173350_data_LFP_V3_count <- GSE173350_data_obj_LFP_V3@assays$RNA@counts

GSE173350_data_LFP_V3_tmp.data = list("dpf24h" = GSE173350_data_LFP_V3_count[,GSE173350_data_obj_LFP_V3$orig.ident %in% "24h"],
"dpf36h" = GSE173350_data_LFP_V3_count[,GSE173350_data_obj_LFP_V3$orig.ident %in% "36h"],
"dpf48h" = GSE173350_data_LFP_V3_count[,GSE173350_data_obj_LFP_V3$orig.ident %in% "48hpf"])

GSE173350_data_LFP_V3_tmp.data = lapply(GSE173350_data_LFP_V3_tmp.data, function(x) as.matrix(x))

##
GSE173350_data_LFP_V3_tmp.data_count_list <- Gene_names_change_fun(GSE173350_data_LFP_V3_tmp.data, rownames(GSE173350_data_LFP_V3_count))

##

V3_rawData_all <- c( GSE173350_data_LFP_V3_tmp.data_count_list)

for(trt in names(V3_rawData_all)){
    tmp = V3_rawData_all[[trt]]
    ngenes = rowSums(tmp > 0)
    tmp = tmp[ngenes >= 3,]
    V3_rawData_all[[trt]] = tmp
}

v3_obj.lists <- list()
    for(trt in names(V3_rawData_all)){
        v3_obj.lists[[trt]] = getObjs(V3_rawData_all[[trt]], regs.out = c("nCount_RNA"), resL = 1, npcas = 30, n.cells = 0, trt = trt)
    }

ldm = function(x) lapply(x, function(x) dim(x))

ldm(v3_obj.lists) # 2dpf 34cells

## Olig2 V3 integrated
Olig2_V3_inteagrated_obj = getIntegratedResults(the.obj.list = v3_obj.lists, the_samples = c("dpf24h", "dpf36h","dpf48h"), vars_regress = c("nCount_RNA"))
DimPlot(Olig2_V3_inteagrated_obj) + FeaturePlot(Olig2_V3_inteagrated_obj, c("nkx2.9", "nkx2.2a","fev","lmx1bb" , "sox19a", "dld", "dla", "dlb", "isl1", "foxa", "sim1a"))

saveRDS(Olig2_V3_inteagrated_obj, file = "~/Olig2_V3_inteagrated_obj.rds")
#
Olig2_V3_inteagrated_obj$Cell_labels <- "LFP"
Olig2_V3_inteagrated_obj$Cell_labels[Olig2_V3_inteagrated_obj$integrated_snn_res.2 %in% c(0, 1, 7, 8, 9, 11)] <- "LFP"
Olig2_V3_inteagrated_obj$Cell_labels[Olig2_V3_inteagrated_obj$integrated_snn_res.2 %in% c(4)] <- "MN"
Olig2_V3_inteagrated_obj$Cell_labels[Olig2_V3_inteagrated_obj$integrated_snn_res.2 %in% c(5)] <- "NPC"
Olig2_V3_inteagrated_obj$Cell_labels[Olig2_V3_inteagrated_obj$integrated_snn_res.2 %in% c(2)] <- "V3_pre"
Olig2_V3_inteagrated_obj$Cell_labels[Olig2_V3_inteagrated_obj$integrated_snn_res.2 %in% c(3, 6)] <- "V3"
Olig2_V3_inteagrated_obj$Cell_labels[Olig2_V3_inteagrated_obj$integrated_snn_res.2 %in% c(10)] <- "ISN_pre"
##
#

DimPlot(Olig2_V3_inteagrated_obj)
DimPlot(Olig2_V3_inteagrated_obj, group.by = "Cell_labels", label = T) + DimPlot(Olig2_V3_inteagrated_obj, group.by = "orig.ident", label = F)
FeaturePlot(Olig2_V3_inteagrated_obj, c("nkx2.2a", "nkx2.9"), split.by = "orig.ident")

## Publication plot:
pdf("~/Olig2_V3_integration_DimPlot.pdf", height = 10, width = 14)

DimPlot(Olig2_V3_inteagrated_obj, group.by = "Cell_labels", label = T, pt.size = 3, label.size = 10) + NoLegend() + theme(text = element_text(size = 30)) + ggtitle("") + scale_color_manual(values = c("LFP" = '#57C3F3', "NPCs" = "#53A85F","ISN_pre" = "#91D0BE", "ISN" = '#476D87', "MN" = '#E0D4CA', "V3_pre" = "#E59CC4", "V3" = "#E95C59")) 


dev.off()


library(gridExtra)
library(ggeasy)

nkx2.9_plot <- FeaturePlot(Olig2_V3_inteagrated_obj, "nkx2.9", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

fev_plot <- FeaturePlot(Olig2_V3_inteagrated_obj, "fev", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
sim1a_plot  <- FeaturePlot(Olig2_V3_inteagrated_obj, "sim1a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
sox1a_plot  <- FeaturePlot(Olig2_V3_inteagrated_obj, "sox1a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()


pdf("~/Olig2_V3_integration_FeaturePlot_nkx2.9_fev_sim1a_sox1a.pdf", width = 20* 2, height = 15*2)
grid.arrange(nkx2.9_plot, fev_plot, sim1a_plot, sox1a_plot, ncol=2, nrow=2)

dev.off()

##
nkx2.2a_plot <- FeaturePlot(Olig2_V3_inteagrated_obj, "nkx2.2a", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()

foxa_plot <- FeaturePlot(Olig2_V3_inteagrated_obj, "foxa", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
dld_plot  <- FeaturePlot(Olig2_V3_inteagrated_obj, "dld", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
dlb_plot  <- FeaturePlot(Olig2_V3_inteagrated_obj, "dlb", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
neurod4_plot  <- FeaturePlot(Olig2_V3_inteagrated_obj, "neurod4", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()
mnx1_plot  <- FeaturePlot(Olig2_V3_inteagrated_obj, "mnx1", pt.size = 3) + theme(text = element_text(size = 50), legend.key.height= unit(2, 'cm')) + easy_remove_axes()


pdf("~/Olig2_V3_integration_FeaturePlot_nkx2.2a_foxa_dld_dlb_neurod4_mnx1.pdf", width = 20* 3, height = 15*2)
grid.arrange(nkx2.2a_plot, foxa_plot, dld_plot, dlb_plot, neurod4_plot, mnx1_plot, ncol=3, nrow=2)

dev.off()

### Cell Label to find markers
Idents(Olig2_V3_inteagrated_obj) <- Olig2_V3_inteagrated_obj$Cell_labels

Olig2_V3_inteagrated_obj.all.mrkrs = FindAllMarkers(Olig2_V3_inteagrated_obj, only.pos = T)

Olig2_V3_integrated_CellType.mrkrs.list <- list()
#for(i in sort(as.character(unique(Cell_NP_obj.all.mrkrs$cluster)))){
for(i in as.character(unique(sort(Olig2_V3_inteagrated_obj$Cell_labels)))){
 
    Olig2_V3_integrated_CellType.mrkrs.list[[i]] <- Olig2_V3_inteagrated_obj.all.mrkrs[Olig2_V3_inteagrated_obj.all.mrkrs$cluster == i, ]

}



writexl::write_xlsx(Olig2_V3_integrated_CellType.mrkrs.list, "~/Olig2_V3_CellType_markers.xlsx")



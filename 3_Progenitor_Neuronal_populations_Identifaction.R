## neuron and progenitor neuron cells
# 1dpf
dpf1_obj <- readRDS("~/Rdata/dpf1_obj.rds")

q1 <- dpf1_obj

DimPlot(q1, label = T, group.by = "RNA_snn_res.1")

FeaturePlot(q1, c("nkx2.9", "tph2","nkx2.2a", "nkx2.2b", "ddc", "fev"))

Idents(q1) = q1$RNA_snn_res.1

q1$Integrated_markers <- "No"
q1$Integrated_markers[colnames(Sox1a_samples_sel_1dpf)] <- "Yes"
#
DimPlot(q1, label = T, group.by = "RNA_snn_res.1") + DimPlot(q1, label = T, group.by = "Integrated_markers")

## filter 1dpf cells
q1.idx = as.character(q1$RNA_snn_res.1) %in% c(0, 12, 1, 8, 7, 14, 19, 4, 17, 11, 10, 3, 9, 13)

DimPlot(q1, label = T, group.by = "RNA_snn_res.1") + DimPlot(q1, label = T, group.by = "Integrated_markers") + DimPlot(subset(q1, RNA_snn_res.1 %in% c(0, 12, 1, 8, 7, 14, 19, 4, 17, 11, 10, 3, 9, 13)), label = T, group.by = "RNA_snn_res.1")

q1_sel <- subset(q1, RNA_snn_res.1 %in% c(0, 12, 1, 8, 7, 14, 19, 4, 17, 11, 10, 3, 9, 13)) # 19, 4, 17, 11,

FeaturePlot(q1, c("sox2", "elavl3", "cldnb")) 
saveRDS(q1_sel, file = "~/Rdata/dpf1_obj_sel.rds")
##
# 2dpf:
q2 <- readRDS("~/Rdata/dpf2_obj.rds")

DimPlot(q2, label = T, group.by = "RNA_snn_res.1")

Idents(q2) = q2$RNA_snn_res.1
q2$Integrated_markers <- "No"
q2$Integrated_markers[colnames(Sox1a_samples_sel_2dpf)] <- "Yes"

DimPlot(q2, label = T, group.by = "RNA_snn_res.1") + DimPlot(q2, label = T, group.by = "Integrated_markers")

## filter 2dpf cells
q2.idx = as.character(q2$RNA_snn_res.1) %in% c(0, 7, 9, 10, 6)

DimPlot(q2, label = T, group.by = "RNA_snn_res.1") + DimPlot(q2, label = T, group.by = "Integrated_markers") + DimPlot(subset(q2, RNA_snn_res.1 %in% c(0, 7, 9, 10, 6)), label = T, group.by = "RNA_snn_res.1")

q2_sel <- subset(q2, RNA_snn_res.1 %in% c(0, 7, 9, 10, 6))
saveRDS(q2_sel, file = "~/Rdata/dpf2_obj_sel.rds")
## 3dpf:
q3 <- readRDS("~/Rdata/dpf3_obj.rds")

DimPlot(q3, label = T, group.by = "RNA_snn_res.1")

Idents(q3) = q3$RNA_snn_res.1
q3$Integrated_markers <- "No"
q3$Integrated_markers[colnames(Sox1a_samples_sel_3dpf)] <- "Yes"

DimPlot(q3, label = T, group.by = "RNA_snn_res.1") + DimPlot(q3, label = T, group.by = "Integrated_markers")

## filter 3dpf cells
q3.idx = as.character(q3$RNA_snn_res.1) %in% c(0, 5, 11, 14, 15)

DimPlot(q3, label = T, group.by = "RNA_snn_res.1") + DimPlot(q3, label = T, group.by = "Integrated_markers") + DimPlot(subset(q3, RNA_snn_res.1 %in% c(0, 5, 11, 14, 15)), label = T, group.by = "RNA_snn_res.1")

q3_sel <- subset(q3, RNA_snn_res.1 %in%  c(0, 5, 11, 14, 15))

saveRDS(q3_sel, file = "~/Rdata/dpf3_obj_sel.rds")

## 5dpf

q5 <- readRDS("~/Rdata/dpf5_obj.rds")

DimPlot(q5, label = T, group.by = "RNA_snn_res.1")

Idents(q5) = q5$RNA_snn_res.1
q5$Integrated_markers <- "No"
q5$Integrated_markers[colnames(Sox1a_samples_sel_5dpf)] <- "Yes"

DimPlot(q5, label = T, group.by = "RNA_snn_res.1") + DimPlot(q5, label = T, group.by = "Integrated_markers")

## filter 5dpf cells
q5.idx = as.character(q5$RNA_snn_res.1) %in% c(0, 7, 3, 11, 22, 26)

DimPlot(q5, label = T, group.by = "RNA_snn_res.1") + DimPlot(q5, label = T, group.by = "Integrated_markers") + DimPlot(subset(q5, RNA_snn_res.1 %in% c(0, 7, 3, 11, 22, 26)), label = T, group.by = "RNA_snn_res.1")

q5_sel <- subset(q5, RNA_snn_res.1 %in%  c(0, 7, 3, 11, 22, 26))

saveRDS(q5_sel, file = "~/Rdata/dpf5_obj_sel.rds")

### neuron and progenitor ratio
library(reshape2)
library(ggplot2)

## N / N + NP, NP / N + NP
df_p_n_rate <- data.frame(Time = c("1dpf", "2dpf", "3dpf", "5dpf"), Progenitors = c(0.687709, 0.1076426, 0.09393184, 0.0514246), Neurons = c(0.312291, 0.8923574, 0.9060682, 0.9485754) )
df_p_n_rate_melt <- melt(df_p_n_rate)
colnames(df_p_n_rate_melt) <- c("Time", "Cell_Type", "Ratio")
p1_ratio <- ggplot(data=df_p_n_rate_melt, aes(x=Time, y=Ratio, group=Cell_Type, color=Cell_Type)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+
  theme_minimal() + xlab("") + ylab("Ratio(Neuron + Progenitor cells)") +ylim(0, 1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pdf("~/Neuron_Progenitor_ratio2.pdf")
png("~/Neuron_Progenitor_ratio2.png")

print(p1_ratio)
dev.off()

#### now integrate them

new.counts = list(
    "dpf1" = q1@assays$RNA@counts[, q1.idx],
    "dpf2" = q2@assays$RNA@counts[, q2.idx],
    "dpf3" = q3@assays$RNA@counts[, q3.idx],
    "dpf5" = q5@assays$RNA@counts[, q5.idx] 
)

setwd("~/New_integration/")
{
    dir.create("Step1", showWarnings = F)
    setwd("Step1")
    dir.create("Step1multiRDAfiles", showWarnings = F)
    t1 = Sys.time()
    
    obj.lists = foreach(trt = names(new.counts)) %dopar% {
        obj.list = getObjs(new.counts[[trt]], regs.out = c("nCount_RNA"), resL = 1, npcas = 30, n.cells = 0, trt = trt)
        return(obj.list)
    }
    names(obj.lists) = names(new.counts)
    saveRDS(obj.lists, file = "Step1multiRDAfiles/obj.lists.step2.rds")
}

{
all_samples = names(obj.lists)
comparisons = list("all_samples" = all_samples)
    #
    setCPU(400/length(comparisons))
    the.results = foreach(comp = names(comparisons)) %dopar%{
      cat(comparisons[[comp]]," ", as.character(Sys.time()),"\n")
      aa = getIntegratedResults(the.obj.list = obj.lists, the_samples = comparisons[[comp]], vars_regress = c("nCount_RNA"))
      return(aa)
    }
    names(the.results) = names(comparisons)
    #saveRDS(the.results, file = paste0("Step1multiRDAfiles/the.results.rds"))
    saveRDS(the.results, file = paste0("Step1multiRDAfiles/Progenitor_Neuronal_population_integration.rds"))
}



q6 = the.results$all_samples
DimPlot(q6, label =T)

# marker filter cell in all cells
Sox1a_samples$filter_labels <- "No"
Sox1a_samples$filter_labels[colnames(q6)] <- "Yes"

DimPlot(Sox1a_samples, label =T, group.by = "filter_labels")


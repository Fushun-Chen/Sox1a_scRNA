{

setCPU = function(x) plan("multicore", workers = as.integer(x))

ldm = function(x) lapply(x, function(x) dim(x))

NL = Seurat::NoLegend()

getUniqs = function(xx){
  idx2 = duplicated(xx)
  for(i in seq(1, length(idx2),1)) if(idx2[i]) xx[i] = paste0(xx[i],"#",i)
  return(xx)
}

getMarkers = function(the.obj, the.idents, onlyPos = T){
  the.mrkrs = foreach(the.ident = names(table(the.idents))) %dopar% {
    print(the.ident)
    nw.idents = the.idents
    nw.idents[nw.idents != the.ident] = "rem.cells"
    Idents(the.obj) = nw.idents
    mrkrs = FindMarkers(the.obj, ident.1 = the.ident, ident.2 = "rem.cells", only.pos = onlyPos)
    mrkrs$cluster = the.ident
    mrkrs$gene = rownames(mrkrs)
    return(mrkrs)
  }
  names(the.mrkrs) = names(table(the.idents))
  
  tmp.mrkrs = the.mrkrs[[1]]
  tmp.mrkrs = tmp.mrkrs[order(tmp.mrkrs$avg_log2FC, decreasing = T),]
  for(i in seq(2, length(the.mrkrs))){
    tmp2 = the.mrkrs[[i]]
    tmp2 = tmp2[order(tmp2$avg_log2FC, decreasing = T),]
    tmp.mrkrs = rbind(tmp.mrkrs, tmp2)
  }
  tmp.mrkrs = tmp.mrkrs[order(tmp.mrkrs$avg_log2FC, decreasing = T),]
  return(tmp.mrkrs)
}


getTopN = function(theDfr, topN, clusterIds = NULL, sortIds = F){
  genes = c()
  the.idents = names(table(theDfr$cluster))
  if(sortIds){ 
    the.idents = as.numeric(the.idents)
    the.idents = the.idents[order(the.idents, decreasing = F)]
    the.idents = as.character(the.idents)
  }
  if(!is.null(clusterIds)) the.idents = the.idents[the.idents %in% clusterIds]
  for(idx in the.idents){
    aa = base::subset(theDfr, cluster == idx)
    genes = c(genes, aa$gene[1:topN])
  }
  genes = genes[!duplicated(genes)]
  return(genes)
}

# a common pipe for clustering
getObjs = function(the.data, regs.out = c("nCount_RNA"), resL = 0.5, npcas = 30, n.cells = 3, trt = NULL, topVar = 2000){
  t1 = Sys.time()
  #the.data = rawData[[trt]]
  raw.data = the.data 
  rwns = rownames(raw.data)
  mtgenes = grep("mt-", rwns, value = T)
  nUMI = colSums(raw.data)
  nMit = colSums(raw.data[rownames(raw.data) %in% mtgenes,])
  mito.perc = nMit/nUMI
  obj.list= CreateSeuratObject(raw.data, project = trt, min.cells = n.cells)
  obj.list$mito.perc = mito.perc
  obj.list = NormalizeData(obj.list)
  obj.list = FindVariableFeatures(obj.list, selection.method = "vst", nfeatures = topVar)
  all.genes <- rownames(obj.list)
  obj.list <- ScaleData(obj.list, features = all.genes, vars.to.regress = c("nCount_RNA")) #"mito.perc"
  obj.list <- RunPCA(obj.list, features = VariableFeatures(object = obj.list), npcs = npcas)
  
  #ElbowPlot(obj.list, ndims = npcas)
  obj.list <- FindNeighbors(obj.list, dims = 1:npcas)
  obj.list <- FindClusters(obj.list, resolution = resL)
  
  # If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
  obj.list <- RunUMAP(obj.list, dims = 1:npcas)#, umap.method = "umap-
  return(obj.list)
}


getIntegratedResults <- function(the.obj.list = NULL, the_samples = NULL, num_dims = 30, vars_regress = c("nCount_RNA"), resoLs = c(0.5, 1.0, 1.5, 2.0)){
  print("Find Anchors ...")
  obj.combined <- FindIntegrationAnchors(object.list = the.obj.list[the_samples], dims = 1:num_dims, verbose = FALSE)
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

#ldm(rawData)

getObjs = function(the.data, regs.out = c("nCount_RNA"), resL = 0.5, npcas = 30, n.cells = 3, trt = NULL, topVar = 2000){
  t1 = Sys.time()
  #the.data = rawData[[trt]]
  raw.data = the.data 
  rwns = rownames(raw.data)
  mtgenes = grep("mt-", rwns, value = T)
  nUMI = colSums(raw.data)
  nMit = colSums(raw.data[rownames(raw.data) %in% mtgenes,])
  mito.perc = nMit/nUMI
  obj.list= CreateSeuratObject(raw.data, project = trt)#, min.cells = n.cells)
  obj.list$mito.perc = mito.perc
  obj.list = NormalizeData(obj.list)
  obj.list = FindVariableFeatures(obj.list, selection.method = "vst", nfeatures = topVar)
  all.genes <- rownames(obj.list)
  obj.list <- ScaleData(obj.list, features = all.genes, vars.to.regress = c("nCount_RNA")) #"mito.perc"
  obj.list <- RunPCA(obj.list, features = VariableFeatures(object = obj.list), npcs = npcas)
  
  #ElbowPlot(obj.list, ndims = npcas)
  obj.list <- FindNeighbors(obj.list, dims = 1:npcas)
  obj.list <- FindClusters(obj.list, resolution = resL)
  
  # If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
  obj.list <- RunUMAP(obj.list, dims = 1:npcas)#, umap.method = "umap-
  return(obj.list)
}

{
  library("Seurat")
  library("Matrix")
  #library("tidyverse")
  library("foreach")
  library("doMC")
  registerDoMC(cores = 400)
  library("dplyr")
  library("future")
  plan("multicore", workers = 50)
  options(future.globals.maxSize = 16000*1024^2)
}  


readMatrix = function(thePath, trt){
    fls = dir(thePath)
    print(fls)
    print(thePath)
    ftr.fl = fls[!fls %in% c("matrix.mtx","barcodes.tsv")]
    print(ftr.fl)
    mtx = readMM(paste0(thePath,"/matrix.mtx"))
    brc = read.delim(paste0(thePath,"/barcodes.tsv"), header = F)
    brc = brc$V1
    brc = gsub("-1$","",brc)
    brc = gsub("_1$","",brc)
    brc = paste(trt, brc, sep = "_")
    print(head(brc))
    ftr = read.delim(paste0(thePath,"/",ftr.fl), header = F)
    #rwn = ftr$V2
    #idx = duplicated(rwn)
    #for(i in seq(1, length(idx),1)) if(idx[i]) rwn[i] = paste0(rwn[i],"#",i)
    #print(head(rwn))
    rownames(mtx) = getUniqs(ftr$V2) #rwn
    colnames(mtx) = brc
    return(mtx)
}

dir.create("20221121", showWarnings = F)
setwd("20221121")
}

{
if(!dir.exists("rdaFiles/lab3676_4174.rawData.rds")){
    dir.create("rdaFiles", showWarnings = F)
    print("Reading the data ..")
    rawData = foreach(trt = c("1dpf","2dpf","3dpf","5dpf")) %dopar% {
        #the.data = Seurat::Read10X(paste0("filtered/",trt,"_filtered_feature_bc_matrix"))
        the.data = readMatrix(paste0("../rawData/",trt,"_filtered_feature_bc_matrix"), trt)
        the.data = as.matrix(the.data)
        nUMI = colSums(the.data)
        nMit = colSums(the.data[grep("mt-",rownames(the.data)),])
        mito.perc = nMit/nUMI
        nGene = colSums(the.data > 0)
        #
        the.data = the.data[, mito.perc <= 0.20]
        nUMI = colSums(the.data)
        the.data = the.data[, nUMI >= 500]
        nGene = colSums(the.data > 0)
        the.data = the.data[, nGene >= 200]
        the.data = the.data[-grep("^ERCC-", rownames(the.data)),]
        #saveRDS(the.data, file = paste0("rdaFiles/",trt,".the.data.rds"))
        return(the.data)
    }
    names(rawData) = c("dpf1","dpf2","dpf3","dpf5")
    saveRDS(rawData, file = paste0("rdaFiles/lab3676_4174.rawData.rds"))
}else{
    print("Loading the data ...")
    rawData = paste0("rdaFiles/lab3676_4174.rawData.rds")
}


{
obj.lists = foreach(trt = names(rawData)) %dopar% {
    obj.list = getObjs(rawData[[trt]], regs.out = c("nCount_RNA"), resL = 1, npcas = 30, n.cells = 0, trt = trt)
    return(obj.list)
  }
  names(obj.lists) = names(rawData)
  saveRDS(obj.lists, file = "rdaFiles/lab3676_lab4174_obj.lists.all.dre.rds")
}


q1 = obj.lists$dpf1

for(resL in c(0.5,1,1.5,2,2.5,3)) q1 <- FindClusters(q1, resolution = resL)

q1.1 = subset(q1, RNA_snn_res.0.5 %in% seq(4,19))

q1.1 = subset(q1.1, RNA_snn_res.2 %in% as.character(q1.1$RNA_snn_res.2)[!as.character(q1.1$RNA_snn_res.2) %in% c(16,1,8,9,13,24,29,31)])

q1.1 = subset(q1.1, RNA_snn_res.1 %in% as.character(q1.1$RNA_snn_res.1)[!as.character(q1.1$RNA_snn_res.1) %in% c(0,2,3)])

Idents(q1.1) = as.character(q1.1$RNA_snn_res.2.5)

q1.1 = subset(q1.1, RNA_snn_res.2.5 %in% as.character(q1.1$RNA_snn_res.2.5)[!as.character(q1.1$RNA_snn_res.2.5) %in% c(19)])
Idents(q1.1) = as.character(q1.1$RNA_snn_res.3)

qq = q1.1@assays$RNA@scale.data
qq = subset(qq, rownames(qq) %in% grep("^hbbe1.1", rownames(qq), value = T))

nwidx = as.character(q1.1$RNA_snn_res.0.5)
names(nwidx) = colnames(q1.1)
qqq = subset(qq, rownames(qq) %in% "hbbe1.1")

nwidx[names(nwidx) %in% colnames(qqq)[qqq >= 0]] = "hbb"
nwidx[nwidx != "hbb"] = "other"
q1.1$hbb = nwidx
q1.1 = subset(q1.1, hbb %in% "other")

keep.cells = colnames(q1.1)
### analyze the data.
rawData$dpf1 = rawData$dpf1[,colnames(rawData$dpf1) %in% keep.cells]

# stp.data = list("obj.lists" = obj.lists, "rawData" = rawData, "keep.cells" = keep.cells, "dpf1.obj" = q1.1)
# saveRDS(stp.data, file = "rdaFiles/spt.data.rds")



{
    theGenes = readRDS("../../theGenes_all_ensVersions.rds")
    theGenes = lapply(theGenes, function(x){ x = x$ens2sym$dre; colnames(x) = c("dreENS","dreSYM"); return(x)})
    #
    tmp.data = readRDS("../../GSE173350/GSE173350_appel.int.RDS")
    tmp_counts = tmp.data@assays$RNA@counts
    rwns = rownames(tmp_counts)
    
    tmp.data = list("dpf24h" = tmp_counts[,tmp.data$orig.ident %in% "24h"],
    "dpf36h" = tmp_counts[,tmp.data$orig.ident %in% "36h"],
    "dpf48h" = tmp_counts[,tmp.data$orig.ident %in% "48hpf"])
    
    tmp.data = lapply(tmp.data, function(x) as.matrix(x))
    
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
    rem.genes = rem.genes[!(rem.genes$symb %in% new.rwns$symb),]
    
    new.rwns = new.rwns[!duplicated(new.rwns$dreENS),]
    
    new.ensids = read.delim("../rawData/1dpf_filtered_feature_bc_matrix/features.tsv", header = F)
    colnames(new.ensids) = c("dreENS","dreSYM","noo")
    new.rwns = new.rwns[new.rwns$dreENS %in% new.ensids$dreENS,]
    idx = match(new.rwns$dreENS, new.ensids$dreENS)
    new.rwns$newSYM = new.ensids$dreSYM[idx]
    new.rwns$symb = NULL
    
    idx = match(new.rwns$dreSYM, rownames(tmp_counts))
    new.counts = lapply(tmp.data, function(x) x[idx,])
    
    new.ensids = new.ensids[,c(1:2)]
    #not.dups = new.ensids$dreSYM 
    #print(head(not.dups))
    #idx = duplicated(not.dups)
    #for(i in seq(1,length(not.dups),1)){
    #  if(idx[i]) not.dups[i] = paste(not.dups[i],i, sep = "#")
    #}
    #new.ensids$not.dups = not.dups
    new.ensids$not.dups = getUniqs(new.ensids$dreSYM)
    idx = match(new.rwns$dreENS,new.ensids$dreENS)
    for(trt in names(new.counts)) rownames(new.counts[[trt]]) = new.ensids$not.dups[idx]
}


dre.rawData = c(rawData, new.counts)
rawData = dre.rawData
# new.ensids: use as common genes; to retrieve ensembl ids.
# dre.rawData: use for integration

#lapply(dre.rawData, function(x) dim(x))
ldm(dre.rawData)

for(trt in names(rawData)){
    tmp = rawData[[trt]]
    ngenes = rowSums(tmp > 0)
    tmp = tmp[ngenes >= 3,]
    rawData[[trt]] = tmp
}


setCPU(as.integer(400/length(rawData)))
###############################################################
{
    dir.create("Step1", showWarnings = F)
    setwd("Step1")
    dir.create("Step1multiRDAfiles", showWarnings = F)
    t1 = Sys.time()
    
    obj.lists = foreach(trt = names(rawData)) %dopar% {
        obj.list = getObjs(rawData[[trt]], regs.out = c("nCount_RNA"), resL = 1, npcas = 30, n.cells = 0, trt = trt)
        return(obj.list)
    }
    names(obj.lists) = names(dre.rawData)
    saveRDS(obj.lists, file = "Step1multiRDAfiles/obj.lists.step2.rds")
}

{
all_samples = names(obj.lists)
comparisons = list("all_samples" = all_samples,
	"lab3676" = c("dpf1","dpf5"),
	"lab4174" = c("dpf2","dpf3"),
	"sox1a123" = c("dpf1","dpf2","dpf3"), 
	"sox1aAll" = c("dpf1","dpf2","dpf3","dpf5"),
	"othrs1" = c("dpf24h","dpf36h","dpf48h"),
	"othrs2" = c("dpf1","dpf5","dpf24h","dpf36h","dpf48h")	
	)
    #
    setCPU(400/length(comparisons))
    the.results = foreach(comp = names(comparisons)) %dopar%{
      cat(comparisons[[comp]]," ", as.character(Sys.time()),"\n")
      aa = getIntegratedResults(the.obj.list = obj.lists, the_samples = comparisons[[comp]], vars_regress = c("nCount_RNA"))
      return(aa)
    }
    names(the.results) = names(comparisons)
    saveRDS(the.results, file = paste0("Step1multiRDAfiles/the.results.rds"))
}
}

###############################################################
setwd("Step1")
the.results = readRDS("Step1multiRDAfiles/the.results.rds")
obj.lists = readRDS("Step1multiRDAfiles/obj.lists.step2.rds")

##




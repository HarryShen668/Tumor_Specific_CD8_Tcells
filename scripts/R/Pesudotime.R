library(Seurat)
library(tidyverse)
library(ggsci)
source("/home/liuzipei/scripts/scRNAseq_2022/scripts/R/projection.plot.R")
source("~/scripts/utils_seurat.R")
ref <- readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref_harmony_symphony.Rds")
ref1 <-  readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref1_harmony_symphony.Rds")

rm(list = ls())
PanCancer <- readRDS("/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs/dataset3_GSE156728_PanCancer_CD8T.rds")
table(PanCancer$cancerType)

PanCancer$sampleID <- PanCancer$patient
PanCancer$Histology <- PanCancer$cancerType


tcr <- read.csv("/picb/lilab5/liuzipei/scRNAseq_2022/data/PanCancer_zhang/PanCancer_filtered_contig_annotations.csv") %>%
  filter(chain == "TRB") %>%
  select(barcode,cdr3,library.id,umis) %>%
  group_by(barcode,library.id) %>%
  slice_max(order_by = umis)

tcr$new.barcode <- paste(str_split(tcr$library.id, "-", simplify = T)[,1],tcr$barcode,sep = "_")
tcr <- tcr[!duplicated(tcr$new.barcode),]
table(duplicated(tcr$new.barcode))
PanCancer$new.barcode <- rownames(PanCancer@meta.data)
table(duplicated(PanCancer$new.barcode))
table(is.na(PanCancer$new.barcode)) 
unique(PanCancer$new.barcode)

meta <- left_join(PanCancer@meta.data,tcr)
rownames(meta) <- rownames(meta$new.barcode)
PanCancer@meta.data <- meta
rownames(PanCancer@meta.data) <- PanCancer@meta.data$new.barcode
predict_virus_cdr3b <- readRDS("~/project/scRNAseq_2022/predict_virus_cdr3b.rds")
table(PanCancer$cdr3 %in% predict_virus_cdr3b)

PanCancer$virus_predict <- ifelse(PanCancer$cdr3 %in% predict_virus_cdr3b, T, F)
table(PanCancer$virus_predict)
clonotype_count <- tibble(cdr3 = PanCancer$cdr3) %>% group_by(cdr3) %>% summarise(clonotype_count = n())
clonotype_count <- clonotype_count[nchar(clonotype_count$cdr3)>0, ]
clonotype_count
PanCancer@meta.data <- left_join(PanCancer@meta.data, clonotype_count) ## 注意给metadata重新赋值后会丢失行名
PanCancer@meta.data
rownames(PanCancer@meta.data) <- PanCancer$new.barcode
query.list <- SplitObject(PanCancer, split.by = "cancerType")
## load darklist
df <- read.csv("/picb/lilab5/liuzipei/Rdata/darklist.csv", header = T)
darklist_df <- as.list(df)
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)
ribo <- darklist_df$Ribosome
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)

for(sa in names(query.list)){
  query.list[[sa]]$pctMT <- PercentageFeatureSet(
    object = query.list[[sa]], pattern = "^MT-"
  )
}
for(sa in names(query.list)){
  query.list[[sa]]$pctRB <- PercentageFeatureSet(
    object = query.list[[sa]], pattern = "^RP[SL]"
  )
}
for(sa in names(query.list)){
  query.list[[sa]]$pctTCR <- PercentageFeatureSet(
    object = query.list[[sa]], pattern = "^TR[ABGD][VCJD]"
  )
}

for(sa in names(query.list)){
  query.list[[sa]]$pctFeature <- query.list[[sa]]$nCount_RNA/query.list[[sa]]$nFeature_RNA
}


for(sa in names(query.list)){
  query.list[[sa]] <- AddModuleScore(
    object = query.list[[sa]],
    features = list(DIG),
    name = 'DIG_SCORE',nbin=15)
}
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
for(sa in names(query.list)){
  query.list[[sa]] <- CellCycleScoring(query.list[[sa]],s.features = s.genes, g2m.features = g2m.genes,nbin=15)
  query.list[[sa]]$CC.Difference <- query.list[[sa]]$S.Score - query.list[[sa]]$G2M.Score
}
query.list.fil <- list()
for(i in 1:length(query.list)){
  na <- names(query.list)[i]
  query.list.fil[[na]] <- subset(query.list[[na]],
                                 subset = pctMT < 10 & pctRB < 45 & pctRB > 5 & pctTCR < 2 &
                                   pctTCR > 0 & nCount_RNA<20000 &nCount_RNA>1000& nFeature_RNA>800&
                                   nFeature_RNA<5000&pctFeature<6&pctFeature>1.6&DIG_SCORE1<15)
  
  
  query.list.fil[[na]] <- SCTransform(
    query.list.fil[[na]],
    assay = "RNA",
    new.assay.name = "SCT",
    variable.features.n = 3000,
    return.only.var.genes = FALSE,
    vars.to.regress = c("CC.Difference","DIG_SCORE1",'pctMT'), method = "glmGamPoi"
  )
}

query.data <- merge(query.list.fil[[1]], query.list.fil[2:length(query.list.fil)])
saveRDS(query.data, "/picb/lilab5/liuzipei/scRNAseq_2022/data/PanCancer20230424.rds")
query.data <- readRDS("/picb/lilab5/liuzipei/scRNAseq_2022/data/PanCancer20230424.rds")


marker <- readRDS("/picb/lilab5/liuzipei/Rdata/three_intersect.rds")
marker
query.data  <- AddModuleScore(query.data , features = list(marker), name = "TR_marker")

query.map <- mapQuery( 
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = query.map , cols = color18, linesize = 0.5, pointsize = 0.5)


ggplot(query.data@meta.data, aes(x = cancerType, y = TR_marker1, fill = cancerType)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('TR_marker') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  ylim(-0.5,1) 
## 依据boxplot大致趋势，选取score>0.2的作为TR候选
## 同时满足以下条件认为是TR
query.map$TR_predict <- ifelse(query.map$TR_marker1 > 0.25 & query.map$celltype_prob >= 0.8 & grepl("Tex", query.map$celltype) & !is.na(query.map$clonotype_count) & query.map$clonotype_count >= 2, TRUE, FALSE)
table(query.map$TR_predict)
unique(query.map$TR_predict)
table(query.map$virus_predict)
## 筛选TR克隆型及细胞
TR_clonotype_id <- query.map@meta.data[query.map$TR_predict,]$cdr3
TR_clonotype_id
TR_cells <- query.data[,query.data$cdr3 %in% TR_clonotype_id] ## TR cells
B_cells <- query.data[,query.data$virus_predict] ## Bystanders cells
  
## 可视化
TR.map <- mapQuery( 
  TR_cells@assays$SCT@scale.data,
  TR_cells@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' ) # return a Seurat object
TR.map <- knnPredict.Seurat(TR.map, ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = TR.map , cols = color18, linesize = 0.5, pointsize = 0.5, title = "TR")

B.map <- mapQuery( 
  B_cells@assays$SCT@scale.data,
  B_cells@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
B.map <- knnPredict.Seurat(B.map, ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = B.map , cols = color18, linesize = 0.5, pointsize = 0.5, title = "Bystander")

## reclustering
######################################################
# ## 高变基因 Top600 and Top2000
######################################################
scRNAlist_filter <- SplitObject(TR_cells, split.by = "cancerType")
for(na in seq(scRNAlist_filter)){
query.list.fil[[na]] <- SCTransform(
  query.list.fil[[na]],
  assay = "RNA",
  new.assay.name = "SCT",
  variable.features.n = 3000,
  return.only.var.genes = FALSE,
  vars.to.regress = c("CC.Difference","DIG_SCORE1",'pctMT'), method = "glmGamPoi")}

Top600_union <- NULL
for(i in seq(scRNAlist_filter)){
Top600_union <- union(x = Top600_union, y = VariableFeatures(scRNAlist_filter[[i]])[1:500])}
Top2000_intersect <- VariableFeatures(scRNAlist_filter[[1]])[1:3000]
for(i in 2:length(scRNAlist_filter)){
Top2000_intersect <- intersect(x = Top2000_intersect, y = VariableFeatures(scRNAlist_filter[[i]])[1:3000])}
var.features <- union(Top600_union, Top2000_intersect)

#var.features <- SelectIntegrationFeatures(object.list = scRNAlist_filter, nfeatures = 1500)
## filter darklist genes
df <- read.csv("/picb/lilab5/liuzipei/Rdata/darklist.csv", header = T)
darklist_df <- as.list(df)
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)
ribo <- darklist_df$Ribosome
MT <- darklist_df$Mitochondria
TCR <- grep(pattern = "^TR[ABDG][VDCJ]", x = var.features, value = TRUE)
#HLA <- grep(pattern = "^HLA", x = rownames(scRNAlist[[i]]), value = TRUE)
BCR <- grep(pattern = "^IGL[VC]", x = var.features, value = TRUE)
darklist <- c(TCR, DIG, ribo, BCR, MT, "MALAT1")
var.features <- setdiff(var.features, darklist)
myobj1 <- merge(scRNAlist_filter[[1]], scRNAlist_filter[2:length(scRNAlist_filter)])
VariableFeatures(myobj1) <- var.features
head(VariableFeatures(object = myobj1),20)

######################################################
# ## PCA and UMAP
######################################################
DefaultAssay(myobj1) <- "SCT"
myobj1 <- ScaleData(myobj1)
myobj1 <- RunPCA(object = myobj1, verbose = FALSE)
ElbowPlot(myobj1)
myobj1 <- RunUMAP(myobj1, reduction="pca", dims=1:30, verbose = F)
myobj1 <-   FindNeighbors(myobj1, reduction = "pca", dims = 1:30, verbose = F)
myobj1 <-   FindClusters(myobj1, resolution = 0.2, verbose = F)
DimPlot(myobj1, reduction = "umap") | DimPlot(myobj1, reduction = "umap", group.by = "meta.cluster") 
DimPlot(myobj1, reduction = "umap", group.by = "cancerType") | DimPlot(myobj1, reduction = "umap", group.by = "patient") 

######################################################
# ## Harmony and Clustering
######################################################
## Harmony
DefaultAssay(myobj1) <- "SCT"
library(harmony)
myobj1 <- RunHarmony(myobj1, group.by.vars = "patient", plot_convergence = TRUE, verbose = F)
myobj1 <- RunUMAP(myobj1, reduction="harmony", dims=1:15, verbose = F)
myobj1 <- FindNeighbors(myobj1, reduction = "harmony", dims = 1:15, verbose = F)
myobj1 <- FindClusters(myobj1, resolution = 0.12, verbose = F)
DimPlot(myobj1, reduction="umap", label=T)
myobj1 <- FindClusters(myobj1, resolution = c(seq(0.05,0.09,0.01),seq(0,1,0.1)), verbose = F)

myobj1 <- FindClusters(myobj1, resolution = c(seq(0.05,0.09,0.01),seq(0,1,0.1)), verbose = F)
clustree(myobj1, prefix = 'SCT_snn_res.') + coord_flip()
DimPlot(myobj1, reduction="umap", label=T)
myobj1 <- FindNeighbors(myobj1, reduction = "harmony", dims = 1:15, verbose = F)
myobj1 <- FindClusters(myobj1, resolution = 1.3, graph.name = "SCT_snn")
DimPlot(myobj1, reduction="umap", label=T)

clustree(myobj1, prefix = 'SCT_snn_res.') + coord_flip() 
DimPlot(myobj1, reduction="umap", group.by="SCT_snn_res.1.3", label=T) | DimPlot(myobj1, reduction="umap", group.by='meta.cluster', label=T, label.size = 2.8)
DimPlot(myobj1, reduction="umap", group.by="SCT_snn_res.0.6", label=T) | clustree(myobj1, prefix = 'SCT_snn_res.') + coord_flip() 
DimPlot(myobj1, reduction="umap", group.by="seurat_clusters", label=T)
DimPlot(myobj1, reduction="umap", group.by='meta.cluster', label=T, label.size = 2.8)
saveRDS(myobj1, file.path(PWD, "myobj1.rds"))

# pdf1 <- DimPlot(myobj1, reduction="umap", group.by = "orig.ident", label=T)
# pdf2 <- DimPlot(myobj1, reduction="umap", group.by='meta.cluster', label=T, label.size = 1.8) + NoLegend()
# ggsave(pdf1, dpi = 300, height = 4, width = 7.5,
#        file = file.path(PWD, 'pdf1.pdf'))
# ggsave(pdf2, dpi = 300, height = 4, width = 5.5,
#        file = file.path(PWD, 'pdf2.pdf'))
#myobj1 <- readRDS(file =  file.path(PWD, "myobj1.rds"))

pdf3 <- DimPlot(myobj1, group.by="seurat_clusters", reduction="umap", label=T) + NoLegend()
ggsave(pdf3, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'pdf3.pdf'))

######################################################
# ## DE analysis
######################################################
DefaultAssay(myobj1) <- "RNA"
myobj1 <- NormalizeData(myobj1, verbose = F)
myobj1 <- FindVariableFeatures(myobj1, selection.method = "vst", nfeatures = 2000, verbose = F)
all.genes <- rownames(myobj1)
myobj1 <- ScaleData(myobj1, feature = all.genes, verbose = F)
cluster_markers <- FindAllMarkers(myobj1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5_markers, file = file.path(PWD, 'top10_markers.csv'), row.names=TRUE, quote = FALSE)
DoHeatmap(myobj1, features = top5_markers$gene)
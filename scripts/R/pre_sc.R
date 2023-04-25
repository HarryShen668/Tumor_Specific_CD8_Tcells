########################################################################
# Prerare Date for scRNAseq
########################################################################

library(Seurat)
library(tidyverse)

ProjectDir <- "/picb/lilab5/liuzipei/Rdata/Integration_Harmony/" ## 设置项目目录
setwd(ProjectDir) 
DataDir <- "/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/GEX/"
SampleList <- list.files(DataDir)

seuList <- list()
seuList <- lapply(SampleList, function(sample){
  ## Create Object
  mtx <- Read10X(file.path(DataDir,paste0(sample,"/outs/filtered_feature_bc_matrix/")))
  seuobj <- CreateSeuratObject(mtx, project = "", min.cells = 10, min.features = 300)
  seuobj <- RenameCells(object = seuobj, new.names = paste0(sample,'_', colnames(seuobj)))
  
  ## Add metadata
  if(grep(sample) %in% c("KSA")){
    type <- "KSA"
  }
  else if(grep(sample) %in% c("QVD")){
    type <- "QVD"
  }
  else{
    tyep <- "PD1pos"
  }
  seuList[[sample]]$sampleID <- sample
  seuList[[sample]]$type <- type
})


data_list <- lapply(datasets, FUN = function(dataset){
  dataset_data <- readRDS(file.path(dir, dataset))
})
if(length(data_list) == 1){
  merge.data <- data_list[[1]]
}else{
  merge.data <- merge(data_list[[1]], y = data_list[2:length(data_list)])
}
scRNAlist <- SplitObject(merge.data, split.by = "patient")  ## 按patient拆分数据

######################################################
# ## QC
######################################################
## Add percent.mt 、percent.ribo and cellcyclescore to metadata
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
for(i in seq(scRNAlist)){
  scRNAlist[[i]]$percent.mt <- PercentageFeatureSet(object = scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]]$percent.ribo <- PercentageFeatureSet(object = scRNAlist[[i]], pattern = "^RP[LS]")
  scRNAlist[[i]] <- CellCycleScoring(scRNAlist[[i]], s.features = s.genes, g2m.features = g2m.genes, nbin=12)
  scRNAlist[[i]]$CC.Difference <- scRNAlist[[i]]$S.Score - scRNAlist[[i]]$G2M.Score}

## BasicQC plot
metaData.big <- NULL
vec <- c("cancerType","patient","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo")
for(i in seq(scRNAlist)){
  metaData.big <- rbind(metaData.big, scRNAlist[[i]]@meta.data[vec])}
show.nCount <- ggplot(metaData.big, aes(x = patient, y = nCount_RNA, fill = cancerType)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('nCount_RNA') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
show.nFeature <- ggplot(metaData.big, aes(x = patient, y = nFeature_RNA, fill = cancerType)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('nFeature_RNA') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
show.percent.mt <- ggplot(metaData.big, aes(x = patient, y = percent.mt, fill = cancerType)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('percent.mt') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
show.percent.ribo <- ggplot(metaData.big, aes(x = patient, y = percent.ribo, fill = cancerType)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('percent.ribo') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
basicQC <- show.nCount + show.nFeature + show.percent.mt + show.percent.ribo
basicQC
# ggsave(basicQC, dpi = 300, height = 8, width = 9,
#        file = file.path(PWD, 'basicQC.pdf'))

######################################################
# ## Subset + SCTransform
######################################################
scRNAlist_filter <- lapply(scRNAlist, FUN = function(sample){
  sample_data <- subset(sample, subset = percent.mt < 7 & nCount_RNA < 40000 & nCount_RNA > 2000)
  sample_data <- SCTransform(
    sample,
    assay = "RNA",
    new.assay.name = "SCT",
    variable.features.n = 3000,
    method="glmGamPoi",
    return.only.var.genes = FALSE,
    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
    verbose = FALSE)
})

######################################################
# ## 高变基因 Top600 and Top2000
######################################################
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



# ######################################################
# # ## Clustering per dataset
# ######################################################
# for(i in seq(scRNAlist_filter)){
#   scRNAlist_filter[[i]] <- RunPCA(scRNAlist_filter[[i]], verbose = FALSE)
#   scRNAlist_filter[[i]] <- RunUMAP(scRNAlist_filter[[i]], dims = 1:30, verbose = FALSE)
#   scRNAlist_filter[[i]] <- FindNeighbors(scRNAlist_filter[[i]], dims = 1:15, k.param = 10)
#   scRNAlist_filter[[i]] <- FindClusters(scRNAlist_filter[[i]], resolution = 50)
# }
# saveRDS(scRNAlist_filter, file = file.path(PWD,"scRNAlist_filter_final_zzm.rds"))
# 
# #########################################################
# # ## Minicluster
# #########################################################
# Minicluster <- NULL
# for(i in seq(scRNAlist_filter)){
#   scRNAlist_filter[[i]] <- NormalizeData(object = scRNAlist_filter[[i]], assay = "RNA")
#   Minicluster[[i]] <- AverageExpression(scRNAlist_filter[[i]], return.seurat = TRUE)
#   Minicluster[[i]]$dataset <- as.factor(i)}
# 
# 
# Intersect_genes <- rownames(Minicluster[[1]])
# for(i in 2:length(Minicluster)){
#   Intersect_genes <- intersect(x = Intersect_genes, y = rownames(Minicluster[[i]]))}
# length(Intersect_genes)
# myobj1 <- merge(x = Minicluster[[1]], y = Minicluster[2:length(Minicluster)])
# ## 只保留所有数据集均包含的基因
# myobj1 <- subset(myobj1, features = Intersect_genes)
# Idents(myobj1) <- myobj1$dataset

######################################################
# ## PCA and UMAP
######################################################
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
myobj1 <- readRDS("../Integration_CCA/myobj1_.rds")

DefaultAssay(myobj1) <- "SCT"
library(harmony)
myobj1 <- RunHarmony(myobj1, group.by.vars = "patient", plot_convergence = TRUE, verbose = F)
myobj1 <- RunUMAP(myobj1, reduction="harmony", dims=1:15, verbose = F)
myobj1 <- FindNeighbors(myobj1, reduction = "harmony", dims = 1:15, verbose = F)
myobj1 <- FindClusters(myobj1, resolution = c(seq(0.05,0.09,0.01),seq(0,1,0.1)), verbose = F)
clustree(myobj1, prefix = 'SCT_snn_res.') + coord_flip()


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

###########################################################
# ## 可视化
###########################################################
EX_mrk <- c('PDCD1',"CXCL13",'ENTPD1','TIGIT','LAG3','CTLA4','TOX','HAVCR2','LAYN')
Naive_mrk <-c('CCR7', 'LEF1', 'SELL', "TCF7")
EMRA_mrk <- c("KLRG1","CX3CR1","FCGR3A","FGFBP2")
MAIT_mrk <- c("SLC4A10","KLPB1","ZBTB16","NCR3","RORC")
ISG_mrk <- c("IFIT1")
K_mrk <- c("TYROBP","KIR2DL4")
EM_mrk <- c("GZMK")
RM_mrk <- c("ZNF683","ITGAE")
M_mrk <- c("IL7R")
new_obj <- myobj1
FeaturePlot(new_obj, features = EX_mrk)
FeaturePlot(new_obj, features = Naive_mrk)
FeaturePlot(new_obj, features = EMRA_mrk)
FeaturePlot(new_obj, features = MAIT_mrk)
FeaturePlot(new_obj, features = ISG_mrk)
FeaturePlot(new_obj, features = K_mrk)
FeaturePlot(new_obj, features = EM_mrk)
FeaturePlot(new_obj, features = RM_mrk)
FeaturePlot(new_obj, features = M_mrk)


myobj1 <- AddModuleScore(
  object = myobj1,
  features = list(EX_mrk),
  name = 'EX_sig')
FeaturePlot(myobj1, features = EX_mrk)
VlnPlot(myobj1, features = EX_mrk, pt.size = 0)
RidgePlot(myobj1, features = EX_mrk)
DotPlot(myobj1, features = EX_mrk) + RotatedAxis()
DotPlot.pdf <- DotPlot(myobj1, features = EX_mrk) + RotatedAxis()
ggsave(DotPlot.pdf, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'DotPlot.pdf'))
FeaturePlot.pdf <- FeaturePlot(myobj1, features = "EX_sig1")
ggsave(FeaturePlot.pdf, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'FeaturePlot.pdf'))
VlnPlot(myobj1, features = "EX_sig1", pt.size = 0)
VlnPlot.pdf <- VlnPlot(myobj1, features = "EX_sig1", pt.size = 0)
ggsave(VlnPlot.pdf, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'VlnPlot.pdf'))
RidgePlot(myobj1, features = "EX_sig1")
DotPlot(myobj1, features = "EX_sig1") + RotatedAxis()
save(list = ls(all.names = TRUE), file = file.path(PWD,"Harmony_Integrtion.Rdata"))




#####################################################
## Annotation by hand
#####################################################
FeaturePlot(myobj1, features = c("MAL","IL7R","RPS12","CD52","CXCR5","GZMK","CX3CR1","TYROBP",
                                 "KIR2DL4","ZNF683","PDCD1","CXCL13","myl12a","TCF7","IFIT1","SLC4A10","NME1"))
new.cluster.ids <- c("Tm.IL7R","Tem.GZMK","Tm.IL7R","Tex.CXCL13","Temra.FGFBP2","Trm.LDCRAD4","Tm","Tm","Tem.GZMK",
                     "Tn.CCR7","Tex.CXCL13","HSP.HSPA6","ISG.IFI6","Tm","Prol.STMN1","Temra.FGFBP2","MAIT.KLRB1", "Tk.TYROBP",
                     "Temra.FGFBP2","Prol.STMN1")
names(new.cluster.ids) <- levels(myobj1)
myobj1 <- RenameIdents(myobj1, new.cluster.ids)
DimPlot(myobj1, reduction = "umap", label = T,
        cols = color20 <- c('#d62e2d', '#e9787a', '#cd9c9b', '#684797', '#3377a9', '#96c3d8',
                            '#67a59b', '#70b866', '#6a9a52', '#a5d38f', '#ab9697', '#f19294',
                            '#f5b375', '#da8f6f', '#e0c880', '#f47d2f', '#ec8d63', '#e45d61',
                            '#a65a35', '#4b9d47'))
DimPlot(myobj1, reduction = "umap", group.by = "meta.cluster",
        cols = color20 <- c('#d62e2d', '#e9787a', '#cd9c9b', '#684797', '#3377a9', '#96c3d8',
                            '#67a59b', '#70b866', '#6a9a52', '#a5d38f', '#ab9697', '#f19294',
                            '#f5b375', '#da8f6f', '#e0c880', '#f47d2f', '#ec8d63', '#e45d61',
                            '#a65a35', '#4b9d47'))
saveRDS(myobj1, file = file.path(PWD,"myobj_final2.rds"))
myobj1 <- readRDS(file = file.path(PWD,"myobj_final2.rds"))
DimPlot(myobj1)
myobj1$meta.cluster
DimPlot(myobj1, group.by = "meta.cluster")
unique(myobj1$seurat_clusters)
new_obj <- myobj1[,!(myobj1$seurat_clusters %in% c("6","7"))]
DimPlot(new_obj, reduction = "umap")
DimPlot(new_obj, reduction = "umap", group.by = "meta.cluster")
saveRDS(new_obj, file = file.path(PWD, "myobj_final_filter.rds"))

devtools::install_github("sajuukLyu/ggunchull", type = "source")
devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
options(repr.plot.width=13,repr.plot.height=10)
clusterCornerAxes(object = seuItg.dft,
                  reduction = 'umap', 
                  noSplit = F,clusterCol = 'SCT_snn_res.1.2',
                  #                  groupFacet = 'diseaseStatus',
                  aspect.ratio = 1,
                  relLength = 0.5,
                  axes = 'mul',cellLabel = T,cellLabelSize = 5)+scale_fill_npg()

# DimPlot(seuItg.dft, reduction = "umap", label=T, label.size = 2.5) + NoLegend()
# pdf4 <- DimPlot(seuItg.dft, reduction = "umap", label=T, label.size = 2.5) + NoLegend()
# ggsave(pdf4, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'pdf4.pdf'))

#####################################################
## ROGUE
#####################################################
#devtools::install_github("PaulingLiu/ROGUE")
library(ROGUE)
library(tidyverse)
library(ggplot2)

expr <- myobj1@assays$RNA@counts
expr <- as.data.frame(expr)
expr[1:5, 1:4]
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)
ent.res <- SE_fun(expr)
head(ent.res)
SEplot(ent.res)
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value
rogue.res <- rogue(expr, labels = myobj1$meta.cluster, samples = myobj1$orig.ident, platform = "UMI", span = 0.6)
rogue.res
rogue.boxplot(rogue.res)
rogue.boxplot <- rogue.boxplot(rogue.res)
ggsave(rogue.boxplot, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'rogue.boxplot.pdf'))

######################################################
# ## singleR annotation
######################################################
# library(SingleR)
# #ref <- ImmGenData()
# #ref <- DatabaseImmuneCellExpressionData()
# ref <- MonacoImmuneData()
# anno.cluster <- SingleR(test = myobj1@assays$RNA@data, ref = ref,
#                         labels = ref$label.fine, clusters = myobj1@active.ident)
# table(anno.cluster$labels)
# # celltype <- data.frame(seuItg.dft$seurat_clusters, anno.cluster$labels)
# # table(celltype[,1:2])


# ######################################################
# # ## scGSEA Example
# ######################################################
# #https://www.jianshu.com/p/00f069cb239c
# # library(devtools)
# # install_github('immunogenomics/presto')
# library(presto)
# # install.packages("msigdbr")
# library(msigdbr)
# library(fgsea)
# library(dplyr)
# library(ggplot2)
# library(tidyverse)
# DefaultAssay(seuItg.dft) <- "RNA"
# de.genes <- wilcoxauc(seuItg.dft,'seurat_clusters')
# head(de.genes)
# dplyr::count(de.genes, group)
# ## 选择group0进行GSEA分析
# cluster0.genes <- de.genes %>% filter(group==0) %>% arrange(desc(auc)) %>% select(feature, auc)
# cluster0.genes <- deframe(cluster0.genes)
# head(cluster0.genes)
# ## MsigDB中的C7免疫基因集
# gene_set <- msigdbr(species = "Homo sapiens", category = "C7") ##选择GSEA的目标基因集
# head(gene_set)
# fgsea_sets <- gene_set %>% split(x=.$gene_symbol, f=.$gs_name) ## 生成参考目标基因集list
# fgsea_sets
# fgseaRes <- fgsea(fgsea_sets, stats = cluster0.genes, nperm = 1000)
# #fgseaRes <- fgseaMultilevel(fgsea_sets, stats = cluster0.genes, scoreType = "pos")
# fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
# fgseaResTidy %>% dplyr::select(-leadingEdge,-ES,-nMoreExtreme)  %>% arrange(padj) %>% head()
# 
# ggplot(fgseaResTidy %>% filter(padj<0.008) %>% head(n=20),aes(reorder(pathway,NES),NES)) +
#   geom_col(aes(fill= NES < 7.5)) +
#   coord_flip() +
#   labs(x="Pathway",y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA") +
#   theme_minimal()
# plotEnrichment(fgsea_sets[["NAKAYA_PBMC_FLUMIST_AGE_18_50YO_7DY_DN"]],ranks) + labs(title = "NAKAYA PBMC FLUMIST AGE 18 50YO 7DY DN")

######################################################
# ## scGSEA Analysis
######################################################
# library(ggplot2)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# DefaultAssay(seuItg.dft) <- "RNA"
# cluster_markers <- FindAllMarkers(seuItg.dft, only.pos = TRUE,min.pct=0.25,logfc.threshold = 0.25)
# geneList <- cluster_markers$avg_log2FC
# rownames(cluster_markers)
# names(geneList) <- rownames(cluster_markers)
# geneList <- sort(geneList,decreasing = T)
# head(geneList)

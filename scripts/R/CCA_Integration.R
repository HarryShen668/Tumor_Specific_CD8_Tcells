########################################################################
# Created on 2023-03-05
# Use SCTransform and rpca For scRNAseq Analysis
########################################################################

library(gghalves)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggsci)
library(patchwork)
library(clustree)
rm(list=ls())
setwd("/picb/lilab5/liuzipei/Rdata/Integration_CCA/") ## 设置当前工作目录
PWD = getwd()

######################################################
# ## Prepare data 
######################################################
dir <- '/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs'
datasets = c(
#  'dataset1_GSE176021_NSCLC_CD8T.rds',
#  'dataset2_GSE180268_HNSCC_CD8T.rds',
  'dataset3_GSE156728_PanCancer_CD8T.rds')
#  'dataset4_Rosenberg_allT.rds')
data_list <- lapply(datasets, FUN = function(dataset){
  dataset_data <- readRDS(file.path(dir, dataset))
})
if(length(data_list) == 1){
  merge.data <- data_list[[1]]
}else{
  merge.data <- merge(data_list[[1]], y = data_list[2:length(data_list)])
}
table(merge.data$patient)
filtered_p <- names(table(merge.data$patient)[table(merge.data$patient) > 1000])
merge.data <- subset(merge.data, patient %in% filtered_p)
table(merge.data$patient)
merge.data <- subset(merge.data, patient != "OV.P20190304")
table(merge.data$cancerType)
scRNAlist <- SplitObject(merge.data, split.by = "cancerType") 

# ######################################################
# # ## QC
# ######################################################
# ## Add percent.mt 、percent.ribo and cellcyclescore to metadata
df <- read.csv("/picb/lilab5/liuzipei/Rdata/darklist.csv", header = T)
darklist_df <- as.list(df)
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)
for(i in seq(scRNAlist)){
  s.genes <- intersect(cc.genes$s.genes,rownames(scRNAlist[[i]]))
  g2m.genes <- intersect(cc.genes$g2m.genes,rownames(scRNAlist[[i]]))
  scRNAlist[[i]]$percent.mt <- PercentageFeatureSet(object = scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]]$percent.ribo <- PercentageFeatureSet(object = scRNAlist[[i]], pattern = "^RP[LS]")
  scRNAlist[[i]] <- CellCycleScoring(scRNAlist[[i]], s.features = s.genes, g2m.features = g2m.genes)
  scRNAlist[[i]]$CC.Difference <- scRNAlist[[i]]$S.Score - scRNAlist[[i]]$G2M.Score
  scRNAlist[[i]] <- AddModuleScore(
    object = scRNAlist[[i]],
    features = list(DIG),
    name = 'DIG_sig')}

## BasicQC plot
metaData.big <- NULL
vec <- c("dataset","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo")
for(i in seq(scRNAlist)){
  metaData.big <- rbind(metaData.big, scRNAlist[[i]]@meta.data[vec])}

show.nCount <- ggplot(metaData.big, aes(x = dataset, y = nCount_RNA, fill = dataset)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('nCount_RNA') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

show.nFeature <- ggplot(metaData.big, aes(x = dataset, y = nFeature_RNA, fill = dataset)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('nFeature_RNA') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

show.percent.mt <- ggplot(metaData.big, aes(x = dataset, y = percent.mt, fill = dataset)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('percent.mt') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

show.percent.ribo <- ggplot(metaData.big, aes(x = dataset, y = percent.ribo, fill = dataset)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('percent.ribo') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

basicQC <- show.nCount + show.nFeature + show.percent.mt + show.percent.ribo
basicQC
ggsave(basicQC, dpi = 300, height = 8, width = 9,
       file = file.path(PWD, 'basicQC.pdf'))

# ######################################################
# # ## Filtering
# ######################################################
scRNAlist_filter <- lapply(scRNAlist, FUN = function(sample){
  sample_data <- subset(sample, subset = percent.mt < 7 & nCount_RNA < 40000 & nCount_RNA > 2000)
  sample_data <- SCTransform(
    sample,
    assay = "RNA",
    new.assay.name = "SCT",
    variable.features.n = 3000,
    method="glmGamPoi",
    return.only.var.genes = FALSE,
    vars.to.regress = c("percent.mt","CC.Difference","patient"),
    verbose = FALSE)
})
saveRDS(scRNAlist_filter, file.path(PWD, "scRNAlist_filter_new.rds"))
# ######################################################
# # ## Top600 and Top2000
# ######################################################
## 取大概1500个高变基因。 Top800取并集,Top2500取交集
# Top600_union <- NULL
# for(i in seq(scRNAlist_filter)){
#   Top600_union <- union(x = Top600_union, y = VariableFeatures(scRNAlist_filter[[i]])[1:600])
# }
# Top2000_intersect <- VariableFeatures(scRNAlist_filter[[1]])[1:2000]
# for(i in 2:length(scRNAlist_filter)){
#   Top2000_intersect <- intersect(x = Top2000_intersect, y = VariableFeatures(scRNAlist_filter[[i]])[1:2000])
# }
# red_features <- union(Top600_union, Top2000_intersect)
# red_features <- setdiff(red_features, darklist)
var.features <- SelectIntegrationFeatures(object.list = scRNAlist_filter, nfeatures = 2000)
df <- read.csv("/picb/lilab5/liuzipei/Rdata/darklist.csv", header = T)
darklist_df <- as.list(df)
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)
ribo <- darklist_df$Ribosome
MT <- darklist_df$Mitochondria
TCR <- grep(pattern = "^TR[ABDG][VDCJ]", x = var.features, value = TRUE)
#HLA <- grep(pattern = "^HLA", x = rownames(scRNAlist[[i]]), value = TRUE)
BCR <- grep(pattern = "^IGL[VC]", x = var.features, value = TRUE)
darklist <- c(TCR, DIG, ribo, BCR, MT, "MALAT1")

var.features <- readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/red_features.Rds")
var.features <- setdiff(var.features, darklist)

######################################################
# ## Integration 
######################################################
SCTanchor_list<-list()
for(i in seq(scRNAlist_filter)){
  SCTanchor_list[[i]] <- scRNAlist_filter[[i]]@assays$SCT@SCTModel.list$model1@feature.attributes}


SCTanchor2 <- Reduce(intersect, x = lapply(X = SCTanchor_list, FUN = rownames))
SCTanchor2 <- setdiff(SCTanchor2, darklist)
scRNAlist_filter <- PrepSCTIntegration(object.list = scRNAlist_filter, anchor.features = var.features)
scRNAlist_filter <- lapply(X = scRNAlist_filter, FUN = function(x){x <- RunPCA(x, verbose=FALSE)})
Anchors.default <- FindIntegrationAnchors(
  object.list = scRNAlist_filter,
  reference=c(1,2,3),
  reduction="rpca",
  dims=1:50,
  normalization.method = 'SCT',
  anchor.features = var.features,
  k.anchor = 5,
  k.filter = 300,
  k.score = 30)
seuItg.dft <- IntegrateData(
  anchorset = Anchors.default,
  normalization.method = 'SCT', features.to.integrate=SCTanchor2)

saveRDS(seuItg.dft, file.path(PWD, "seuItg.dft.rds"))
#seuItg.dft <- readRDS(file =  file.path(PWD, "seuItg.dft.rds"))

######################################################
# ## Reduce and Clustree 
######################################################
DefaultAssay(seuItg.dft) <- "integrated"
red_features <- var.features
head(red_features, 50)
seuItg.dft <- ScaleData(seuItg.dft, features = rownames(seuItg.dft))
seuItg.dft <- RunPCA(seuItg.dft, features = red_features ,npcs = 30, verbose = FALSE)
seuItg.dft <- RunUMAP(seuItg.dft, reduction = "pca", dims = 1:30)
DimPlot(seuItg.dft, reduction="umap", label=T)
DimPlot(seuItg.dft, reduction="umap", group.by='meta.cluster', label=T, label.size = 3) + NoLegend()
Tzhang <- subset(seuItg.dft, dataset == "3")
DimPlot(Tzhang, reduction="umap", group.by='meta.cluster', label=T, label.size = 3) + NoLegend()
Tex_list = list(Tex=colnames(seuItg.dft)[seuItg.dft$meta.cluster %in% c("CD8.c11.Tex.PDCD1","CD8.c14.Tex.TCF7","CD8.c12.Tex.CXCL13","CD8.c13.Tex.myl12a")])
DimPlot(seuItg.dft, cells.highlight = list(Tex=colnames(seuItg.dft)[seuItg.dft$meta.cluster %in% c("CD8.c11.Tex.PDCD1","CD8.c14.Tex.TCF7","CD8.c12.Tex.CXCL13","CD8.c13.Tex.myl12a")]))
pdf1 <- DimPlot(seuItg.dft, reduction="umap", label=T)
pdf2 <- DimPlot(Tzhang, reduction="umap", group.by='meta.cluster', label=T, label.size = 1.0) + NoLegend()
DimPlot(seuItg.dft, cells.highlight = Tex_list, sizes.highlight = 0.5) + NoLegend()
pdf_Tex <- DimPlot(seuItg.dft, cells.highlight = Tex_list, sizes.highlight = 0.5) + NoLegend()
ggsave(pdf1, dpi = 300, height = 4, width = 7.5, file = file.path(PWD, 'pdf1.pdf'))
ggsave(pdf2, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'pdf2.pdf'))
ggsave(pdf_Tex, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'pdf_Tex.pdf'))
#seuItg.dft <- RunTSNE(seuItg.dft, reduction = "pca", dims = 1:30)
#DimPlot(seuItg.dft, reduction="tsne", label=T)
#DimPlot(seuItg.dft, reduction="tsne", group.by='meta.cluster', label=T)
seuItg.dft <- FindNeighbors(seuItg.dft, reduction = "pca", dims = 1:30)
seuItg.dft <- FindClusters(seuItg.dft, resolution = c(seq(0.05,0.09,0.01),seq(0,1,0.1)))
clustree(seuItg.dft, prefix = 'integrated_snn_res.') + coord_flip()
#saveRDS(seuItg.dft, file = file.path(dirPjtHome, 'seuItg.dft_dimentionred.Rds'))

######################################################
# ## FindCluster
######################################################
DefaultAssay(seuItg.dft) <- "integrated"
seuItg.dft <- FindNeighbors(object=seuItg.dft, dims=1:15)
seuItg.dft <- FindClusters(object=seuItg.dft, resolution = 0.5)
DimPlot(seuItg.dft, group.by="seurat_clusters", reduction="umap", label=T) + NoLegend()
DimPlot(seuItg.dft, group.by="seurat_clusters", reduction="umap", label=T, cells.highlight = Tex_list) + NoLegend()
DimPlot(seuItg.dft, group.by="integrated_snn_res.0.5", reduction="umap", label=T) + NoLegend()
DimPlot(seuItg.dft, group.by="meta.cluster", reduction="umap", label=T) + NoLegend()
#DimPlot(seuItg.dft, group.by="seurat_clusters", reduction="tsne", label=T) + NoLegend()
pdf3 <- DimPlot(seuItg.dft, group.by="seurat_clusters", reduction="umap", label=T) + NoLegend()
ggsave(pdf3, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'pdf3.pdf'))

###########################################################
# ## Annotation
###########################################################
EX_mrk<-c('PDCD1',"CXCL13",'ENTPD1','TIGIT','LAG3','CTLA4','TOX','HAVCR2')
seuItg.dft <- AddModuleScore(
  object = seuItg.dft,
  features = list(EX_mrk),
  name = 'EX_sig')
FeaturePlot(seuItg.dft, features = EX_mrk)
VlnPlot(seuItg.dft, features = EX_mrk, pt.size = 0)
RidgePlot(seuItg.dft, features = EX_mrk)
DotPlot(seuItg.dft, features = EX_mrk) + RotatedAxis()
DotPlot.pdf <- DotPlot(seuItg.dft, features = EX_mrk) + RotatedAxis()
ggsave(DotPlot.pdf, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'DotPlot.pdf'))

FeaturePlot(seuItg.dft, features = "EX_sig1")
VlnPlot(seuItg.dft, features = "EX_sig1", pt.size = 0)
VlnPlot.pdf <- VlnPlot(seuItg.dft, features = "EX_sig1", pt.size = 0)
ggsave(VlnPlot.pdf, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'VlnPlot.pdf'))
RidgePlot(seuItg.dft, features = "EX_sig1")
DotPlot(seuItg.dft, features = "EX_sig1") + RotatedAxis()

save(list = ls(all.names = TRUE), file = file.path(PWD,"CCA_Integrtion.Rdata"))
#saveRDS(seuItg.dft, file.path(PWD, "seuItg.dft.final.rds"))
#seuItg.dft <- readRDS(file.path(PWD, "seuItg.dft.final.rds"))

######################################################
# ## ROGUE
######################################################
#devtools::install_github("PaulingLiu/ROGUE")
# library(ROGUE)
# library(tidyverse)
# library(ggplot2)
# 
# expr <- seuItg.dft@assays$RNA@counts
# expr <- as.data.frame(expr)
# expr[1:5, 1:4]
# expr <- matr.filter(expr, min.cells = 10, min.genes = 10)
# ent.res <- SE_fun(expr)
# head(ent.res)
# SEplot(ent.res)
# rogue.value <- CalculateRogue(ent.res, platform = "UMI")
# rogue.value
# rogue.res <- rogue(expr, labels = seuItg.dft$seurat_clusters, samples = seuItg.dft$orig.ident, platform = "UMI", span = 0.6)
# rogue.res
# rogue.boxplot <- rogue.boxplot(rogue.res)
# ggsave(rogue.boxplot, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'rogue.boxplot.pdf'))

######################################################
# ## Annotation by hand
######################################################
# FeaturePlot(seuItg.dft, features = c("MAL","IL7R","RPS12","CD52","CXCR5","GZMK","CX3CR1","TYROBP",
#                                      "KIR2DL4","ZNF683","PDCD1","CXCL13","myl12a","TCF7","IFIT1","SLC4A10","NME1"))
# new.cluster.ids <- c("0.Tm.IL7R+ZNF683", "1.Tem.CXCR5+GZMK", "2.Tm.NME1+RPS12+CD52", "3.Tex.CXCL13+PDCD1+TCF7+myl12a",
#                      "4.Tn.MAL", "5.Temra.CX3CR1", "6.Temra.CX3CR1+Tk.TYROBP", "7.Tk.KIR2DL4+Tm.IL7R", "8.MAIT.SLC4A10", "9.NULL", "10.ISG.IFIT1")
# names(new.cluster.ids) <- levels(seuItg.dft)
# seuItg.dft <- RenameIdents(seuItg.dft, new.cluster.ids)
# DimPlot(seuItg.dft, reduction = "umap", label=T, label.size = 2.5) + NoLegend()
# pdf4 <- DimPlot(seuItg.dft, reduction = "umap", label=T, label.size = 2.5) + NoLegend()
# ggsave(pdf4, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'pdf4.pdf'))

######################################################
# ## singleR annotation
######################################################
library(SingleR)
#ref <- ImmGenData()
#ref <- DatabaseImmuneCellExpressionData()
# ref <- MonacoImmuneData()
# anno.cluster <- SingleR(test = seuItg.dft@assays$RNA@data, ref = ref,
#                         labels = ref$label.fine, clusters = seuItg.dft@active.ident)
# table(anno.cluster$labels)
# celltype <- data.frame(seuItg.dft$seurat_clusters, anno.cluster$labels)
# table(celltype[,1:2])

######################################################
# ## DE analysis
######################################################
# DefaultAssay(seuItg.dft) <- "RNA"
# cluster_markers <- FindAllMarkers(seuItg.dft, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# # saveRDS(cluster_markers, file.path(PWD, "cluster_markers.rds"))
# # #cluster_markers <- readRDS(file.path(PWD, "cluster_markers.rds"))
# seuItg.dft <- ScaleData(seuItg.dft, features = rownames(seuItg.dft))
# top5 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# DoHeatmap(seuItg.dft, features = top5$gene)
# marker.plot <- VlnPlot(seuItg.dft, features = top2$gene[1:20], group.by = "seurat_clusters", pt.size = 0)
# marker.plot
# top5_markers <- cluster_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 5, order_by = avg_log2FC)
# write.csv(top5_markers, file = file.path(PWD, 'top5_markers.csv'), row.names=TRUE, quote = FALSE)

######################################################
# ## scGSEA Example
######################################################
#https://www.jianshu.com/p/00f069cb239c
# library(devtools)
# install_github('immunogenomics/presto')
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


##### merge Tm cluster #####
# unique(harmony_cluster$meta.cluster)
# merge.clsuer <- NULL
# for(i in harmony_cluster$meta.cluster){
#   if((i == "CD8.c01.Tn.MAL")){
#     merge.clsuer = append(merge.clsuer,"C1_Tn") 
#   }
#   if((i == "CD8.c02.Tm.IL7R") | (i == "CD8.c03.Tm.RPS12") |  (i ==  "CD8.c04.Tm.CD52") |  (i ==  "CD8.c05.Tem.CXCR5") |  (i ==  "CD8.c06.Tem.GZMK") | (i == "CD8.c10.Trm.ZNF683") | (i == "CD8.c17.Tm.NME1")){
#     merge.clsuer = append(merge.clsuer,"C2_Tm")
#   }
#   if((i == "CD8.c07.Temra.CX3CR1")){
#     merge.clsuer = append(merge.clsuer,"C3_Temra") 
#   }
#   if((i == "CD8.c08.Tk.TYROBP") | (i == "CD8.c09.Tk.KIR2DL4")){
#     merge.clsuer = append(merge.clsuer,"C4_Tk") 
#   }
#   if((i == "CD8.c11.Tex.PDCD1") | (i == "CD8.c12.Tex.CXCL13") | (i == "CD8.c13.Tex.myl12a") | (i == "CD8.c14.Tex.TCF7")){
#     merge.clsuer = append(merge.clsuer,"C5_Tex") 
#   }
#   if((i == "CD8.c15.ISG.IFIT1")){
#     merge.clsuer = append(merge.clsuer,"C6_ISG")
#   }
#   if((i == "CD8.c16.MAIT.SLC4A10")){
#     merge.clsuer = append(merge.clsuer,"C7_MAIT")
#   }}
# DimPlot(harmony_cluster, reduction = "umap", label=T, label.size = 3)
# pdf5 <- DimPlot(harmony_cluster, reduction = "umap", label=T, label.size = 3)
# ggsave(pdf5, dpi = 300, height = 4, width = 5.5, file = file.path(PWD, 'pdf5.pdf'))

# ## Get old/new cluster cell rate
# newcluster <- harmony_cluster@active.ident
# cluster0 <- names(newcluster[newcluster == "0.Tm.IL7R+ZNF683"])
# cluster1 <- names(newcluster[newcluster == "1.Tem.CXCR5+GZMK"])
# cluster2 <- names(newcluster[newcluster == "2.Tm.NME1+RPS12+CD52"])
# cluster3 <- names(newcluster[newcluster == "3.Tex.CXCL13+PDCD1+TCF7+myl12a"])
# cluster4 <- names(newcluster[newcluster == "4.Tn.MAL"])
# cluster5 <- names(newcluster[newcluster == "5.Temra.CX3CR1"])
# cluster6 <- names(newcluster[newcluster == "6.Temra.CX3CR1+Tk.TYROBP"])
# cluster7 <- names(newcluster[newcluster == "7.Tk.KIR2DL4+Tm.IL7R"])
# cluster8 <- names(newcluster[newcluster == "8.MAIT.SLC4A10"])
# cluster9 <- names(newcluster[newcluster == "9.NULL"])
# cluster10 <- names(newcluster[newcluster == "10.ISG.IFIT1"])
# oldcluster <- harmony_cluster$meta.cluster
# unique(oldcluster)
# oldcluster0 <- c(names(oldcluster[oldcluster == "CD8.c02.Tm.IL7R"]), names(oldcluster[oldcluster == "CD8.c10.Trm.ZNF683"]))
# oldcluster1 <- c(names(oldcluster[oldcluster == "CD8.c05.Tem.CXCR5"]), names(oldcluster[oldcluster == "CD8.c06.Tem.GZMK"]))
# oldcluster2 <- c(names(oldcluster[oldcluster == "CD8.c17.Tm.NME1" ]), names(oldcluster[oldcluster == "CD8.c03.Tm.RPS12"]), names(oldcluster[oldcluster == "CD8.c04.Tm.CD52"]))
# oldcluster3 <- c(names(oldcluster[oldcluster == "CD8.c12.Tex.CXCL13" ]), names(oldcluster[oldcluster == "CD8.c11.Tex.PDCD1"]), names(oldcluster[oldcluster == "CD8.c14.Tex.TCF7"]), names(oldcluster[oldcluster == "CD8.c13.Tex.myl12a"]))
# oldcluster7 <- c(names(oldcluster[oldcluster == "CD8.c09.Tk.KIR2DL4"]), names(oldcluster[oldcluster == "CD8.c02.Tm.IL7R"]))
# 
# intersect(cluster0, oldcluster0)
# intersect(cluster1, oldcluster1)
# intersect(cluster2, oldcluster2)
# intersect(cluster3, oldcluster3)
# intersect(cluster7, oldcluster7)
# length(intersect(cluster0, oldcluster0)) / length(cluster0)
# length(intersect(cluster1, oldcluster1)) / length(cluster1)
# length(intersect(cluster2, oldcluster2)) / length(cluster2)
# length(intersect(cluster3, oldcluster3)) / length(cluster3)
# length(intersect(cluster7, oldcluster7)) / length(cluster7)

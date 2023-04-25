#################################################################
# 2023-02-20
# Mapping and annotating query datasets
#################################################################

library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

## loda query data
dataset <-  readRDS("/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs/dataset3_GSE156728_PanCancer_CD8T.rds")
query.data <- subset(dataset, patient %in% c("PACA.P20190225","RC.P20190919","THCA.P20190621","UCEC.P20190717","OV.P20190304"))
table(query.data$patient)
query.list <- SplitObject(query.data, split.by = "patient")
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# for(i in seq(query.list)){
#   query.list[[i]]$percent.mt <- PercentageFeatureSet(object = query.list[[i]], pattern = "^MT-")
#   query.list[[i]]$percent.ribo <- PercentageFeatureSet(object = query.list[[i]], pattern = "^RP[LS]")}
  # query.list[[i]] <- CellCycleScoring(query.list[[i]], s.features = s.genes, g2m.features = g2m.genes, nbin=15)
  # query.list[[i]]$CC.Difference <- query.list[[i]]$S.Score - query.list[[i]]$G2M.Score}
# query.list <- lapply(query.list, FUN = function(sample){
#   sample_data <- subset(sample, subset = percent.mt < 7 & nCount_RNA < 40000 & nCount_RNA > 2000)
#   sample_data <- SCTransform(
#     sample,
#     assay = "RNA",
#     new.assay.name = "SCT",
#     variable.features.n = 3000,
#     method="glmGamPoi",
#     return.only.var.genes = FALSE,
#     vars.to.regress = c("percent.mt"),
#     # vars.to.regress = c("percent.mt", "CC.Difference"),
#     verbose = FALSE)
# })
for (i in 1:length(query.list)){
  query.list[[i]] <- subset(query.list[[i]], subset = percent.mt < 7 & nCount_RNA < 40000 & nCount_RNA > 2000)
  query.list[[i]] <- NormalizeData(query.list[[i]], verbose = FALSE)
  query.list[[i]] <- FindVariableFeatures(query.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)}
query.data <- merge(query.list[[1]],query.list[2:5])
df <- read.csv("/picb/lilab5/liuzipei/Rdata/darklist.csv", header = T)
darklist_df <- as.list(df)
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)
ribo <- darklist_df$Ribosome
MT <- darklist_df$Mitochondria
TCR <- grep(pattern = "^TR[ABDG][VDCJ]", x = rownames(query.data), value = TRUE)
#HLA <- grep(pattern = "^HLA", x = rownames(scRNAlist[[i]]), value = TRUE)
BCR <- grep(pattern = "^IGL[VC]", x = rownames(query.data), value = TRUE)
darklist <- c(TCR, DIG, ribo, BCR, MT, "MALAT1")
query.data <- query.data[-(which(rownames(query.data) %in% darklist)),] ## remove darklist genes

ref <- readRDS("/picb/lilab5/liuzipei/Rdata/Integration_Harmony/myobj_final2.rds") #参考
ref <- RunUMAP(ref, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony", reduction.key = "Uh_", return.model = T)
ref$celltype <- Idents(ref)
# anchors <- FindTransferAnchors(reference = ref,
#                                reference.reduction = "harmony",
#                                query = query.data,
#                                k.filter = NA )
# query.map <- MapQuery(anchorset = anchors,
#                    query = query.data, 
#                    reference = ref, 
#                    refdata = list(celltype = "celltype"),
#                    reduction.model = "umap.harmony")
# table(query.map$predicted.celltype == query.map$meta.cluster)
# 
# p1 <- DimPlot(query.map, reduction = "ref.umap", group.by = "predicted.celltype", label = T) + NoLegend()
# p2 <- DimPlot(query.map, reduction = "ref.umap", group.by = "celltype", label = T) + NoLegend()
# p1 | p2
ref[[]]
ref[["harmony2"]] <- CreateDimReducObject(ref[["harmony"]]@cell.embeddings,
                                          key = "harmony2_",
                                          loadings = ref[["pca"]]@feature.loadings,
                                          assay = "RNA")
anchors  <- FindTransferAnchors(reference = ref,
                               reference.reduction = "harmony2",
                               query = query.data ,
                               k.filter = NA)
query.map <- MapQuery(anchorset = anchors,
                      query = query.data, 
                      reference = ref, 
                      refdata = list(celltype = "celltype"),
                      reduction.model = "umap.harmony")

my.levels <- c("Prol.STMN1","HSP.HSPA6","MAIT.KLRB1","ISG.IFI6","Tn.CCR7",
               "Tk.TYROBP","Temra.FGFBP2","Tex.CXCL13","Tm","Tm.IL7R",
               "Trm.LDCRAD4","Tem.GZMK")
ref@meta.data$celltype <- factor(ref@meta.data$celltype, levels = my.levels)
query.map@meta.data$predicted.celltype <- factor(query.map@meta.data$predicted.celltype, levels = my.levels)
colors1 <- c("Tm"='#d62e2d', "Tm.IL7R"='#e9787a', "Trm.LDCRAD4"='#cd9c9b',
             "Tem.GZMK"='#684797', "Tem.GZMK"='#3377a9', " HSP.HSPA6"='#96c3d8',
             "Tn.CCR7"='#67a59b', "MAIT.KLRB1"='#70b866', "Tex.CXCL13"='#6a9a52',
             "ISG.IFI6"='#a5d38f', "Prol.STMN1"='#ab9697', "Tk.TYROBP"='#f19294',
             "Temra.FGFBP2"='#f5b375')

p1 <- DimPlot(query.map, reduction = "ref.umap", group.by = "predicted.celltype", label = F, cols = colors1) + NoLegend()
p2 <- DimPlot(ref, group.by = "celltype", label = T, cols =colors1) + NoLegend()            
p1|p2

p2           
query.map$meta.cluster
p4 <- DimPlot(query.map, reduction = "ref.umap", group.by = "meta.cluster", label = F) + NoLegend()
                                                                                                                                                                          

## mapping
pancreas.anchors <- FindTransferAnchors(reference = ref, query = query.data,
                                        dims = 1:30, reference.reduction = "pca" )
predictions <- TransferData(anchorset = pancreas.anchors, refdata = ref$celltype,
                            dims = 1:30)
query.data <- AddMetaData(query.data, metadata = predictions)
table(query.data$predicted.id == query.data$celltype)

Idents(seuItg.dft) <- seuItg.dft$meta.cluster
unique(seuItg.dft$meta.cluster)
new.cluster.ids <- c("Tem","Temra","Tm","MAIT","Tn",
                     "Tem","Tk","Tm","Tm","Trm",
                     "Tm","ISG","Tk","Tex","Tex",
                     "Tex","Tex")
names(new.cluster.ids) <- levels(seuItg.dft)
seuItg.dft <- RenameIdents(seuItg.dft, new.cluster.ids)

f <- subset(query.data, predicted.id != query.data$celltype)
table(f$celltype)
table(f$celltype) / table(query.data$celltype)
table(f$predicted.id)

t <- subset(query.data, predicted.id == query.data$celltype)
table(t$celltype)
table(query.data$celltype)
table(t$predicted.id) 



p1 <- VlnPlot(query.map, c("CXCL13",'ENTPD1','CCR7',"GZMK","ZNF683","IL7R"),
        group.by = "predicted.celltype",pt.size = 0, cols = colors1)
p2 <- VlnPlot(ref, c("CXCL13",'ENTPD1','CCR7',"GZMK","ZNF683","IL7R"),
        group.by = "celltype",pt.size = 0, cols = colors1)
p1
p2F

p3 <- VlnPlot(query.map, c("KLRG1","FCGR3A","FGFBP2","SLC4A10","IFIT1","STMN1"),
              group.by = "predicted.celltype",pt.size = 0, cols = colors1)
p4 <- VlnPlot(ref, c("KLRG1","FCGR3A","FGFBP2","SLC4A10","IFIT1","STMN1"),
              group.by = "celltype",pt.size = 0, cols = colors1)
p3
p4


levels(ref)
levels(query.map)
levels = my.levers
levels(query.map) <- ref.levers


VlnPlot(query.map, c('PDCD1',"CXCL13",'ENTPD1','TIGIT','LAG3','CTLA4','TOX','HAVCR2','LAYN'),
              group.by = "predicted.celltype",pt.size = 0, cols = colors1)
VlnPlot(query.data, c('CCR7', 'LEF1', 'SELL', "TCF7"),
        group.by = "predicted.id")
VlnPlot(query.data, c("KLRG1","CX3CR1","FCGR3A","FGFBP2"),
        group.by = "predicted.id")
VlnPlot(query.data, c("SLC4A10","KLPB1","ZBTB16","NCR3","RORC"),
        group.by = "predicted.id")
VlnPlot(query.data, c("IFIT1"),
        group.by = "predicted.id")
VlnPlot(query.data, c("GZMK"),
        group.by = "predicted.id")
VlnPlot(query.data, c("ZNF683","ITGAE"),
        group.by = "predicted.id")
VlnPlot(query.data, c("IL7R"),
        group.by = "predicted.id")

install.packages("symphony")


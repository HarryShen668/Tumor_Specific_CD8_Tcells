library(Seurat)
library(tidyverse)

rm(list = ls())
sce.ro <- readRDS("/home/shenhaoyu/xinda/dataset/Rosenberg_2022/5.clonetype/sce.lab_v0.Rds")
sce.ro <- sce.ro[,sce.ro$sampleID %in% c("4283","4317","4322","4323","4324","4325","4385","4394","4400","4421")]
unique(sce.ro$sampleID)
sce.lung <- readRDS("/home/shenhaoyu/xinda/dataset/NSCLC_Smith/1.expr/sce.lab_v0.Rds")
sce.lung$Histology <- "lung"
sce.wu <- readRDS("/picb/lilab5/liuzipei/Wu_2021/query.data.rds")
sce.wu$Histology <- "melanoma"
sceList <- list(sce.ro, sce.lung, sce.wu)

comb.sce.all <- merge(sceList[[1]], sceList[2:length(sceList)])
tumor_reactive_marker <- c("CXCL13","CTLA4","PDCD1","ENTPD1","LAG3","TIGIT","LAYN","HAVCR2")
comb.sce.all <- AddModuleScore(comb.sce.all, 
                               features = list(tumor_reactive_marker),
                               assay = "RNA",
                               name = "tumor_reactive_score")
comb.sce.all$antigen <- ifelse(is.na(comb.sce.all$antigen), "unknown", comb.sce.all$antigen)
comb.sce.all$TRT <-  ifelse(is.na(comb.sce.all$TRT ), "unknown", comb.sce.all$TRT)
comb.sce.all$Bystander <-  ifelse(is.na(comb.sce.all$Bystander ), "unknown", comb.sce.all$Bystander)
comb.sce.all$type <- ifelse(comb.sce.all$antigen == "tumor.specific" | comb.sce.all$TRT == "Tumor Reactive", "Tumor Reactive",
                            ifelse(comb.sce.all$antigen == "bystander" | comb.sce.all$Bystander == "Bystander", "Bystander", "unknown"))
unique(comb.sce.all$type)

ggplot(comb.sce.all@meta.data, mapping = aes(x = reorder(type,-tumor_reactive_score1), y = tumor_reactive_score1, fill = type)) +
  geom_boxplot() +
  ylim(-5,5)zz +
  labs(x = "", y = "TumorReactiveScore") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  scale_fill_aaas()

ref <- readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref_harmony_symphony.Rds")
ref1 <-  readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref1_harmony_symphony.Rds")
source("/home/liuzipei/projection.plot.R")
source("~/scripts/utils_seurat.R")
saveRDS(comb.sce.all, "~/project/scRNAseq_2022/comb.before")

comb.sce.all <- merge(comb_scc_bcc, comb_scc_bcc)
query.data <- comb.sce.all
query.map <- mapQuery( 
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)

query.map$clonotype_id <- ifelse(is.na(query.map$raw_clonotype_id), query.map$final.clonotype.family, query.map$raw_clonotype_id)
unique(query.map$clonotype_id)
clonotype_count <- tibble(clonotype_id = query.map$clonotype_id) %>% group_by(clonotype_id) %>% summarise(clonotype_count = n())
df <- left_join(query.map@meta.data, clonotype_count)
rownames(df) <- colnames(query.map)
unique(df$clonotype_count)
df$TR_predict <- ifelse(df$tumor_reactive_score1 > 0 & df$celltype_prob >= 0.8 & grepl("Tex", df$celltype) & df$clonotype_count >= 2, TRUE, FALSE)
df$virus_predict <- ifelse(df$type == "Bystander", T, F)
table(df$TR_predict)
table(df$virus_predict)


TR_clonotype_id <- df[df$TR_predict,]$clonotype_id
TR_clonotype_id
TR_memory_cells <- rownames(df[df$clonotype_id %in% TR_clonotype_id & df$clonotype_id != "0" & df$celltype_prob >= 1 & df$celltype %in% c("CD8.c03.Tm","CD8.c02.Tm.IL7R","CD8.c10.Trm.ZNF683","CD8.c05.Tem.CXCR5","CD8.c06.Tem.GZMK"),])
TR_memory_cells
TR_memory <- query.data[,TR_memory_cells]
virus_memory_cells <- rownames(df[df$virus_predict & df$celltype_prob >= 1 & df$celltype %in% c("CD8.c03.Tm","CD8.c02.Tm.IL7R","CD8.c10.Trm.ZNF683","CD8.c05.Tem.CXCR5","CD8.c06.Tem.GZMK"),])
virus_memory <- query.data[,virus_memory_cells]
table(TR_memory$Histology)
table(TR_memory$type)
table(virus_memory$Histology)
table(virus_memory$type)
TR_memory$predict <- "Tumor Reactive"
virus_memory$predict <- "Bystander"

DefaultAssay(TR_memory) <- "RNA"
DefaultAssay(virus_memory) <- "RNA"
TR_memory@assays$SCT <- NULL
virus_memory@assays$SCT <- NULL
comb.memory <- merge(TR_memory,virus_memory)
sceList <- SplitObject(comb.memory, split.by = "Histology")


marker.List <- lapply(sceList, function(cancer){
  DefaultAssay(cancer) <- "RNA"
  cancer <- NormalizeData(cancer, verbose = F)
  TCR <- grep(pattern = "^TR[ABGD][VCJD]", x = rownames(cancer@assays$RNA@data))
  BCR <- grep(pattern = "^IG[LKH][VCJD]", x = rownames(cancer@assays$RNA@data))
  darklist <- c(TCR,BCR)
  cancer@assays$RNA@data <- cancer@assays$RNA@data[-(darklist),]
  all.genes <- rownames(cancer)
  cancer <- ScaleData(cancer, feature = all.genes, verbose = F)
  Idents(cancer) <- cancer$predict
  markers <- FindAllMarkers(cancer, only.pos = TRUE,  min.pct = 0, logfc.threshold = 0)
})
de.List <- lapply(marker.List, function(marker){
  de.genes <- marker %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25 & cluster == "Tumor Reactive") %>%
    arrange(p_val) %>%
#    head(200) %>%
    select(gene)
  de.genes <- de.genes$gene
  TCR <- grep(pattern = "^TR[ABGD][VCJD]", x = de.genes, value = TRUE)
  BCR <- grep(pattern = "^IG[LKH][VCJD]", x = de.genes, value = TRUE)
  de.genes <- setdiff(de.genes, c(TCR,BCR)) %>%
    unique()
})
marker.List <- marker.List[2:3]
intersect(de.List[[1]], de.List[[2]])


projection.plot(ref, query = query.map, cols = color18, linesize = 0.5, pointsize = 0.5)
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(query.map , reduction = 'umap', group.by = 'celltype', shuffle = TRUE, cols = color18, label = T)
FeaturePlot(query.map, features = "tumor_reactive_score1", max.cutoff = 5, min.cutoff = -5)
query.map$tumor_reactive_score1
unique(query.map$celltype)




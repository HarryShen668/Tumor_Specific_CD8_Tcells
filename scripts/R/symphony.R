library(Seurat)
library(tidyverse)
library(harmony)

rm(list=ls())
source('~/scripts/utils_seurat.R')
color18 <- c('#d62e2d', '#e9787a', '#cd9c9b', '#684797', '#3377a9', '#96c3d8',
  '#67a59b', '#70b866', '#6a9a52', '#a5d38f', '#ab9697', '#f19294',
  '#f5b375', '#da8f6f', '#e0c880', '#f47d2f', '#ec8d63', '#e45d61',
  '#a65a35', '#4b9d47')
#source('/home/shenhaoyu/code/fig/color/single.dimplot.R')
ref <- readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/seuList_v1_harmony_map.Rds")
#ref$celltype <- Idents(ref)
ref$celltype <- ref$meta.cluster

ref <- RunHarmony.Seurat(ref, group.by.vars = "cancerType", plot_convergence = TRUE, verbose = F,project.dim = FALSE)

ref[['umap']] <- RunUMAP2(Embeddings(ref, 'harmony')[, 1:30],
                          assay='RNA', verbose=FALSE, umap.method='uwot', return.model=TRUE)

DimPlot(ref, reduction = "umap", group.by = 'meta.cluster',cols=color18)

#saveRDS(ref,'/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref_harmony_symphony.Rds')
#ref <- readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref_harmony_symphony.Rds")

ref1 <- buildReferenceFromSeurat(ref, assay = 'SCT',
                                 verbose=TRUE, save_umap=TRUE, save_uwot_path='~/cache_symphony_sct.uwot')

ref1$normalization_method = 'SCTransform'
saveRDS(ref1,'/home/liuzipei/project/ref1_harmony_symphony.Rds')
ref1 <-  readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref1_harmony_symphony.Rds")

### query.data
dataset <-  readRDS("/home/shenhaoyu/xinda/dataset/Rosenberg_2022/2.mapping/query.map.cd8.Rds")

#ref <- readRDS("/picb/lilab5/liuzipei/Rdata/Integration_Harmony/myobj_final_filter2.rds")
filtered_p <- names(table(dataset$patient)[table(dataset$patient) < 1000 & table(dataset$patient) > 500])
query.data <- subset(dataset, patient %in% filtered_p)
table(dataset$celltype)
DimPlot()

rm(dataset)


query.list <- SplitObject(query.data, split.by = "patient")
## load darklist
df <- read.csv("/picb/lilab5/liuzipei/Rdata/darklist.csv", header = T)
darklist_df <- as.list(df)
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)
ribo <- darklist_df$Ribosome
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)
# 创建包含所有Seurat对象中存在的基因和模块定义中包含的所有基因的列表
DIG <- intersect(rownames(query.data), DIG)

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

query.data <- merge(query.list.fil[[1]], query.list.fil[2:7])

query.data <- readRDS("/home/shenhaoyu/xinda/dataset/Rosenberg_2022/2.mapping/query.map.cd8.Rds")
query.data$
query.map <- mapQuery(
    query.data@assays$SCT@scale.data,
    query.data@meta.data,
    ref1,
    vars = 'sampleID',
    do_normalize = FALSE,
    return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
query.map@reductions$umap@cell.embeddings

options(repr.plot.height = 4, repr.plot.width = 10)
(DimPlot(ref, reduction = 'umap', group.by = 'meta.cluster', shuffle = TRUE, cols = color18) + labs(title = 'Original Reference (Clusters)')) +
  (DimPlot(query.map, reduction = 'umap', group.by = 'celltype', shuffle = TRUE, cols = color18) + labs(title = 'Mapped Query (Predicted Clusters)'))

source("/home/liuzipei/projection.plot.R")
projection.plot(ref, query = query.map, cols = color18, linesize = 0.5, pointsize = 0.5, title = "Rosenberg")

sce.ro <- readRDS("/home/shenhaoyu/xinda/dataset/Rosenberg_2022/5.clonetype/sce.lab_v0.Rds")
sce.lung <- readRDS("/home/shenhaoyu/xinda/dataset/NSCLC_Smith/1.expr/sce.lab_v0.Rds")

table(sce.lung$TRT)
table(sce.lung$Bystander)

filter.clones <- tibble(
  trb = sce.lung$raw_clonotype_id,
  type = sce.lung$TRT,
  Bystander = sce.lung$Bystander) %>%
  group_by(trb, type) %>%
  dplyr::summarise(count = n(), Bystander) %>%
  ungroup() %>%
  filter((type == "Tumor Reactive" & count >= 3) | (Bystander == "Bystander" & count >= 3))

filter.ro.clones <- tibble(
  trb = sce.ro$raw_clonotype_id,
  type = sce.ro$TRT,
  Bystander = sce.ro$Bystander) %>%
  group_by(trb, type) %>%
  dplyr::summarise(count = n(),Bystander) %>%
  ungroup() %>%
  filter((type == "Tumor Reactive" & count >= 3) | (Bystander == "Bystander" & count >= 3))

sce.ro <- sce.ro[,sce.ro$raw_clonotype_id %in% filter.ro.clones$trb]
sce.lung.pro <- sce.lung[,sce.lung$raw_clonotype_id %in% filter.clones$trb]
comb.sce <- merge(sce.ro, sce.lung.pro)

comb.sce.r <- subset(comb.sce, TRT == "Tumor Reactive")
table(comb.sce.r$TRT)
options(warn = -1)
comb.sce.r.map <- mapQuery(
  comb.sce.r@assays$SCT@scale.data,
  comb.sce.r@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
comb.sce.r.map  <- knnPredict.Seurat(comb.sce.r.map , ref1, label_transfer = 'celltype', confidence = TRUE)
comb.sce.r.filter.celltype <- tibble(
  barcode = comb.sce.r.map$new_barcode,
  celltype = comb.sce.r.map$celltype,
  prob = comb.sce.r.map$celltype_prob) %>%
  filter(prob >= 0.6) %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  filter( count >= 10)

comb.sce.r.clonotype <- tibble(
  barcode = comb.sce.r.map$new_barcode,
  clonotype = comb.sce.r.map$raw_clonotype_id,
  celltype = comb.sce.r.map$celltype,
  prob = comb.sce.r.map$celltype_prob) %>%
  filter(prob >= 0.6) %>%
  group_by(celltype, clonotype) %>%
  dplyr::summarise(count = n()) %>%
  arrange(desc(count))

source("~/projection.plot.R")
clonotype1 <- comb.sce.r.map[,comb.sce.r.map$raw_clonotype_id %in% c("MD043-011_clonotype13")]
clonotype2 <- comb.sce.r.map[,comb.sce.r.map$raw_clonotype_id %in% c("MD043-011_clonotype15")]
clonotype3<- comb.sce.r.map[,comb.sce.r.map$raw_clonotype_id %in% c("MD043-011_clonotype12")]

projection.plot(ref, query = clonotype1, cols = color18, linesize = 0.8, pointsize = 1.5, title = "clonotype1")
projection.plot(ref, query = clonotype2, cols = color18, linesize = 0.8, pointsize = 1.5, title = "clonotype2")
projection.plot(ref, query = clonotype3, cols = color18, linesize = 0.8, pointsize = 1.5, title = "clonotype3")


rep(x = "black", 17)
DimPlot(ref, reduction = "umap", label = FALSE, 
        group.by = "celltype", repel = TRUE, pt.size = NULL, 
        cols = rep(x = "black", 17)) + geom_point(df, mapping = aes(x = umap_1, y = umap_2, color = raw_clonotype), alpha = 0.6, 
                                      size = 1, shape = 17) + theme(aspect.ratio = 1) 
  # geom_density_2d(data = data.frame(query@reductions$umap@cell.embeddings),
  #                 mapping = aes(x = umap_1, y = umap_2), color = "black",
  #                 n = 200, h = 2, size = linesize) + ggtitle(title)


DimPlot(ref, reduction = "umap", label = FALSE, 
             group.by = labels.col, repel = TRUE, pt.size = ref.size, 
             cols = cols_use) +
  geom_point(df, mapping = aes(x = umap_1, y = umap_2, fill = raw_clonotype_id), alpha = 0.6, 
             size = pointsize, shape = 17, color = "red") 

df1 <- data.frame(clonotype@reductions$umap@cell.embeddings)
df1$barcode <- rownames(df1) 
  # geom_density_2d(data = data.frame(query@reductions$umap@cell.embeddings), 
  #                 mapping = aes(x = umap_1, y = umap_2), color = "black", 
  #                 n = 200, h = 2, size = linesize) + ggtitle(title) + 
  theme(aspect.ratio = 1)

comb.sce.r.map.filter <- comb.sce.r.map[,comb.sce.r.map$celltype %in% comb.sce.r.filter.celltype$celltype]
projection.plot(ref, query = comb.sce.r.map.filter, cols = color18, linesize = 0.8, pointsize = 1.5, title = "Tumor Reactive")

celltype_stat.ref <- ref@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  mutate(ref_frequence = count/length(colnames(ref))) %>%
  select(celltype,ref_frequence)
celltype_stat.r <- comb.sce.r.map@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  mutate(frequence = count/length(colnames(comb.sce.r.map))) %>%
  left_join(celltype_stat.ref) %>%
  mutate(relative_freq = frequence/ref_frequence)
sum1 <- sum(celltype_stat.r$relative_freq)
celltype_stat.r$relative_freq <- celltype_stat.r$relative_freq/sum1

ggplot(data = celltype_stat.r, mapping = aes(x = reorder(celltype,-relative_freq), y = relative_freq)) +
  geom_bar(stat = "identity", fill = "#9b2226") +
  geom_text(aes(label=count),vjust=-0.15,size=3) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "Frequence") +
  ggtitle("Tumor Reactive") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

comb.sce.b <- subset(comb.sce, TRT == "non Tumor Reactive")
comb.sce.b.map <- mapQuery(
  comb.sce.b@assays$SCT@scale.data,
  comb.sce.b@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
comb.sce.b.map  <- knnPredict.Seurat(comb.sce.b.map , ref1, label_transfer = 'celltype', confidence = TRUE)
comb.sce.b.filter.celltype <- tibble(
  barcode = comb.sce.b.map$new_barcode,
  celltype = comb.sce.b.map$celltype,
  prob = comb.sce.b.map$celltype_prob) %>%
  filter(prob >= 0.6) %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  filter( count >= 10)
comb.sce.b.map.filter <- comb.sce.b.map[ ,comb.sce.b.map$celltype %in% comb.sce.b.filter.celltype$celltype]
projection.plot(ref, query = comb.sce.b.map.filter, cols = color18, linesize = 0.8, pointsize = 1.5, title = "Bystander")
table(comb.sce.b.map.filter$celltype)

celltype_stat.b <- comb.sce.b.map@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  mutate(frequence = count/length(colnames(comb.sce.b.map))) %>%
  left_join(celltype_stat.ref) %>%
  mutate(relative_freq = frequence/ref_frequence)
sum2 <- sum(celltype_stat.b$relative_freq)
celltype_stat.b$relative_freq <- celltype_stat.b$relative_freq/sum2

ggplot(data = celltype_stat.b, mapping = aes(x = reorder(celltype,-relative_freq), y = relative_freq)) +
  geom_bar(stat = "identity", fill = "#457b9d") +
  geom_text(aes(label=count),vjust=-0.15,size=3) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "Frequence") +
  ggtitle("Bystander") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

tumor_specific <- sce.ro[,sce.ro$TRT == "Tumor Reactive"]
table(tumor_specific$TRT)

tumor_specific.map <- mapQuery(
  tumor_specific@assays$SCT@scale.data,
  tumor_specific@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
tumor_specific.map  <- knnPredict.Seurat(tumor_specific.map , ref1, label_transfer = 'celltype', confidence = TRUE)

source("/home/liuzipei/projection.plot.R")
tumor_specific.filter.celltype <- tibble(
  barcode = tumor_specific.map$new_barcode,
  celltype = tumor_specific.map$celltype,
  prob = tumor_specific.map$celltype_prob) %>%
  filter(prob >= 0.6) %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  filter( count >= 10)
tumor_specific.map.filter <- tumor_specific.map[ ,tumor_specific.map$celltype %in% tumor_specific.filter.celltype$celltype]
projection.plot(ref, query = tumor_specific.map.filter , cols = color18, linesize = 0.8, pointsize = 1.5, title = "treatment naive")
table(tumor_specific.map.filter$celltype)

celltype_stat.treatment.naive <- tumor_specific.map@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  mutate(frequence = count/length(colnames(tumor_specific.map))) %>%
  left_join(celltype_stat.ref) %>%
  mutate(relative_freq = frequence/ref_frequence)
sum3 <- sum(celltype_stat.treatment.naive$relative_freq)
celltype_stat.treatment.naive$relative_freq <- celltype_stat.treatment.naive$relative_freq/sum3

ggplot(data = celltype_stat.treatment.naive, mapping = aes(x = reorder(celltype,-relative_freq), y = relative_freq)) +
  geom_bar(stat = "identity", fill = "#9b2226") +
  geom_text(aes(label=count),vjust=-0.15,size=3) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "Frequence") +
  ggtitle("treatment naive") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

bystander <- sce.ro[,sce.ro$Bystander == "Bystander"]
table(bystander$Bystander)

bystander.map <- mapQuery(
  bystander@assays$SCT@scale.data,
  bystander@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
bystander.map  <- knnPredict.Seurat(bystander.map , ref1, label_transfer = 'celltype', confidence = TRUE)


options(repr.plot.height = 4, repr.plot.width = 10)
(DimPlot(ref, reduction = 'umap', group.by = 'meta.cluster', shuffle = TRUE, cols = color18, label = T) + labs(title = 'Original Reference (Clusters)')) 
(DimPlot(bystander.map , reduction = 'umap', group.by = 'celltype', shuffle = TRUE, cols = color18) + labs(title = 'Mapped Query (Predicted Clusters)'))

source("/home/liuzipei/projection.plot.R")
projection.plot(ref, query = bystander.map , cols = color18, linesize = 0.5, pointsize = 0.8)

comb.sce <- merge(sce.lung.pro, sce.ro)
comb.sce$TRT
comb.sce <- subset(comb.sce, TRT == "Tumor Reactive")
DefaultAssay(comb.sce) <- "RNA"
cxcl13_neg <- subset(comb.sce, CXCL13 == 0)

table(cxcl13_neg$Bystander)

cxcl13_neg.map <- mapQuery(
  cxcl13_neg@assays$SCT@scale.data,
  cxcl13_neg@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
cxcl13_neg.map  <- knnPredict.Seurat(cxcl13_neg.map , ref1, label_transfer = 'celltype', confidence = TRUE)


options(repr.plot.height = 4, repr.plot.width = 10)
(DimPlot(ref, reduction = 'umap', group.by = 'meta.cluster', shuffle = TRUE, cols = color18, label = T) + labs(title = 'Original Reference (Clusters)')) 
(DimPlot(cxcl13_neg.map , reduction = 'umap', group.by = 'celltype', shuffle = TRUE, cols = color18) + labs(title = 'Mapped Query (Predicted Clusters)'))

source("/home/liuzipei/projection.plot.R")
projection.plot(ref, query = cxcl13_neg.map , cols = color18, linesize = 0.5, pointsize = 0.8, title = "CXCL13-")

table(sce.ro$TRT)
sce.ro.r <- subset(sce.ro, TRT == "Tumor Reactive")
sce.ro.r.map <- mapQuery(
  sce.ro.r@assays$SCT@scale.data,
  sce.ro.r@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
sce.ro.r.map  <- knnPredict.Seurat(sce.ro.r.map , ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = sce.ro.r.map , cols = color18, linesize = 0.8, pointsize = 1.5, title = "Tumor Reactive")

sce.ro.b <- subset(sce.ro, TRT == "non Tumor Reactive")
sce.ro.b.map <- mapQuery(
  sce.ro.b@assays$SCT@scale.data,
  sce.ro.b@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
sce.ro.b.map  <- knnPredict.Seurat(sce.ro.b.map , ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = sce.ro.b.map , cols = color18, linesize = 0.8, pointsize = 1.5, title = "Bystander")
DimPlot(ref, reduction = "umap", group.by = 'meta.cluster',cols=color18, label = T)

table(sce.lung.pro$TRT)
sce.lung.pro.r <- subset(sce.lung.pro, TRT == "Tumor Reactive")
sce.lung.pro.r.map <- mapQuery(
  sce.lung.pro.r@assays$SCT@scale.data,
  sce.lung.pro.r@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
sce.lung.pro.r.map  <- knnPredict.Seurat(sce.lung.pro.r.map , ref1, label_transfer = 'celltype', confidence = TRUE)
sce.lung.pro.r.filter.celltype <- tibble(
  barcode = sce.lung.pro.r.map$new_barcode,
  celltype = sce.lung.pro.r.map$celltype,
  prob = sce.lung.pro.r.map$celltype_prob) %>%
  filter(prob >= 0.6) %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  filter( count >= 20)
sce.lung.pro.r.map.filter <- sce.lung.pro.r.map[ ,sce.lung.pro.r.map$celltype %in% sce.lung.pro.r.filter.celltype$celltype]
projection.plot(ref, query = sce.lung.pro.r.map.filter , cols = color18, linesize = 0.8, pointsize = 1.5, title = "post-ICB")
table(sce.lung.pro.r.map.filter$celltype)

celltype_stat.post <- sce.lung.pro.r.map@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  mutate(frequence = count/length(colnames(sce.lung.pro.r.map))) %>%
  left_join(celltype_stat.ref) %>%
  mutate(relative_freq = frequence/ref_frequence)
sum4 <- sum(celltype_stat.post$relative_freq)
celltype_stat.post$relative_freq <- celltype_stat.post$relative_freq/sum4

ggplot(data = celltype_stat.post, mapping = aes(x = reorder(celltype,-relative_freq), y = relative_freq)) +
  geom_bar(stat = "identity", fill = "#457b9d") +
  geom_text(aes(label=count),vjust=-0.15,size=3) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "Frequence") +
  ggtitle("Post-ICB") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
# query.map@meta.data <- query.map@meta.data %>%
#   mutate(meta.cluster=ifelse(meta.cluster %in% c('CD8.c03.Tm.RPS12','CD8.c04.Tm.CD52','CD8.c17.Tm.NME1'),'CD8.c03.Tm',meta.cluster))
# query.map@meta.data<-query.map@meta.data %>%
#   mutate(meta.cluster=ifelse(meta.cluster %in% c('CD8.c08.Tk.TYROBP','CD8.c09.Tk.KIR2DL4'),'CD8.c08.Tk',meta.cluster))
# query.map@meta.data<-query.map@meta.data %>%
#   mutate(meta.cluster=ifelse(celltype %in% c('CD8.c18.cycling'),'CD8.c18.cycling',meta.cluster))
# 
# query.map<-subset(query.map,meta.cluster %in% c('CD8.c13.Tex.myl12a'),invert=TRUE)
# 
# query.map<-subset(query.map,celltype_prob>=.4)
# 
# 
# table(query.map$celltype == query.map$meta.cluster)
# 
# 
# f <- subset(query.map, celltype != query.map$meta.cluster)
# table(f$meta.cluster) / table(query.map$meta.cluster)[-14]
# DefaultAssay(query.map) <- "SCT"
# VlnPlot(query.map, c("TCF7", "CXCL13", "IL7R", "ENTPD1", "GZMK", "ZNF683"), group.by = "meta.cluster")
# VlnPlot(query.map, c("TCF7", "CXCL13", "IL7R", "ENTPD1", "GZMK", "ZNF683"), group.by = "celltype")
# 
# 
# # 计算混淆矩阵
# cm <- confusionMatrix(factor(query.map@meta.data$meta.cluster,levels = c("CD8.c01.Tn.MAL", "CD8.c02.Tm.IL7R",
#                                                                          "CD8.c03.Tm","CD8.c05.Tem.CXCR5","CD8.c06.Tem.GZMK","CD8.c07.Temra.CX3CR1",
#                                                                          "CD8.c08.Tk","CD8.c10.Trm.ZNF683","CD8.c11.Tex.PDCD1",
#                                                                          "CD8.c12.Tex.CXCL13","CD8.c14.Tex.TCF7"
#                                                                          ,"CD8.c15.ISG.IFIT1","CD8.c16.MAIT.SLC4A10","CD8.c18.cycling"),ordered = TRUE),
#                       factor(query.map@meta.data$celltype,levels = c("CD8.c01.Tn.MAL", "CD8.c02.Tm.IL7R",
#                                                                      "CD8.c03.Tm","CD8.c05.Tem.CXCR5","CD8.c06.Tem.GZMK","CD8.c07.Temra.CX3CR1",
#                                                                      "CD8.c08.Tk","CD8.c10.Trm.ZNF683","CD8.c11.Tex.PDCD1",
#                                                                      "CD8.c12.Tex.CXCL13","CD8.c14.Tex.TCF7"
#                                                                      ,"CD8.c15.ISG.IFIT1","CD8.c16.MAIT.SLC4A10","CD8.c18.cycling"),ordered = TRUE))
# 
# 
# # 将混淆矩阵转换为比例矩阵
# cm_prop <- prop.table(cm$table, margin = 1) # margin = 1 表示按行计算比例
# 
# 
# # 载入 pheatmap 包
# library(pheatmap)
# # 绘制混淆矩阵热图
# pheatmap(cm_prop,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          fontsize = 10,
#          color = colorRampPalette(c("white", "blue"))(100))


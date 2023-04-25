library(Seurat)
library(tidyverse)

rm(list=ls())
#source('/home/shenhaoyu/code/fig/color/single.dimplot.R')
source('/home/dongdanyue/script/system_config/color.R')
color18<-c(color15[c(1,3,5,7,9,11,13)],'#17BECFFF')

######################################################
# ## Read data
######################################################
dirPjtHome <- '/home/,/xinda/'
dir_raw<-'/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/GEX/'
dirOut <- '/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/'
tcr_folder <- '/picb/lilab5/liuzipei/Rosenberg_2022/vdj_out/outs/'
dirList = list.files("/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/GEX")
seuList <- list()
for(i in 1:length(dirList)){
  sa <- dirList[i]
  lab_ID<-sa
  if(grepl("KSA", sa)){
    type <- 'KSA'
  }
  else if(grepl("QVD", sa)){
    type<-'QVD'
  }
  else if(grepl("PD1pos", sa)){
    type <-"PD1pos"
  }
  mtx<-Read10X(file.path(dir_raw, sa))
  if(class(mtx) == "list"){
    mtx = mtx[["Gene Expression"]]
  }
  seuList[[lab_ID]]<-CreateSeuratObject(mtx,project = "7.HNSCC_Ahmed",min.cells = 10,min.features = 300)
  seuList[[lab_ID]] <-  RenameCells(
    object = seuList[[lab_ID]],
    new.names = paste0(sa,'_', colnames(x = seuList[[lab_ID]]))
    )
  seuList[[lab_ID]]@meta.data$sampleID <- lab_ID
  seuList[[lab_ID]]@meta.data$type <- type
}
merge.data <- merge(seuList[[1]], y = seuList[2:length(seuList)])
seuList <- SplitObject(merge.data, split.by = "sampleID")

######################################################
# ## data QC
######################################################
dirOutQC <- file.path(dirOut, 'QC')
## load darklist
df <- read.csv("/picb/lilab5/liuzipei/Rdata/darklist.csv", header = T)
darklist_df <- as.list(df)
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)
ribo <- darklist_df$Ribosome
DIG <- c(darklist_df$Heat.shock.protein, darklist_df$Dissociation)
DIG <- intersect(rownames(merge.data), DIG)
rm(merge.data)

if(! file.exists(dirOutQC)){ dir.create(dirOutQC) }
for(sa in names(seuList)){
  seuList[[sa]]$pctMT <- PercentageFeatureSet(
    object = seuList[[sa]], pattern = "^MT-"
  )
}
for(sa in names(seuList)){
  seuList[[sa]]$pctRB <- PercentageFeatureSet(
    object = seuList[[sa]], pattern = "^RP[SL]"
  )
}
for(sa in names(seuList)){
  seuList[[sa]]$pctTCR <- PercentageFeatureSet(
    object = seuList[[sa]], pattern = "^TR[ABGD][VCJD]"
  )
}
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
for(sa in names(seuList)){
  seuList[[sa]] <- CellCycleScoring(seuList[[sa]],s.features = s.genes, g2m.features = g2m.genes)
  seuList[[sa]]$CC.Difference <- seuList[[sa]]$S.Score - seuList[[sa]]$G2M.Score
}

for(sa in names(seuList)){
  seuList[[sa]]$pctFeature <- seuList[[sa]]$nCount_RNA/seuList[[sa]]$nFeature_RNA
}

options(warn = -1)
for(sa in names(seuList)){
  seuList[[sa]] <- AddModuleScore(
    object = seuList[[sa]],
    features = list(DIG),
    name = 'DIG_SCORE')
}
metaData.big <- rbind(seuList[[1]]@meta.data,
                      seuList[[2]]@meta.data,
                      seuList[[3]]@meta.data,
                      seuList[[4]]@meta.data,
                      seuList[[5]]@meta.data,
                      seuList[[6]]@meta.data,
                      seuList[[7]]@meta.data,
                      seuList[[8]]@meta.data,
                      seuList[[9]]@meta.data,
                      seuList[[10]]@meta.data,
                      seuList[[11]]@meta.data,
                      seuList[[12]]@meta.data,
                      seuList[[13]]@meta.data,
                      seuList[[14]]@meta.data,
                      seuList[[15]]@meta.data,
                      seuList[[16]]@meta.data)


show.nCount <- ggplot(metaData.big, aes(x = sampleID, y = nCount_RNA, fill = sampleID)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('nCount_RNA') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+ylim(0,5000)

show.nFeature <- ggplot(metaData.big, aes(x = sampleID, y = nFeature_RNA, fill = sampleID)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('nFeature_RNA') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+ylim(0,5000)

show.pctMT <- ggplot(metaData.big, aes(x = sampleID, y = pctMT, fill = sampleID)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('pctMT') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+ylim(0,25)

show.pctRB <- ggplot(metaData.big, aes(x = sampleID, y = pctRB, fill = sampleID)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('pctRB') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+ylim(0,10)

show.pctTCR <- ggplot(metaData.big, aes(x = sampleID, y = pctTCR, fill = sampleID)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('pctTCR') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

show.pctFeature <- ggplot(metaData.big, aes(x = sampleID, y = pctFeature, fill = sampleID)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('pctFeature') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+ylim(0,6)

show.DIG_SCORE <- ggplot(metaData.big, aes(x = sampleID, y = DIG_SCORE1, fill = sampleID)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle('DIG_SCORE') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))


basicQC <- show.nCount + show.nFeature + show.pctMT + show.pctRB + show.pctTCR

p1 <- ggplot(metaData.big, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point() +
  theme_bw()

p2 <- ggplot(metaData.big, aes(x = nCount_RNA, y = pctMT)) +
  geom_point() +
  theme_bw()

p3 <- ggplot(metaData.big, aes(x = nFeature_RNA, y = pctMT)) +
  geom_point() +
  theme_bw()

p4 <- ggplot(metaData.big, aes(x = nCount_RNA, y = pctRB)) +
  geom_point() +
  theme_bw()

p5 <- ggplot(metaData.big, aes(x = nFeature_RNA, y = pctRB)) +
  geom_point() +
  theme_bw()

p6 <- ggplot(metaData.big, aes(x = nCount_RNA, y = pctTCR)) +
  geom_point() +
  theme_bw()

p7 <- ggplot(metaData.big, aes(x = nFeature_RNA, y = pctTCR)) +
  geom_point() +
  theme_bw()

show.int.nCount.nFeature <- p1 + facet_wrap( ~ sampleID, nrow = 2)
show.int.nCount.petMT <- p2 + facet_wrap( ~ sampleID, nrow = 2)
show.int.nFeature.perMT <- p3 + facet_wrap( ~ sampleID, nrow = 2)
show.int.nCount.perRB <- p4 + facet_wrap( ~ sampleID, nrow = 2)
show.int.nFeature.perRB <- p5 + facet_wrap( ~ sampleID, nrow = 2)
show.int.nCount.perTCR <- p6 + facet_wrap( ~ sampleID, nrow = 2)
show.int.nFeature.perTCR <- p7 + facet_wrap( ~ sampleID, nrow = 2)

advancedQC <- show.int.nCount.nFeature / show.int.nCount.petMT / show.int.nFeature.perMT /show.int.nCount.perRB/show.int.nFeature.perRB/show.int.nCount.perTCR/show.int.nFeature.perTCR
basicQC
ggsave(basicQC, file = file.path(dirOutQC, 'basicQC.pdf'),
       dpi = 300, width = 13.5, height = 5)
ggsave(advancedQC, file = file.path(dirOutQC, 'advancedQC.pdf'),
       dpi = 300, width = 10, height = 17.5)

histogram<-ggplot(metaData.big,aes(x=nCount_RNA,y=..count..))+geom_histogram(binwidth=500,colour="#9987ce")+scale_x_continuous(breaks=c(15000,30000,45000))+geom_vline(xintercept =25000,color="red",linetype="dashed")+geom_vline(xintercept =800,color="red",linetype="dashed")+theme_classic()+xlim(c(0,50000))+facet_wrap(~sampleID,scales="free")
ggsave(histogram, file = file.path(dirOutQC, 'histogram.pdf'),
       dpi = 300, width = 15, height = 15)

seuList = seuList[c(-4,-5,-6)]
seuList.fil <- list()
for(i in 1:length(seuList)){
  na <- names(seuList)[i]
  
  seuList.fil[[na]] <- subset(seuList[[na]],
                              subset = pctMT < 15 &pctRB < 50 & pctRB > 1&
                                pctTCR < 2& pctTCR > 0 & nCount_RNA<25000 &
                                nCount_RNA>1000& nFeature_RNA>400& nFeature_RNA<3500&
                                pctFeature<5&DIG_SCORE1<5)
  
  seuList.fil[[na]] <- SCTransform(
    seuList.fil[[na]],
    assay = "RNA",
    new.assay.name = "SCT",
    variable.features.n = 3000,
    return.only.var.genes = FALSE,
    vars.to.regress = c("nCount_RNA","CC.Difference",'pctMT','DIG_SCORE1'), method = "glmGamPoi")
}

######################################################
# ## Filtering
######################################################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
proliferation<- c('ZWINT', 'E2F1', 'FEN1', 'FOXM1', 'H2AFZ', 'HMGB2', 'MCM2', 'MCM3', 'MCM4', 'MCM5',
                  'MCM6', 'MKI67', 'MYBL2', 'PCNA', 'PLK1', 'CCND1', 'AURKA', 'BUB1', 'TOP2A', 'TYMS',
                  'DEK', 'CCNB1','CCNE1')

## filter darklist genes
scRNAlist_filter <- NULL
darklist <- NULL
for(i in seq(seuList.fil)){
  na <- names(seuList.fil)[i]
  TCR <- grep(pattern = "^TR[ABGD][VCJD]", x = rownames(seuList.fil[[i]]), value = TRUE)
  #  HLA <- grep(pattern = "^HLA", x = rownames(scRNAlist[[i]]), value = TRUE)
  BCR <- grep(pattern = "^IG[LKH][VCJD]", x = rownames(seuList.fil[[i]]), value = TRUE)
  darklist[[i]] <- c(TCR, DIG, ribo, BCR, "MALAT1",s.genes,g2m.genes,proliferation)
  scRNAlist_filter[[na]] <- seuList.fil[[na]]
  }
#  scRNAlist_filter[[na]] <- subset(seuList.fil[[na]], (CD3G>1|CD3D>1|CD3E>1)&(CD8A>2|CD8B>2),slot = "counts")}

darklist<-Reduce(union, darklist)
saveRDS(scRNAlist_filter,  file = file.path('/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/', 'seuList_v0_filter.Rds'))
#scRNAlist_filter <- readRDS(file.path('/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/', 'seuList_v0_filter.Rds'))
scRNAlist_filter <- readRDS(file = file.path('/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/', 'seuList_v0_filter.Rds'))

######################################################
# ## TCR clonetype
######################################################
tcr <- read.csv("/picb/lilab5/liuzipei/scRNAseq_2022/TCR/7.HNSCC_Ahmed_all_filtered_contig_annotations.csv", sep=",")%>%
  select(new_barcode,raw_clonotype_id) %>%
  unique()
for(i in 1:length(scRNAlist_filter)){
  na <- names(scRNAlist_filter)[i]
  meta<-scRNAlist_filter[[na]]@meta.data
  meta<-meta %>% mutate(new_barcode=rownames(meta))
  meta<-left_join(meta,tcr) %>% mutate(raw_clonotype_id=paste0(na,'_',raw_clonotype_id))
  rownames(meta)<-meta$new_barcode
  scRNAlist_filter[[na]]@meta.data<-meta
  scRNAlist_filter[[na]]<-subset(scRNAlist_filter[[na]],subset = raw_clonotype_id %in% c(paste0(na,'_','None'),paste0(na,'_','NA')),invert=TRUE)
}
meta <- scRNAlist_filter[[1]]@meta.data
meta<-meta %>% mutate(new_barcode=rownames(meta))

meta<-left_join(meta,tcr) %>% mutate(raw_clonotype_id=paste0(na,'_',raw_clonotype_id))
saveRDS(scRNAlist_filter,  file = file.path('/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/', 'seuList_v1_filter.Rds'))

scRNAlist_filter[[1]]

tcr.family <- read.csv("~/CD8_TCR_clones_tumor_specificity.csv", header = T)
tcr.family$Tumor.specific
scRNAlist_filter[[1]]$final.clonotype.family
scRNAlist_filter1 <- lapply(scRNAlist_filter, function(seuobj){
  if(seuobj$final.clonotype.family %in% tcr.family$Tumor.specific){
    seuobj$antigen <- tumor.specific
  }
  else if(seuobj$final.clonotype.family %in% tcr.family$Viral){
    seuobj$antigen <- bystander
  }
  else{
    seuobj$antigen <- NA
  }
})
scRNAlist_filter1 <- lapply(scRNAlist_filter, function(seuobj){
  seuobj@meta.data <- seuobj@meta.data %>%
    mutate(antigen=ifelse(final.clonotype.family %in% tcr.family$Tumor.specific, "tumor.specific",
                          ifelse(final.clonotype.family %in% tcr.family$Viral, "bystander", NA)))
  return(seuobj)
})
table(scRNAlist_filter1[[1]]$antigen)



######################################################
# ## MAPPING
######################################################
ref <- readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/seuList_v1_harmony_map.Rds")
ref <- RunHarmony.Seurat(ref, group.by.vars = "cancerType", plot_convergence = TRUE, verbose = F,project.dim = FALSE)
ref[['umap']] <- RunUMAP2(Embeddings(ref, 'harmony')[, 1:30],
                          assay='RNA', verbose=FALSE, umap.method='uwot', return.model=TRUE)
ref$celltype <- ref$meta.cluster
ref1 <- buildReferenceFromSeurat(ref, assay = 'SCT',
                                 verbose=TRUE, save_umap=TRUE, save_uwot_path='~/cache_symphony_sct.uwot')

ref1$normalization_method = 'SCTransform'
saveRDS(ref1,'~/project/ref1_harmony_symphony.Rds')
source('~/scripts/utils_seurat.R')
color18 <- c('#d62e2d', '#e9787a', '#cd9c9b', '#684797', '#3377a9', '#96c3d8',
             '#67a59b', '#70b866', '#6a9a52', '#a5d38f', '#ab9697', '#f19294',
             '#f5b375', '#da8f6f', '#e0c880', '#f47d2f', '#ec8d63', '#e45d61',
             '#a65a35', '#4b9d47')
ref1 <- readRDS("~/project/ref1_harmony_symphony.Rds")


query.data <- merge(scRNAlist_filter1[[1]], scRNAlist_filter1[2:10])
scRNAlist_filter1[[10]]$antigen
query.data$antigen
t <- query.data[,query.data$antigen == "tumor.specific"]
t.map <- mapQuery(
  t@assays$SCT@scale.data,
  t@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
t.map <- knnPredict.Seurat(t.map, ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = t.map, cols = color18, linesize = 0.5, pointsize = 1, title = "Tumor.specific")

b <- query.data[,query.data$antigen == "bystander"]
b.map <- mapQuery(
  b@assays$SCT@scale.data,
  b@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
b.map <- knnPredict.Seurat(b.map, ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = b.map, cols = color18, linesize = 0.5, pointsize = 1, title = "bystander")

saveRDS(query.data, "/picb/lilab5/liuzipei/Wu_2021/query.data.rds")
options(warn = -1)
query.map <- mapQuery(
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)

options(repr.plot.height = 4, repr.plot.width = 10)
(DimPlot(ref, reduction = 'umap', group.by = 'meta.cluster', shuffle = TRUE, cols = color18) + labs(title = 'Original Reference (Clusters)')) +
  (DimPlot(query.map, reduction = 'umap', group.by = 'celltype', shuffle = TRUE, cols = color18) + labs(title = 'Mapped Query (Predicted Clusters)'))

source("/home/liuzipei/projection.plot.R")
projection.plot(ref, query = query.map, cols = color18, linesize = 0.5, pointsize = 0.5, title = "Wu")

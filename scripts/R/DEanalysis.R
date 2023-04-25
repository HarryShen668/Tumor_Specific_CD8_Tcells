library(Seurat)
library(tidyverse)

rm(list = ls())
sce.ro <- readRDS("/home/shenhaoyu/xinda/dataset/Rosenberg_2022/5.clonetype/sce.lab_v0.Rds")
sce.lung <- readRDS("/home/shenhaoyu/xinda/dataset/NSCLC_Smith/1.expr/sce.lab_v0.Rds")
sce.wu <- readRDS("/picb/lilab5/liuzipei/Wu_2021/query.data.rds")

filter.clones <- tibble(
  trb = sce.lung$raw_clonotype_id,
  type = sce.lung$TRT,
  Bystander = sce.lung$Bystander) %>%
  group_by(trb, type) %>%
  dplyr::summarise(count = n(), Bystander = first(Bystander)) %>%
  ungroup() %>%
  filter((type == "Tumor Reactive" & count >= 3) | (Bystander == "Bystander" ))

filter.ro.clones <- tibble(
  trb = sce.ro$raw_clonotype_id,
  type = sce.ro$TRT,
  Bystander = sce.ro$Bystander) %>%
  group_by(trb, type) %>%
  dplyr::summarise(count = n(), Bystander = first(Bystander)) %>%
  ungroup() %>%
  filter((type == "Tumor Reactive" & count >= 3) | (Bystander == "Bystander" ))

filter.wu.clones <- tibble(
  tcr = sce.wu$final.clonotype.family,
  type = sce.wu$antigen) %>%
  group_by(tcr,type) %>%
  dplyr::summarise(count = n()) %>%
  ungroup() %>%
  filter((type == "tumor.specific" & count >= 3) | (type == "bystander" ))

sce.ro <- sce.ro[,sce.ro$raw_clonotype_id %in% filter.ro.clones$trb]
sce.lung.pro <- sce.lung[,sce.lung$raw_clonotype_id %in% filter.clones$trb]
sce.wu <- sce.wu[,sce.wu$final.clonotype.family %in% filter.wu.clones$tcr]

sce.lung.pro$Histology <- "lung cancer"
sce.wu$Histology <- "melanoma"
sce.wu$antigen

comb.sce <- merge(sce.lung.pro, sce.ro)
TCR <- grep(pattern = "^TR[ABGD][VCJD]", x = rownames(comb.sce), value = TRUE)
BCR <- grep(pattern = "^IG[LKH][VCJD]", x = rownames(comb.sce), value = TRUE)
comb.sce@meta.data <- comb.sce@meta.data %>% mutate(antigen=ifelse(TRT=='Tumor Reactive','tumor.specific',"bystander"))
comb.sce <- merge(comb.sce, sce.wu)
table(comb.sce$antigen)
unique(comb.sce$Histology)
table(comb.sce$Histology)
sceList <- SplitObject(comb.sce, split.by = "Histology")
sceList.temp <- list() 
sceList.temp[["lung cancer"]] <- sceList[["lung cancer"]]
sceList.temp[["melanoma"]] <- sceList[["melanoma"]]
sceList.temp[["colorectal cancer"]] <- merge(sceList[["Rectal"]],sceList[["Colon"]])
sceList <- sceList.temp

marker.List <- lapply(sceList, function(cancer){
  DefaultAssay(cancer) <- "RNA"
  cancer <- NormalizeData(cancer, verbose = F)
  TCR <- grep(pattern = "^TR[ABGD][VCJD]", x = rownames(cancer@assays$RNA@data))
  BCR <- grep(pattern = "^IG[LKH][VCJD]", x = rownames(cancer@assays$RNA@data))
  darklist <- c(TCR,BCR)
  cancer@assays$RNA@data <- cancer@assays$RNA@data[-(darklist),]
  all.genes <- rownames(cancer)
  cancer <- ScaleData(cancer, feature = all.genes, verbose = F)
  Idents(cancer) <- cancer$antigen
  markers <- FindAllMarkers(cancer, only.pos = TRUE,  min.pct = 0.1, logfc.threshold = 0.2)
})

marker.List[[1]]
marker.List[[2]]
marker.List[[3]]

de.List <- lapply(marker.List, function(marker){
  de.genes <- marker %>%
    filter(p_val < 0.05 & abs(avg_log2FC) > 0.25 & cluster == "tumor.specific") %>%
    arrange(p_val) %>%
    head(200) %>%
    select(gene)
  de.genes <- de.genes$gene
  TCR <- grep(pattern = "^TR[ABGD][VCJD]", x = de.genes, value = TRUE)
  BCR <- grep(pattern = "^IG[LKH][VCJD]", x = de.genes, value = TRUE)
  de.genes <- setdiff(de.genes, c(TCR,BCR)) %>%
    unique()
})

length(rownames(sceList[[1]]))

options(warn = -1)
library(VennDiagram)
library(RColorBrewer)
library(ggvenn)
ggvenn(data = list(
  `lung cancer` = de.List$`lung cancer`,
  melanoma = de.List$melanoma,
  `colorectal cancer` = de.List$`colorectal cancer`),
  show_percentage = FALSE,
  text_size = 6,
  fill_color =  c("#6930c3","#5e60ce","#720026"))


three_intersect <- intersect(de.List[[3]],intersect(de.List[[1]], de.List[[2]]))
three_intersect
two_intersect <- unique(c(intersect(de.List[[1]], de.List[[2]]),intersect(de.List[[1]], de.List[[3]]),intersect(de.List[[3]], de.List[[2]])))
two_intersect
saveRDS(three_intersect, "/picb/lilab5/liuzipei/Rdata/three_intersect.rds")
saveRDS(two_intersect, "/picb/lilab5/liuzipei/Rdata/two_intersect.rds")

clusters <- readRDS("/home/shenhaoyu/xinda/dataset/TCR.clone.classification/degene.cluster.rds")          
clusters[[3]]

counts <- c()
count_list <- lapply(clusters, function(cluster){
  count <- length(intersect(cluster, two_intersect))
  counts <- c(counts,count)
  counts
})

length(two_intersect) -40-13-28
# 绘制条形图
library(ggplot2)
library(ggsci)
count_list
# 创建示例数据集
df <- data.frame(
  cluster = c(paste0("C",c(2,4,6)),"NA"),
#  count = c(13,16,6,2,0,5,6),
  count = c(13,40,28,17),
  x = "DEgenes"
)

# ggplot(data = df, aes(x = "", y = count, fill = cluster)) + 
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar(theta = "y") +
#   geom_text(aes(label = paste0(cluster, ": ", count)), position = position_stack(vjust = 0.5)) +
#   theme_void() +
#   scale_fill_npg()

# 计算每个元素所占的百分比
df$percent <- df$count / sum(df$count) * 100

# 绘制堆叠条形图
ggplot(df, aes(x = x, y = percent, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.05) +
  coord_flip() +
  theme_void() +
  theme(legend.position = c(0.9,0.7)) +
  scale_fill_manual(values = c("#9b2226","#bb3e03","#005f73","#b7b7a4"))

  

## 超几何假设检验
#phyper(194, 1850+195, 15220-1850-195, 195+596, lower.tail = F)
library(gmp)
enrich_pvalue <- function(N, A, B, k)
{
  m <- A + k
  n <- B + k
  i <- k:min(m,n)
  
  as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}

lung_melanoma <- enrich_pvalue(20149, 140+17, 119+38, 22+21)
lung_colorectal <- enrich_pvalue(20149, 140+22, 124+38,21+17)
colorectal_melanoma <- enrich_pvalue(20149, 124+17, 119+22, 38+21)

lung_melanoma
lung_colorectal
colorectal_melanoma

Tem <- readRDS("/home/shenhaoyu/xinda/dataset/TCR.clone.classification/celltpye_marker/degene.list.ALL.TEM.rds")
Tex <- readRDS("/home/shenhaoyu/xinda/dataset/TCR.clone.classification/celltpye_marker/degene.list.ALL.Tex.rds")
Tm <- readRDS("/home/shenhaoyu/xinda/dataset/TCR.clone.classification/celltpye_marker/degene.list.ALL.TM.rds")
Trm <- readRDS("/home/shenhaoyu/xinda/dataset/TCR.clone.classification/celltpye_marker/degene.list.ALL.TrM.rds")


c2 <- clusters[[2]]
c4 <- clusters[[1]]
c6 <- clusters[[3]]
de.tumor <- c(c2,c4,c6)
de.tumor

class1 <- de.tumor[!(de.tumor %in% Tem) & !(de.tumor %in% Tex) & !(de.tumor %in% Tm) & !(de.tumor %in% Trm)]
class2 <- c(de.tumor[de.tumor %in% Tem & !(de.tumor %in% Tex) & !(de.tumor %in% Tm) & !(de.tumor %in% Trm)],
            de.tumor[!de.tumor %in% Tem & (de.tumor %in% Tex) & !(de.tumor %in% Tm) & !(de.tumor %in% Trm)],
            de.tumor[!de.tumor %in% Tem & !(de.tumor %in% Tex) & (de.tumor %in% Tm) & !(de.tumor %in% Trm)],
            de.tumor[!de.tumor %in% Tem & !(de.tumor %in% Tex) & !(de.tumor %in% Tm) & (de.tumor %in% Trm)])
class2 <- c(de.tumor[de.tumor %in% Tem & !(de.tumor %in% Tex) & !(de.tumor %in% Tm) & !(de.tumor %in% Trm)],
            de.tumor[!de.tumor %in% Tem & (de.tumor %in% Tex) & !(de.tumor %in% Tm) & !(de.tumor %in% Trm)],
            de.tumor[!de.tumor %in% Tem & !(de.tumor %in% Tex) & (de.tumor %in% Tm) & !(de.tumor %in% Trm)],
            de.tumor[!de.tumor %in% Tem & !(de.tumor %in% Tex) & !(de.tumor %in% Tm) & (de.tumor %in% Trm)])
class3 <- de.tumor[!(de.tumor %in% c(class1, class2))]
class3
class <- list(class1,class2,class3)
saveRDS(class, "/picb/lilab5/liuzipei/scRNAseq_2022/class.rds")

class4 <- c(de.tumor[de.tumor %in% Tem &  (de.tumor %in% Tex) & (de.tumor %in% Tm) & (de.tumor %in% Trm)])

df <- data.frame(
  cluster = c(paste0("Class",c(1,2,3))),
  #  count = c(13,16,6,2,0,5,6),
  count = c(length(class1),length(class2),length(class3)),
  x = "DEgenes"
)

# ggplot(data = df, aes(x = "", y = count, fill = cluster)) + 
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar(theta = "y") +
#   geom_text(aes(label = paste0(cluster, ": ", count)), position = position_stack(vjust = 0.5)) +
#   theme_void() +
#   scale_fill_npg()

# 计算每个元素所占的百分比
df$percent <- df$count / sum(df$count) * 100

# 绘制堆叠条形图
ggplot(df, aes(x = x, y = percent, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.05) +
  coord_flip() +
  theme_void() +
  theme(legend.position = c(0.9,0.7)) +
  scale_fill_manual(values = c("#b7b7a4","#9b2226","#005f73"))




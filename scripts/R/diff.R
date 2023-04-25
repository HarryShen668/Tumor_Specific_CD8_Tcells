library(Seurat)
library(tidyverse)

rm(list = ls())
sce.ro <- readRDS("/home/shenhaoyu/xinda/dataset/Rosenberg_2022/5.clonetype/sce.lab_v0.Rds")
sce.lung <- readRDS("/home/shenhaoyu/xinda/dataset/NSCLC_Smith/1.expr/sce.lab_v0.Rds")


table(sce.lung$TRT)
table(sce.lung$Bystander)

filter.clones <- tibble(
  trb = sce.lung$raw_clonotype_id,
  type = sce.lung$TRT,
  Bystander = sce.lung$Bystander) %>%
  group_by(trb, type) %>%
  summarise(count = n(), Bystander) %>%
  ungroup() %>%
  filter((type == "Tumor Reactive" & count >= 3) | (Bystander == "Bystander" ))

filter.ro.clones <- tibble(
  trb = sce.ro$raw_clonotype_id,
  type = sce.ro$TRT,
  Bystander = sce.ro$Bystander) %>%
  group_by(trb, type) %>%
  summarise(count = n(),Bystander) %>%
  ungroup() %>%
  filter((type == "Tumor Reactive" & count >= 3) | (Bystander == "Bystander" ))

sce.ro <- sce.ro[,sce.ro$raw_clonotype_id %in% filter.ro.clones$trb]
sce.lung.pro <- sce.lung[,sce.lung$raw_clonotype_id %in% filter.clones$trb]



comb.sce <- merge(sce.lung.pro, sce.ro)
comb.sce <- ScaleData(comb.sce,assay='SCT',split.by = 'sampleID',scale.max = 3,scale.min = -3)
comb.sce <- ScaleData(comb.sce,assay='RNA',split.by = 'sampleID',scale.max = 3,scale.min = -3)
comb.sce <- ScaleData(comb.sce,assay='RNA',scale.max = 3)
comb.sce <- ScaleData(comb.sce,assay='SCT',scale.max = 3)


comb.sce@meta.data <- comb.sce@meta.data %>%
  mutate(type=ifelse(TRT=='Tumor Reactive','Tumor Reactive',"Bystander"))

tmp.rank <- tibble(
  cell = colnames(comb.sce),
  type = comb.sce$type) %>%
  arrange(type)

use.genes <- c("CXCL13","LAYN","ENTPD1","ITGAE","TOX","PKM","GAPDH","HIF1A",
               "GZMB","GZMA","GZMH","CD74","PFN1","EVL","PSME1","EZH2","FASLG","IFNG","IL7R")
#use.genes <- degene
tmp.matr <- comb.sce@assays$SCT@scale.data[use.genes,tmp.rank$cell]
colnames(tmp.matr) <- NULL

library(ComplexHeatmap)
library(circlize)
ha = HeatmapAnnotation(
  type = tmp.rank$type,
  col = list(
    cluster = c("Tumor Reactive" = "red","Bystander" = "#009ACD")),
  simple_anno_size = unit(0.4, "cm")
)
col_fun <- colorRamp2(seq(from = -3, to = 3, length.out = 10), rev(RColorBrewer::brewer.pal(10, "RdBu")))
ComplexHeatmap::Heatmap(tmp.matr, cluster_columns = F, cluster_rows = F,
                        col = col_fun,
                        top_annotation = ha)

#----- diff expression
comb.sce <- merge(sce.lung.pro, sce.ro)
comb.sce@meta.data <- comb.sce@meta.data %>% mutate(type=ifelse(TRT=='Tumor Reactive','Tumor Reactive',"Bystander"))
DefaultAssay(comb.sce) <- "RNA"
comb.sce <- NormalizeData(comb.sce)
comb.sce <- ScaleData(comb.sce, do.center = F, do.scale = F,assay ="RNA" )
table(comb.sce$type)
Idents(comb.sce) <- comb.sce$type

library(tidyverse)
library(ggrepel)
markers <- FindAllMarkers(comb.sce, only.pos = T, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0)
marker.gene <- c("CXCL13","ITGAE","CTLA4","GAPDH","ENTPD1","BATF","GZMB","LAYN","BAG3","TNFRSF9","PFN1","CXCR3","MT2A","CALM1","CCR7","HLA-DRB1","CD74","RPS4X","ATP5ME","RPS14","IL21R","PTPRC")
markers %>%
  dplyr::rename(FC = avg_log2FC) %>%
  dplyr::mutate(name = ifelse(gene %in% c(marker.gene), gene, NA)) %>%
  dplyr::mutate(FC = ifelse(cluster != "Tumor Reactive", -FC,FC)) %>%
  dplyr::mutate(sig = ifelse(p_val_adj < 0.01 & abs(FC) > 0.25, "yes", "no")) %>%
  ggplot(aes(FC, -log10(p_val_adj))) +
  geom_point(aes(color = sig), size = 0.9) +
  geom_text_repel(aes(label = name),
                  parse = T,
                  size = 3.5,
                  box.padding = 0.5,
                  segment.color = "black",
                  show.legend = FALSE) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#D6D6D6","#CD2626"))





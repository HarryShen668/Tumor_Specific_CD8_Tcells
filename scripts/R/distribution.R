library(Seurat)
library(tidyverse)
library(ggsci)
source("/home/liuzipei/projection.plot.R")
source("~/scripts/utils_seurat.R")
ref <- readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref_harmony_symphony.Rds")
ref1 <-  readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref1_harmony_symphony.Rds")


## 数据预处理，只保留tumor reactive and bystander
rm(list = ls())
sce.ro <- readRDS("/home/shenhaoyu/xinda/dataset/Rosenberg_2022/5.clonetype/sce.lab_v0.Rds")
sce.ro <- sce.ro[,sce.ro$sampleID %in% c("4283","4317","4322","4323","4324","4325","4385","4394","4400","4421")] ## 去除PD1-POST
unique(sce.ro$sampleID)
sce.lung <- readRDS("/home/shenhaoyu/xinda/dataset/NSCLC_Smith/1.expr/sce.lab_v0.Rds")
sce.wu <- readRDS("/picb/lilab5/liuzipei/Wu_2021/query.data.rds")

filter.clones <- tibble(
  trb = sce.lung$raw_clonotype_id,
  type = sce.lung$TRT,
  Bystander = sce.lung$Bystander) %>%
  group_by(trb, type) %>%
  dplyr::summarise(count = n(), Bystander = first(Bystander)) %>%
  ungroup() %>%
  filter((type == "Tumor Reactive" & count >= 3) | (Bystander == "Bystander"))

filter.ro.clones <- tibble(
  trb = sce.ro$raw_clonotype_id,
  type = sce.ro$TRT,
  Bystander = sce.ro$Bystander) %>%
  group_by(trb, type) %>%
  dplyr::summarise(count = n(), Bystander = first(Bystander)) %>%
  ungroup() %>%
  filter((type == "Tumor Reactive" & count >= 3) | (Bystander == "Bystander"))

filter.wu.clones <- tibble(
  tcr = sce.wu$final.clonotype.family,
  type = sce.wu$antigen) %>%
  group_by(tcr,type) %>%
  dplyr::summarise(count = n()) %>%
  ungroup() %>%
  filter((type == "tumor.specific" & count >= 3) | (type == "bystander" ))
  
sce.ro <- sce.ro[,sce.ro$raw_clonotype_id %in% filter.ro.clones$trb]
sce.ro.r <- sce.ro[,sce.ro$TRT == "Tumor Reactive"]

sce.lung <- sce.lung[,sce.lung$raw_clonotype_id %in% filter.clones$trb]
sce.lung.r <- sce.lung[,sce.lung$TRT == "Tumor Reactive"]
unique(sce.lung.r$sampleID)
sce.lung.r.mpr <- sce.lung.r[,sce.lung.r$sampleID %in% "MD01-005"]
sce.lung.r.nmpr <- sce.lung.r[,sce.lung.r$sampleID %in% c("MD043-011","NY016-014")]

sce.wu <- sce.wu[,sce.wu$final.clonotype.family %in% filter.wu.clones$tcr]
sce.wu.r <- sce.wu[,sce.wu$antigen == "tumor.specific"]
comb.sce.r <- merge(sce.ro.r, sce.lung.r)
comb.sce.r <- merge(comb.sce.r, sce.wu.r)
unique(comb.sce.r$sampleID)
saveRDS(comb.sce.r, "/picb/lilab5/liuzipei/scRNAseq_2022/comb.sce.r.rds")
comb.sce.all <- merge(sce.ro, sce.lung)
comb.sce.all$antigen <- ifelse(comb.sce.all$TRT == "Tumor Reactive", "tumor.specific", "bystander")
comb.sce.all <- merge(comb.sce.all, sce.wu)
unique(comb.sce.all$antigen)
unique(comb.sce.all$clonotype_id)
sce.lung$raw_clonotype_id
sce.wu$final.clonotype.family
comb.sce.all$clonotype = ifelse(!is.na(comb.sce.all$clonotype_id), comb.sce.all$clonotype_id, 
                                ifelse(!is.na(comb.sce.all$raw_clonotype_id), comb.sce.all$raw_clonotype_id,
                                       ifelse(!is.na(comb.sce.all$final.clonotype.family), comb.sce.all$final.clonotype.family, "None")))
unique(comb.sce.all$clonotype)
clonotype_count <- tibble(clonotype = comb.sce.all$clonotype) %>% group_by(clonotype) %>% summarise(clonotype_count = n())
comb.sce.all@meta.data <- left_join(comb.sce.all@meta.data, clonotype_count)
rownames(comb.sce.all@meta.data) <- colnames(comb.sce.all)
unique(comb.sce.all$clonotype_count)
#SplitObject(comb.sce.all,split.by = "clonotype_count")
saveRDS(comb.sce.all, "/picb/lilab5/liuzipei/scRNAseq_2022/comb.sce.all.rds")


## mapping 
# ref <- readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref_harmony_symphony.Rds")
# ref1 <-  readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/ref1_harmony_symphony.Rds")
# source("/home/liuzipei/projection.plot.R")
# source("~/scripts/utils_seurat.R")
# query.data <- sce.lung.r
# query.map <- mapQuery(
#   query.data@assays$SCT@scale.data,
#   query.data@meta.data,
#   ref1,
#   vars = 'sampleID',
#   do_normalize = FALSE,
#   return_type = 'Seurat' )# return a Seurat object
# query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
# query.map1 <- query.map[,query.map$sampleID %in% "MD01-005"]
# query.map2 <- query.map[,query.map$sampleID %in% c("MD01-005","MD043-011","NY016-014")]
# query.map$sampleID
# projection.plot(ref, query = query.map1, cols = color18, linesize = 0.5, pointsize = 0.5, title = "MPR")
# projection.plot(ref, query = query.map2, cols = color18, linesize = 0.5, pointsize = 0.5, title = "NON-MPR")

comb.data <- readRDS("/home/shenhaoyu/xinda/dataset/TCR.clone.classification/comb.data.all.rds")
unique(comb.data$orig.ident)
comb.data <- comb.data[,comb.data$orig.ident != "7.HNSCC_Ahmed"]
unique(comb.data$sampleID)


## proportion of T and B per cluster
query.data <- comb.data
unique(query.data$sampleID)
query.map <- mapQuery(
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)

type_counts <- query.map@meta.data %>% 
  group_by(type) %>%
  summarize(total = sum(n()))

df <- query.map@meta.data %>% 
  group_by(celltype, type) %>%
  summarize(freq = n()) %>%
  mutate(freq_norm = ifelse(type == "Tumor Reactive", freq / type_counts$total[2] , freq / type_counts$total[1]))
df_wider <- tidyr::pivot_wider(df, id_cols = celltype, names_from = type, values_from = freq_norm)  
df_wider$proportion <- df_wider$`Tumor Reactive`/df_wider$Bystander
df_wider
ordered_celltype <- df_wider %>% 
  arrange(desc(proportion)) %>%
  select(celltype) 

df$celltype <- factor(df$celltype, levels = as.vector(ordered_celltype$celltype))
query.map@meta.data
ggplot(data = df) +
  geom_bar(aes(x = celltype, y = freq_norm, fill = type),
           position = "fill", color = "black", stat = "identity") +
  scale_fill_npg() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "") +
  ggtitle("Proportion of Tumor Reactive and Bystander per Cluster") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

## ICB-post and treatment naive 
comb.data$treatment <- ifelse(grepl("MD", comb.data$sampleID) | grepl("NY", comb.data$sampleID)|comb.data$sampleID %in% c("p2","p11"), "post", "pre")
query.data <- comb.data
unique(query.data$sampleID)
query.map <- mapQuery(
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)

type_counts <- query.map@meta.data %>% 
  group_by(treatment) %>%
  summarize(total = sum(n()))
type_counts
df <- query.map@meta.data %>% 
  group_by(celltype, treatment) %>%
  summarize(freq = n()) %>%
  mutate(freq_norm = ifelse(treatment == "pre", freq / type_counts$total[2] , freq / type_counts$total[1]))
df_wider <- tidyr::pivot_wider(df, id_cols = celltype, names_from = treatment, values_from = freq_norm)  
df_wider$proportion <- df_wider$pre/df_wider$post
df_wider
ordered_celltype <- df_wider %>% 
  arrange(desc(proportion)) %>%
  select(celltype) 

df$celltype <- factor(df$celltype, levels = as.vector(ordered_celltype$celltype))
ggplot(data = df) +
  geom_bar(aes(x = celltype, y = freq_norm, fill = treatment),
           position = "fill", color = "black", stat = "identity") +
  scale_fill_npg() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "") +
  ggtitle("Proportion of ICB-post and Treatment-Naive per Cluster") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



## T and B
comb.r <- comb.data[,comb.data$type == "Tumor Reactive"]
comb.b <- comb.data[,comb.data$type == "Bystander"]
unique(comb.data$sampleID)
unique(comb.data$type)
query.data <- comb.r
unique(query.data$sampleID)
query.map <- mapQuery(
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = query.map, cols = color18, linesize = 0.8, pointsize = 0.5, title = "Tumor Reactive")

query.data <- comb.b
query.map <- mapQuery(
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = query.map, cols = color18, linesize = 0.8, pointsize = 0.8, title = "Bystander")

## pre and post 
comb.data$treatment <- ifelse(grepl("MD", comb.data$sampleID) | grepl("NY", comb.data$sampleID)|comb.data$sampleID %in% c("p2","p11"), "post", "pre")
unique(comb.data$treatment)
table(comb.data$treatment)
comb.post <- comb.data[,comb.data$treatment == "post"]
unique(comb.post$sampleID)
comb.pre <- comb.data[,comb.data$treatment == "pre"]
unique(comb.pre$sampleID)
query.data <- comb.pre
query.map <- mapQuery(
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = query.map, cols = color18, linesize = 0.8, pointsize = 0.8, title = "Treatment Naive")

query.data <- comb.post
query.map <- mapQuery(
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
projection.plot(ref, query = query.map, cols = color18, linesize = 0.8, pointsize = 0.8, title = "ICB-Post")


## barplot 
seuList <- list()
seuList[["TR"]] <- comb.r
seuList[["bystander"]] <- comb.b
seuList[["pre"]] <- comb.pre
seuList[["post"]] <- comb.post

celltype_stat.ref <- ref@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(count = n()) %>%
  mutate(ref_frequence = count/length(colnames(ref))) %>%
  select(celltype,ref_frequence)

celltype_stat <- list()
celltype_stat <- lapply(seuList, function(seuobj){
  query.data <- seuobj
  query.map <- mapQuery(
    query.data@assays$SCT@scale.data,
    query.data@meta.data,
    ref1,
    vars = 'sampleID',
    do_normalize = FALSE,
    return_type = 'Seurat' )# return a Seurat object
  query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
  celltype_stat <- query.map@meta.data %>%
    group_by(celltype) %>%
    dplyr::summarise(count = n()) %>%
    mutate(frequence = count/length(colnames(query.map))) %>%
    left_join(celltype_stat.ref) %>%
    mutate(relative_freq = frequence/ref_frequence)
  sum2 <- sum(celltype_stat$relative_freq)
  celltype_stat$relative_freq <- celltype_stat$relative_freq/sum2
  return(celltype_stat)
})


ggplot(data = celltype_stat[[1]], mapping = aes(x = reorder(celltype,-relative_freq), y = relative_freq)) +
  geom_bar(stat = "identity", fill = "#e76f51") +
  geom_text(aes(label=count),vjust=-0.15,size=3) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "Frequence") +
  ggtitle("Tumor Specific") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggplot(data = celltype_stat[[2]], mapping = aes(x = reorder(celltype,-relative_freq), y = relative_freq)) +
  geom_bar(stat = "identity", fill = "#264653") +
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

ggplot(data = celltype_stat[[3]], mapping = aes(x = reorder(celltype,-relative_freq), y = relative_freq)) +
  geom_bar(stat = "identity", fill = "#e76f51") +
  geom_text(aes(label=count),vjust=-0.15,size=3) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "Frequence") +
  ggtitle("Treatment Naive") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggplot(data = celltype_stat[[4]], mapping = aes(x = reorder(celltype,-relative_freq), y = relative_freq)) +
  geom_bar(stat = "identity", fill = "#264653") +
  geom_text(aes(label=count),vjust=-0.15,size=3) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "Frequence") +
  ggtitle("ICB-Post") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
unique(comb.sce.all$antigen)
comb.sce.all.r <- comb.sce.all[,comb.sce.all$antigen == "tumor.specific"]
comb.sce.all.b <- comb.sce.all[,comb.sce.all$antigen == "bystander"]
top20.r <- tibble(clonotype = comb.sce.all.r$clonotype,
       clonotype_count = comb.sce.all.r$clonotype_count) %>%
  group_by(clonotype) %>%
  arrange(desc(clonotype_count)) %>%
  distinct(clonotype, .keep_all = T) %>%
  select(clonotype) %>%
  head(20)
top20.b <- tibble(clonotype = comb.sce.all.b$clonotype,
                  clonotype_count = comb.sce.all.b$clonotype_count) %>%
  group_by(clonotype) %>%
  arrange(desc(clonotype_count)) %>%
  distinct(clonotype, .keep_all = T) %>%
  select(clonotype) %>%
  head(20)

tibble(clonotype = comb.sce.all.r$clonotype,
       clonotype_count = comb.sce.all.r$clonotype_count) %>%
  group_by(clonotype) %>%
  arrange(desc(clonotype_count)) %>%
  distinct(clonotype, .keep_all = T)

tibble(clonotype = comb.sce.all.b$clonotype,
       clonotype_count = comb.sce.all.b$clonotype_count) %>%
  group_by(clonotype) %>%
  arrange(desc(clonotype_count)) %>%
  distinct(clonotype, .keep_all = T)

query.data <- comb.sce.all
query.map <- mapQuery(
  query.data@assays$SCT@scale.data,
  query.data@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
unique(query.map$celltype)
query.map <- query.map[,query.map$clonotype %in% c(top20.b$clonotype,top20.r$clonotype) & !(query.map$celltype %in% c("CD8.c08.Tk","CD8.c01.Tn.MAL","CD8.c15.ISG.IFIT1","CD8.c16.MAIT.SLC4A10"))]
query.map$cluster <- ifelse(query.map$celltype %in% c("CD8.c10.Trm.ZNF683"), "Trm",
                            ifelse(query.map$celltype %in% c("CD8.c06.Tem.GZMK","CD8.c05.Tem.CXCR5"), "Tem", 
                                   ifelse(query.map$celltype %in% c("CD8.c03.Tm","CD8.c02.Tm.IL7R","CD8.c07.Temra.CX3CR1"), "Tm",
                                          ifelse(grepl("Tex",query.map$celltype), "Tex", 
                                                 ifelse(query.map$celltype %in% c("CD8.c18.cycling"), "Cycling", "other")))))
unique(query.map$cluster)
unique(query.map$antigen)
factor(query.map$cluster, levels = c("Trm","Tem","Tm","Tex","Cycling"))


query.map$cluster <- factor(query.map$cluster, levels = c("Cycling","Tex","Tm","Tem","Trm"))

query.map.r <- query.map[,query.map$antigen == "tumor.specific"]
query.map.b <- query.map[,query.map$antigen == "bystander"]


ordered_clonotypes <- query.map.r@meta.data %>%
  group_by(clonotype) %>%
  mutate(trm_prop = sum(cluster == "Trm") / n()) %>%
  arrange(desc(trm_prop)) %>%
  pull(clonotype) %>%
  unique()
query.map.r@meta.data$clonotype <- factor(query.map.r@meta.data$clonotype,
                                          levels = ordered_clonotypes)
ggplot(data = query.map.r@meta.data) +
  geom_bar(mapping = aes(x = clonotype, fill = cluster),
           position = "fill") +
  scale_fill_npg() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "Proportion") +
  ggtitle("Tumor Specific") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

ordered_clonotypes <- query.map.b@meta.data %>%
  group_by(clonotype) %>%
  mutate(trm_prop = sum(cluster == "Trm") / n()) %>%
  arrange(desc(trm_prop)) %>%
  pull(clonotype) %>%
  unique()
query.map.b@meta.data$clonotype <- factor(query.map.b@meta.data$clonotype,
                                          levels = ordered_clonotypes)
ggplot(data = query.map.b@meta.data) +
  geom_bar(mapping = aes(x = clonotype, fill = cluster),
           position = "fill") +
  scale_fill_npg() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "", y = "Proportion") +
  ggtitle("Bystander") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# ggplot(data = query.map@meta.data) +
#   geom_bar(mapping = aes(x = clonotype, fill = cluster),
#            position = "fill") +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
#   scale_fill_npg() +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black"),
#         axis.text.y = element_text(color = "black")) +
#   labs(x = "", y = "Proportion") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   facet_wrap(~antigen)
comb.sce.all <- readRDS("/picb/lilab5/liuzipei/scRNAseq_2022/comb.sce.all.rds")
query.map <- mapQuery(
  comb.sce.all@assays$SCT@scale.data,
  comb.sce.all@meta.data,
  ref1,
  vars = 'sampleID',
  do_normalize = FALSE,
  return_type = 'Seurat' )# return a Seurat object
query.map <- knnPredict.Seurat(query.map, ref1, label_transfer = 'celltype', confidence = TRUE)
query.map <- query.map[,!(query.map$celltype %in% c("CD8.c08.Tk","CD8.c01.Tn.MAL","CD8.c15.ISG.IFIT1","CD8.c16.MAIT.SLC4A10"))]
query.map$cluster <- ifelse(query.map$celltype %in% c("CD8.c10.Trm.ZNF683"), "Trm",
                            ifelse(query.map$celltype %in% c("CD8.c06.Tem.GZMK","CD8.c05.Tem.CXCR5"), "Tem", 
                                   ifelse(query.map$celltype %in% c("CD8.c03.Tm","CD8.c02.Tm.IL7R","CD8.c07.Temra.CX3CR1"), "Tm",
                                          ifelse(grepl("Tex",query.map$celltype), "Tex", 
                                                 ifelse(query.map$celltype %in% c("CD8.c18.cycling"), "Cycling", "other")))))

tibble(clonotype = query.map$clonotype,
       cluster = query.map$cluster,
       type = query.map$antigen) %>%
  group_by(clonotype)


df1 <- tibble(clonotype = query.map$clonotype,
       cluster = query.map$cluster,
       type = query.map$antigen) %>%
  group_by(clonotype) %>%
  mutate(Trm = sum(cluster == "Trm") / n()) %>%
  mutate(Tem = sum(cluster == "Tem") / n()) %>%
  mutate(Tm = sum(cluster == "Tm") / n()) %>%
  mutate(Tex = sum(cluster == "Tex") / n()) %>%
  mutate(Cycling = sum(cluster == "Cycling") / n()) %>%
  unique()
df_long <- tidyr::pivot_longer(df1, cols = c("Trm","Tem","Tm","Tex","Cycling"), names_to = "cell_type", values_to = "proportion")
df_long
saveRDS(df_long, "/picb/lilab5/liuzipei/scRNAseq_2022/data/df_long20220410.rds")


ggplot(df_long, aes(x = type, y = proportion, color = type)) +
  geom_boxplot() +
  labs(x = "", y = "Proportion of cluster per clonotype", color = NULL) +
  facet_wrap(~ cell_type, nrow = 2, scales = "free_y") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(color = "black"))



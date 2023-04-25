#####################################################
# merge public scRNAseq Data
# Create on 2023-01-07
#####################################################

library(Seurat)
library(stringr)
rm(list=ls())

## Yang
setwd("/picb/lilab5/liuzipei/Yang_2022/GEX/")
folders <- list.files("/picb/lilab5/liuzipei/Yang_2022/GEX/")
folders
sceList <- lapply(folders, function(folder){
  CreateSeuratObject(counts = Read10X((file.path(folder,"outs/filtered_feature_bc_matrix"))),
                                      project = folder,min.cells = 3, min.features = 200)
})
folders
# rename <- c(rep("p3",3),rep("p2",5),rep("p1",5),rep("p4",2))
# rename
sce.big <- merge(sceList[[1]], y = sceList[2:length(sceList)],
                 add.cell.ids = folders, project = "Yang_2022")
table(sce.big$orig.ident)
a <- as.data.frame(sce.big@assays$RNA@counts)
saveRDS(sce.big,file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs/","Yang2022_CD8T.rds"))

setwd("/picb/lilab5/liuzipei/Yang_2022/TCR/")
samples <- list.files("/picb/lilab5/liuzipei/Yang_2022/TCR/")
samples
TcrList <- lapply(samples,function(sample){
  df <- read.csv(file.path(sample, "outs/filtered_contig_annotations.csv"), header = T)
  df$sample <- sample
  df$barcode <- paste(df$sample, df$barcode, sep = "_")
  df
})
tcr.big <- TcrList[[1]]
for(i in 2:length(TcrList)){
  tcr.big <- rbind(TcrList[[i]], tcr.big)
}
unique(tcr.big$high_confidence)
tcr.filtered <- tcr.big %>% 
  group_by(barcode,chain) %>%
  filter(umis == max(umis)) %>%
  ungroup()
write.csv(tcr.filtered, file.path("/picb/lilab5/liuzipei/Yang_2022/TCR/", "all_filtered_contig_annotations.csv"),
                quote = FALSE, row.names = FALSE)
tcr <- read.csv(file.path("/picb/lilab5/liuzipei/Yang_2022/TCR/", "all_filtered_contig_annotations.csv"))
tcr$barcode <- ifelse(grepl(pattern = "^4237", tcr$sample), paste0("p1_", tcr$barcode), tcr$barcode)
tcr$barcode <- ifelse(grepl(pattern = "^4234", tcr$sample), paste0("p2_", tcr$barcode), tcr$barcode)
tcr$barcode <- ifelse(grepl(pattern = "^4129", tcr$sample), paste0("p3_", tcr$barcode), tcr$barcode)
tcr$barcode <- ifelse(grepl(pattern = "^4369", tcr$sample), paste0("p4_", tcr$barcode), tcr$barcode)
write.csv(tcr, file = file.path("/picb/lilab5/liuzipei/Yang_2022/TCR/", "all_filtered_contig_annotations.csv"),
          quote = FALSE, row.names = FALSE)

TRA <- tcr[tcr$chain == "TRA",]
TRB <- tcr[tcr$chain == "TRB",]
# 新建空数据框
new_df <- data.frame(barcode = character(),
                     clonotype = character(),
                     stringsAsFactors = FALSE)
# 遍历 TRA 数据框
for (i in 1:nrow(TRA)) {
  # 查找 TRB 中与当前行 barcode 值相等的行
  idx <- which(TRB$barcode == TRA$barcode[i])
  if (length(idx) > 0) {
    # 如果有匹配行，则生成 clonotype，并添加到新数据框中
    clonotype <- paste(TRA$cdr3[i], TRB$cdr3[idx], sep = "_")
    new_df[nrow(new_df) + 1,] <- list(TRA$barcode[i], clonotype)
  }
}
write.csv(new_df, file = file.path("/picb/lilab5/liuzipei/Yang_2022/TCR/", "clonotype.csv"),
          quote = FALSE, row.names = FALSE)
tumor_reactive <- read.csv("../1.csv", header = T)
tumor_reactive_clonotype <- paste(tumor_reactive$CDR3, tumor_reactive$CDR3B, sep = "_")
tumor_reactive_clonotype
df <- new_df[new_df$clonotype %in% tumor_reactive_clonotype, ]
df
write.csv(df, file = file.path("/picb/lilab5/liuzipei/Yang_2022/TCR/", "tumor_reactive_clonotype.csv"),
          quote = FALSE, row.names = FALSE)

##1.tcr
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/1.HCC_zhang")
tcr <- read.csv("TCR_info.csv", header = T, skip = 1)
tcr_filter <- tcr[tcr$Productive..Alpha1. == "TRUE" & tcr$Productive.Beta1. == "TRUE",] 
tcr_filter <- tcr_filter[!(is.na(tcr_filter$Productive..Alpha1.) & is.na(tcr_filter$Productive.Beta1.)),] 
write.csv(tcr_filter, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/TCR","filtered_tcr.csv"),
          quote = FALSE, row.names = FALSE)
##2.tcr
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/2.NSCLC_zhang/")
tcr <- read.csv("TCR_info.csv", header = T, skip = 1)
tcr_filter <- tcr[tcr$Productive..Alpha1. == "TRUE" & tcr$Productive.Beta1. == "TRUE",] 
tcr_filter <- tcr_filter[!(is.na(tcr_filter$Productive..Alpha1.) & is.na(tcr_filter$Productive.Beta1.)),] 
tcr_filter
write.csv(tcr_filter, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/TCR","filtered_tcr.csv"),
          quote = FALSE, row.names = FALSE)

##6.NSCLC_Smith
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor/GEX/")
folders <- list.files("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor/GEX/")
folders
sceList <- lapply(folders, function(folder){
  CreateSeuratObject(counts = Read10X(folder), project = folder,
                     min.cells = 3, min.features = 200)
})
sce.big <- merge(sceList[[1]], y = sceList[2:length(sceList)],
                 add.cell.ids = folders, project = "6.NSCLC_Smith")
table(sce.big$orig.ident)
saveRDS(sce.big,file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs/","dataset1_GSE176021_NSCLC_CD8T.rds"))

## 6.NSCLC_smith part 
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor/GEX/part/")
folders <- list.files("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor/GEX/part/")
folders
sceList <- lapply(folders, function(folder){
  CreateSeuratObject(counts = Read10X(folder), project = folder,
                     min.cells = 3, min.features = 200)
})
sce.big <- merge(sceList[[1]], y = sceList[2:length(sceList)],
                 add.cell.ids = folders, project = "6.NSCLC_Smith")
table(sce.big$orig.ident)
saveRDS(sce.big,file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs/","dataset1_GSE176021_NSCLC_CD8T.rds"))

## 6.NSCLC_Smith TCR
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor/TCR")
samples <- list.files("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor/TCR")
samples
TcrList <- NULL
TcrList <- lapply(samples, FUN = function(sample){
  df <- read.csv(file.path(sample, "filtered_contig_annotations.csv"), header = T)
  df$sample <- sample
  df$new_barcode <- paste(df$sample, df$barcode, sep = "_")
  df
})
tcr.big <- TcrList[[1]]
for(i in 2:length(TcrList)){
  tcr.big <- rbind(TcrList[[i]], tcr.big)
}
unique(tcr.big$sample)
unique(tcr.big$high_confidence)
tcr.filtered <- tcr.big %>% 
  group_by(barcode,chain) %>%
  filter(umis == max(umis)) %>%
  ungroup()
tcr.big
write.csv(tcr.big, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor", "6.NSCLC_Smith_all_filtered_contig_annotations.csv"),
          quote = FALSE, row.names = FALSE)

## 6.NSCLC_Smith TCR
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor/TCR/part")
samples <- list.files("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor/TCR/part")
samples
TcrList <- NULL
TcrList <- lapply(samples, FUN = function(sample){
  df <- read.csv(file.path(sample, "filtered_contig_annotations.csv"), header = T)
  df$sample <- sample
  df$new_barcode <- paste(df$sample, df$barcode, sep = "_")
  df
})
tcr.big <- TcrList[[1]]
for(i in 2:length(TcrList)){
  tcr.big <- rbind(TcrList[[i]], tcr.big)
}
unique(tcr.big$sample)
filter_bool <- (tcr.big$is_cell == "True" & tcr.big$high_confidence == "True" & tcr.big$full_length == "True" & tcr.big$productive == "True")
tcr.big <- tcr.big[filter_bool,]
unique(tcr.big$high_confidence)
unique(tcr.big$is_cell)
unique(tcr.big$full_length)
unique(tcr.big$productive)
tcr.filtered <- tcr.big %>% 
  group_by(barcode,chain) %>%
  filter(umis == max(umis)) %>%
  ungroup()
tcr.big
write.csv(tcr.big, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor", "6.NSCLC_Smith_all_filtered_contig_annotations.csv"),
          quote = FALSE, row.names = FALSE)

## Tumor specific
library(readxl)
df <- read.csv("/picb/lilab5/liuzipei/scRNAseq_2022/TCR/6.NSCLC_Smith_all_filtered_contig_annotations.csv")
ts <- read_excel("/picb/lilab5/liuzipei/scRNAseq_2022/TCR/6.NSCLC_Smith_tcr_class.xlsx")
ts_bool <- ts$`Type of antigen recognized (viral or MANA)` == "MANA"
vir_bool <- ts$`Type of antigen recognized (viral or MANA)` != "MANA"
ts_MANA <- ts[ts_bool,]
ts_vir <- ts[vir_bool,]
ts_MANA$`TCR clonotype (AA)`
ts_vir
df1 <- df[,c("new_barcode","cdr3")]
df1
MANA.data <- merge(ts_MANA, df1, by.x = "TCR clonotype (AA)", by.y = "cdr3")
MANA.data
write.table(MANA.data, file = "/picb/lilab5/liuzipei/scRNAseq_2022/TCR/6.NSCLC_Smith_MANA.txt", sep = "\t", quote = FALSE, row.names = FALSE)
vir.data <- merge(ts_vir, df1, by.x = "TCR clonotype (AA)", by.y = "cdr3")
vir.data
write.table(vir.data, file = "/picb/lilab5/liuzipei/scRNAseq_2022/TCR/6.NSCLC_Smith_Virus.txt", sep = "\t", quote = FALSE, row.names = FALSE)



## 7.HNSCC_Ahmed
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/GEX")
fs <- list.files("./", pattern = "^GSM")
fs
samples <- str_split(fs,"_",simplify = T)[,1]
samples
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  y_len=length(str_split(y[1],'_',simplify = T))
  folder=paste(str_split(y[1],'_',simplify = T)[,2:(y_len-1)], collapse = "_")
  dir.create(folder,recursive = T)
  file.rename(y[1],file.path(folder,"barcodes.tsv.gz"))
  file.rename(y[2],file.path(folder,"features.tsv.gz"))
  file.rename(y[3],file.path(folder,"matrix.mtx.gz"))
})

folders <- list.files("/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/GEX")
folders

aList <- lapply(folders, function(folder){Read10X(folder)})
aList[[5]] <- aList[[5]]$`Gene Expression`
bList <- NULL
for(i in folders){
  bList[[i]] <- Read10X(i)
}
bList[["HPV42_TU_QVD"]] <- bList[["HPV42_TU_QVD"]]$`Gene Expression`
sceList <- lapply(folders, function(folder){
  CreateSeuratObject(counts = bList[[folder]], project = folder,
                     min.cells = 3, min.features = 200)
})
sceList
sce.big <- merge(sceList[[1]], y = sceList[2:length(sceList)],
                 add.cell.ids = folders, project = "7.HNSCC_Ahmed")
sce.big@assays$RNA@counts[,20]
table(sce.big$orig.ident)
saveRDS(sce.big,file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs/","dataset2_GSE180268_HNSCC_CD8T.rds"))

## 7.HNSCC_Ahmed TCR
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/TCR/")
fs <- list.files("./", pattern = "^GSM")
fs
samples <- str_split(fs,"_",simplify = T)[,1]
samples
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  y_len=length(str_split(y[1],'_',simplify = T))
  folder=paste(str_split(y[1],'_',simplify = T)[,2:4], collapse = "_")
  dir.create(folder,recursive = T)
  file.rename(y[1],file.path(folder,"barcodes.tsv.gz"))
  file.rename(y[2],file.path(folder,"features.tsv.gz"))
  file.rename(y[3],file.path(folder,"matrix.mtx.gz"))
})
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/TCR")
samples <- list.files("/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/TCR")
samples
data <- read.csv(samples[1])
for (file_name in samples[-1]) {
  data <- rbind(data, read.csv(file_name))
}
unique(data$new_barcode)
write.csv(data, file = paste("7.HNSCC_Ahmed", "all_filtered_contig_annotations.csv", sep = "_"), row.names = FALSE)
samples <- sapply(strsplit(samples, "_"), function(x){paste(x[1:3], collapse = "_")})
samples
for(sample in samples){
  df <- read.csv(sample, header = T)
  sample <- paste(strsplit(sample,"_")[[1]][1:2], collapse = "_")
  df$sample <- sample
  df$new_barcode <- paste(sample,df$barcode, sep = "_")
  write.csv(df, file = paste(sample, "all_filtered_contig_annotations.csv", sep = "_"), row.names = FALSE)
}

for(sample in samples){
  df <- read.csv(sample, header = T)
  df <- subset(df, df$high_confidence == "True" & df$full_length == "True" & df$productive == "True" & df$chain %in% c("TRA","TRB"))
  sample <- paste(strsplit(sample,"_")[[1]][2:4], collapse = "_")
  write.csv(df, file = paste(sample, "all_filtered_contig_annotations.csv", sep = "_"), row.names = FALSE)
}

write.csv()
sample
paste(str_split("GSM5456908_HPV15_QVD_Tumor_all_contig_annotations.csv","_")[1], collapse = "_")



TcrList <- NULL
TcrList <- lapply(samples, FUN = function(sample){
  df <- read.csv(file.path(sample, "filtered_contig_annotations.csv"), header = T)
  df$sample <- sample
  df$new_barcode <- paste(df$sample, df$barcode, sep = "_")
  df
})
tcr.big <- TcrList[[1]]
for(i in 2:length(TcrList)){
  tcr.big <- rbind(TcrList[[i]], tcr.big)
}
unique(tcr.big$sample)
tcr.big
write.csv(tcr.big, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/6.NSCLC_Smith/tumor", "6.NSCLC_Smith_all_filtered_contig_annotations.csv"),
          quote = FALSE, row.names = FALSE)



## PanCancer_zhang
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/PanCancer_zhang/10×")
samples <- list.files("./")
samples
sceList <- lapply(samples,function(sample){
  folder=paste(str_split(sample,'_',simplify = T)[,2], collapse = "")
  mtx <- read.table(sample, header = T)
  CreateSeuratObject(counts = mtx, project = folder,
                     min.cells = 3, min.features = 200)
})
folders <- NULL
for(i in samples){
  a <- paste(str_split(i,'_',simplify = T)[,2], collapse = "")
  folders <- c(folders,a)
}
samples
folders
metadata.info <- read.table("../GSE156728_metadata.txt", header = T)
head(metadata.info)
row.names(metadata.info) <- metadata.info$cellID
sceList <- lapply(sceList, function(x){
  AddMetaData(x, metadata = metadata.info, col.name = NULL)
})
sce.big <- merge(sceList[[1]], y = sceList[2:length(sceList)],
                 add.cell.ids = folders, project = "PanCancer_zhang")
sce.list <- SplitObject(sce.big, split.by = "patient")
sce.list <- lapply(sce.list,function(sce){
  mtx <- sce@assays$RNA
  CreateSeuratObject(counts = mtx, project = sce$patient,
                     min.cells = 3, min.features = 200)
})
sce.big2 <- merge(sce.list[[1]], y = sce.list[2:length(sce.list)], project = "PanCancer_zhang")
head(sce.big$cancerType)
tail(sce.big$cancerType)
table(sce.big$cancerType)
table(sce.big2$cancerType)
saveRDS(sce.big, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs/","dataset3_GSE156728_PanCancer_CD8T.rds"))

setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/PanCancer_zhang/")
df <- read.table("./GSE156728_10X_VDJ.merge.txt", sep = "\t", header = T)
df <- df[df$high_confidence == TRUE & df$is_cell == T & df$full_length == T & df$productive == "True", ]
df
df1 <- df %>% group_by(barcode, chain, library.id) %>%
  filter(umis == max(umis))
write.csv(df1, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/PanCancer_zhang/", "PanCancer_filtered_contig_annotations.csv"),
          quote = FALSE, row.names = FALSE)

## 9.NSCLC_zhang
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/9.NSCLC_zhang/GEX")
mtx <- readRDS("/picb/lilab5/liuzipei/scRNAseq_2022/data/9.NSCLC_zhang/GEX/GSE179994_all.Tcell.rawCounts.rds")
sce.big <- CreateSeuratObject(counts = mtx, project = "9.NSCLC_zhang", min.cells = 3, min.features = 200)
metadata.info <- read.table("./GSE179994_Tcell.metadata.tsv", header = T)
row.names(metadata.info) <- metadata.info$cellid
sce.big <- AddMetaData(sce.big, metadata = metadata.info, col.name = NULL)
sce.big$celltype
sce.cd8 <- subset(sce.big, celltype == "CD8")
sce.big@assays[["RNA"]]@counts
saveRDS(sce.cd8, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs/","dataset3_GSE179994_NSCLC_CD8T.rds"))

##Rosenberg
setwd("/picb/lilab5/liuzipei/Rosenberg_2022/rename")
samples <- list.files("/picb/lilab5/liuzipei/Rosenberg_2022/rename")
samples
sample = "FrTu4261"
sceList <- lapply(samples,function(sample){
  CreateSeuratObject(counts = Read10X(file.path(sample,"outs/filtered_feature_bc_matrix")),
                     project = sample, min.cells = 3, min.features = 200)
})
sce.big <- merge(sceList[[1]], y = sceList[2:length(sceList)],
                 add.cell.ids = samples, project = "Rosenberg")
table(sce.big$orig.ident)
saveRDS(sce.big, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/data/seuobjs/","dataset4_Rosenberg_allT.rds"))


##Rosenberg TCR
setwd("/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/TCR/")
samples <- list.files("/picb/lilab5/liuzipei/scRNAseq_2022/data/7.HNSCC_Ahmed/TCR/")
str_split(samples,"_",simplify = T)[,c(2,3,4)]
samples <- str_split(fs,"_",simplify = T)[,1]
TcrList <- list()
i = 1
for(sample in samples){
  TcrList[[i]] <- read.csv(file.path(sample), header = T)
  TcrList[[i]]$sample <- sample
  i = i + 1
}
tcr.big <- TcrList[[1]]
for(i in 2:length(TcrList)){
  tcr.big <- rbind(TcrList[[i]], tcr.big)
}
unique(tcr.big$sample.)
tcr.big$new.barcode <- paste(tcr.big$sample, tcr.big$barcode, sep = "_")
write.csv(tcr.big, file = file.path("/picb/lilab5/liuzipei/scRNAseq_2022/TCR", "all_filtered_contig_annotations.csv"),
          quote = FALSE, row.names = FALSE)

# i = 1
# for(sample in samples){
#   TcrList[[i]] <- read.csv(file.path(sample, "outs/filtered_contig_annotations.csv"), header = T)
#   TcrList[[i]]$sample <- sample 
#   i = i + 1
# }
samples
str_split(samples,'_',simplify = T)[,2:4]

paste(str_split(samples,'_',simplify = T)[,2:4], collapse = "_")
tcr.big <- TcrList[[1]]
for(i in 2:length(TcrList)){
  tcr.big <- rbind(TcrList[[i]], tcr.big)
}
unique(tcr.big$sample.)

tcr.big$new.barcode <- paste(tcr.big$sample, tcr.big$barcode, sep = "_")
unique(tcr.big$new.barcode)

write.csv(tcr.big, file = file.path("/picb/lilab5/liuzipei/Rosenberg_2022","all_filtered_contig_annotations.csv"),
          quote = FALSE, row.names = FALSE)



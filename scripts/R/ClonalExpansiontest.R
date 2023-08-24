suppressMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(rjson)
  library(rlist)
  library(tidyverse)
})

# 定义命令行参数
option_list = list(
  make_option(c("-t", "--tumorTCR"),
              type = "character",
              default = NULL,
              help = "rds file containing tumor TCR-seq results for a specific cancer type
              (Barcode,SampleID and a colmn to define clonetype must be included and match with scRNAseq data)",
              metavar = 'character'),
  make_option(c("-r", "--referenceTCR"),
              type = "character",
              default = NULL,
              help = "rds file containing reference-tissue TCR-seq results for a specific cancer type with the same information and matched SampleID and clonetype with tumorTCR",
              metavar = 'character'),
  make_option(c("-c", "--clonetype"),
              type = "character",
              default = NULL,
              help = "Name of the colnm to define a specific clonetype",
              metavar = 'character'),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output directory path for the results",
              metavar = 'character')
)

# 创建选项解析器
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# 检查必要参数是否提供
if (is.null(opt$tumorTCR) || is.null(opt$clonetype) || is.null(opt$out)) {
  print_help(opt_parser)
  stop("All required parameters must be provided", call. = FALSE)
}


tumorTCR <- readRDS(opt$tumorTCR)
referenceTCR <- readRDS(opt$referenceTCR)
clonetype<-opt$clonetype
names(tumorTCR)[names(tumorTCR) == clonetype] <- 'clonetype'
names(referenceTCR)[names(referenceTCR) == clonetype] <- 'clonetype'




referenceTCR<- referenceTCR %>% mutate(tissue='reference Tissue') %>%
  mutate(clonotype_id=paste0(tissue,'_',SampleID,'_',clonetype))

tumorTCR<- tumorTCR %>% mutate(tissue='Tumor') %>%
  mutate(clonotype_id=paste0(tissue,'_',SampleID,'_',clonetype))



calclonefreq<- function(clonetable){
  freq_clone_table <- clonetable %>% group_by(SampleID,tissue,clonotype_id) %>%
    summarise(SampleID,tissue,clonotype_id,CloneSize = n(), Frequency = n()/sum(n()),clonetype) %>%
    group_by(SampleID,tissue) %>%
    summarise(SampleID,tissue,clonotype_id,CloneSize, Frequency = CloneSize/sum(n()),clonetype) %>%
    unique()
  return(freq_clone_table)
}


referenceTCRfreq<-calclonefreq(referenceTCR)
tumorTCRfreq<-calclonefreq(tumorTCR)
matchedclone<-left_join(tumorTCRfreq,referenceTCRfreq,by=c('clonetype','SampleID')) %>%
  dplyr::select(SampleID,clonotype_id.x,CloneSize.x,Frequency.x,clonetype,clonotype_id.y,CloneSize.y,Frequency.y)
colnames(matchedclone)<-c('SampleID','clonotype_id','CloneSize','Frequency','clonetype','Ref.clonotype_id','Ref.CloneSize','Ref.Frequency')
matchedclone$Ref.CloneSize[is.na(matchedclone$Ref.CloneSize)] <- 0
matchedclone$Ref.Frequency[is.na(matchedclone$Ref.Frequency)] <- 0
tumorTCR<-left_join(tumorTCR,matchedclone,by=c('SampleID','clonotype_id','clonetype'))
saveRDS(tumorTCR,paste0(opt$out,'/tcr.matched.rds'))



























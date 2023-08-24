suppressMessages({
  library(GSEABase)
  library(tidyverse)
  library(optparse)
  library(AUCell)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(rjson)
  library(rlist)
})

# 定义命令行参数
option_list = list(
  make_option(c("-d", "--data"),
              type = "character",
              default = NULL,
              help = "rds file containing Seurat object for a specific cancer type.",
              metavar = 'character'),
  make_option(c("-s", "--signature"),
              type = "character",
              default = NULL,
              help = ".txt file containing the signature used for exhaustion score.",
              metavar = 'character'),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output directory path for results",
              metavar = 'character')
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


if (is.null(opt$data) || is.null(opt$signature)|| is.null(opt$out))  {
  print_help(opt_parser)
  stop("All required parameters must be provided", call. = FALSE)
}


seuratObj <- readRDS(opt$data)
seuList <- SplitObject(seuratObj, split.by = "dataset")
signatures <- read.table(opt$signature)
exh.geneset <- GeneSet(exhaustion.genes, setName="exh.genes")
if(length(seuList)>1){
  new.exList <- lapply(seuList, function(sce){
    exp.mtx <- GetAssayData(sce, assay = "RNA", slot = "count")
    cells_rankings <- AUCell_buildRankings(exp.mtx,plotStats=TRUE,splitByBlocks=TRUE,nCores=20)
    cells_rankings
    cells_AUC <- AUCell_calcAUC(exh.geneset, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05,nCores=20)
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assignCells =TRUE,nCores=20)
    ex.Assigned <- cells_assignment$exh.geneset$assignment
    return(ex.Assigned)
    })
  exh<-list()
  for(i in seq(seuList)){
    exh[[i]] <- subset(seuList[[i]], cells = new.exList[[i]])
    seuList[[i]]@meta.data<-seuList[[i]]@meta.data %>% mutate(barcode=rownames(.)) %>%
      mutate(isexhaustion=ifelse(barcode %in% rownames(exh[[i]]@meta.data),TRUE,FALSE))
  }
  sce.merge <- merge(seuList[[1]], y = seuList[2:length(seuList)])
}else{
  exp.mtx <- GetAssayData(seuratObj, assay = "RNA", slot = "count")
  cells_rankings <- AUCell_buildRankings(exp.mtx,plotStats=TRUE,splitByBlocks=TRUE,nCores=20)
  cells_rankings
  cells_AUC <- AUCell_calcAUC(exh.geneset, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05,nCores=20)
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assignCells =TRUE,nCores=20)
  ex.Assigned <- cells_assignment$exh.gene$assignment
  exh.cell <- subset(seuratObj, cells = ex.Assigned)
  seuratObj@meta.data<-seuratObj@meta.data %>% mutate(barcode=rownames(.)) %>%
    mutate(isexhaustion=ifelse(barcode %in% rownames(exh.cell@meta.data),TRUE,FALSE))
  sce.merge <- seuratObj
}
out<-paste0(opt$out,'/sce.rds')
saveRDS(sce.merge,out)




















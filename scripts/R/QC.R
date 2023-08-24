suppressMessages({
  library(scater)
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(rjson)
  library(rlist)
  library(tidyverse)
})

option_list = list(
  make_option(c("-d", "--data"),
              type = "character",
              default = NULL,
              help = "rds file containing Seurat object for a specific cancer type",
              metavar = 'character'),
  make_option(c("-q", "--thresholds"),
              type = "character",
              default = NULL,
              help = "Quality control standards: comma-separated values for
              'maxnCount_RNA,minnCount_RNA,maxnFeature_RNA,minnFeature_RNA,maxpctMT,maxpctRB,minpctRB,minpctTCR,minpctFeature'",
              metavar = 'character'),
  make_option("--dynamic",
              default = FALSE,
              action = "store_true",
              help = "Use dynamic thresholds for quality control"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output directory path for QC results",
              metavar = 'character')
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data) || is.null(opt$out)) {
  print_help(opt_parser)
  stop("All required parameters must be provided", call. = FALSE)
}


if (opt$dynamic!=TRUE){
  qc_values = strsplit(opt$thresholds, ",")[[1]]
  if (length(qc_values) != 9) {
    stop("Incorrect number of quality control standards provided", call. = FALSE)
  }
  qc_values <- as.numeric(qc_values)
}



out<-opt$out
if(! file.exists(out) | !file.exists(paste0(out,'/QCplot'))){
  dir.create(out)
  dir.create(paste0(out,'/QCplot'))
  }

seuratObj <- readRDS(opt$data)
seuList <- SplitObject(seuratObj, split.by = "dataset")
for(sa in names(seuList)){
  dirOutQC <-paste0(out,'/QCplot','/',sa)
  if(! file.exists(dirOutQC)){
    dir.create(dirOutQC)}

  seuList[[sa]]$pctMT <- PercentageFeatureSet(
    object = seuList[[sa]], pattern = "^MT-"
  )

  seuList[[sa]]$pctRB <- PercentageFeatureSet(
    object = seuList[[sa]], pattern = "^RP[SL]"
  )

  seuList[[sa]]$pctTCR <- PercentageFeatureSet(
    object = seuList[[sa]], pattern = "^TR[ABGD][VCJD]"
  )
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  seuList[[sa]] <- CellCycleScoring(seuList[[sa]],s.features = s.genes, g2m.features = g2m.genes)
  seuList[[sa]]$CC.Difference <- seuList[[sa]]$S.Score - seuList[[sa]]$G2M.Score
  seuList[[sa]]$pctFeature <- seuList[[sa]]$nCount_RNA/seuList[[sa]]$nFeature_RNA
  metaData.big <- seuList[[sa]]@meta.data
  show.nCount <- ggplot(metaData.big, aes(x = sampleID, y = nCount_RNA, fill = sampleID)) +
    geom_boxplot() +
    theme_bw() +
    ggtitle('nCount_RNA') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

  show.nFeature <- ggplot(metaData.big, aes(x = sampleID, y = nFeature_RNA, fill = sampleID)) +
    geom_boxplot() +
    theme_bw() +
    ggtitle('nFeature_RNA') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

  show.pctMT <- ggplot(metaData.big, aes(x = sampleID, y = pctMT, fill = sampleID)) +
    geom_boxplot() +
    theme_bw() +
    ggtitle('pctMT') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

  show.pctRB <- ggplot(metaData.big, aes(x = sampleID, y = pctRB, fill = sampleID)) +
    geom_boxplot() +
    theme_bw() +
    ggtitle('pctRB') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

  show.pctTCR <- ggplot(metaData.big, aes(x = sampleID, y = pctTCR, fill = sampleID)) +
    geom_boxplot() +
    theme_bw() +
    ggtitle('pctTCR') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

  show.pctFeature <- ggplot(metaData.big, aes(x = sampleID, y = pctFeature, fill = sampleID)) +
    geom_boxplot() +
    theme_bw() +
    ggtitle('pctFeature') +
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

  ggsave(basicQC, file = file.path(dirOutQC, 'basicQC.pdf'),
         dpi = 300, width = 13.5, height = 5)
  ggsave(advancedQC, file = file.path(dirOutQC, 'advancedQC.pdf'),
         dpi = 300, width = 10, height = 17.5)
  if (opt$dynamic==TRUE){
    attributes <- c("nCount_RNA", "nFeature_RNA", "pctMT", "pctRB", "pctTCR",'pctFeature')
    thresholds <- list()

    for (attribute in attributes) {
      outlier <- isOutlier(seuList[[sa]]@meta.data[[attribute]], nmads = 3, type = c("both"))
      thresholds[[attribute]] <- attr(outlier, "thresholds")
    }


    seuList[[sa]] <- subset(seuList[[sa]],subset = pctMT <= thresholds$pctMT[[2]] &
                              pctRB <= thresholds$pctRB[[2]] & pctRB >= thresholds$pctRB[[1]] &
                              pctTCR <= thresholds$pctTCR[[2]] & pctTCR >= thresholds$pctTCR[[1]] &
                              nCount_RNA<thresholds$nCount_RNA[[2]] & nCount_RNA>thresholds$nCount_RNA[[1]] &
                              nFeature_RNA<thresholds$nFeature_RNA[[2]] & nFeature_RNA>thresholds$nFeature_RNA[[1]] &
                              pctFeature<thresholds$pctFeature[[2]])
  }else{
    seuList[[sa]] <- subset(seuList[[sa]],subset = pctMT <= qc_values[[5]] &
                              pctRB <= qc_values[[7]] & pctRB >= qc_values[[6]] &
                              pctTCR >= qc_values[[8]] &
                              nCount_RNA<qc_values[[1]] & nCount_RNA>qc_values[[2]] &
                              nFeature_RNA<qc_values[[3]] & nFeature_RNA>qc_values[[4]] &
                              pctFeature<qc_values[[9]])

  }

}

if(length(seuList)>1){
  sce.merge <- merge(seuList[[1]], y = seuList[2:length(seuList)])
}else{
  sce.merge <- seuList[[1]]
}

saveRDS(sce.merge,paste0(out,'/sce.rds'))









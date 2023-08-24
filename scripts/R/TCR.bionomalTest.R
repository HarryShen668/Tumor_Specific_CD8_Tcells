suppressMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(rjson)
  library(rlist)
  library(dplyr)
  library(tidyverse)
})

# 定义命令行参数
option_list = list(
  make_option(c("-t", "--tumorTCR"),
              type = "character",
              default = NULL,
              help = "rds file containing tumor TCR-seq results for a specific cancer type generate by ClonalExpansiontest.R",
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
if (is.null(opt$tumorTCR) || is.null(opt$out)) {
  print_help(opt_parser)
  stop("All required parameters must be provided", call. = FALSE)
}


matched.TCR <-readRDS(opt$tumorTCR)

binomial.test <- function(matched.TCR) {
  TCR.Test <- matched.TCR %>% dplyr::select(SampleID,
                                            clonetype,
                                            CloneSize,
                                            Ref.CloneSize,
                                            Ref.Frequency) %>%
    unique() %>%
    group_by(SampleID) %>%
    summarise(
      SampleID,
      clonetype,
      CloneSize,
      Ref.CloneSize,
      fc = CloneSize / Ref.CloneSize,
      sum.CloneSize.tumor = sum(CloneSize),
      sum.CloneSize.Ref = sum(Ref.CloneSize) / sum(Ref.Frequency),
      exp = sum.CloneSize.tumor / sum.CloneSize.Ref) %>%
    as.data.frame()

  TCR.Test$sum.CloneSize.Ref <- ifelse(is.nan(TCR.Test$sum.CloneSize.Ref), 0, TCR.Test$sum.CloneSize.Ref)
  TCR.Test$p.val <- sapply(
    1:nrow(TCR.Test),
    function(i) {
      x1 <- as.numeric(TCR.Test[i, "CloneSize"])
      x2 <- as.numeric(TCR.Test[i, "Ref.CloneSize"])
      y1 <- as.numeric(TCR.Test[i, "sum.CloneSize.tumor"])
      y2 <- as.numeric(TCR.Test[i, "sum.CloneSize.Ref"])
      binom.test(x = x1, n = x1 + x2, p = y1 / (y1 + y2), alternative = "two.sided")$p.value
    }
  )

  return(TCR.Test)
}


TCR.test<-binomial.test(matched.TCR)
matched.TCR<-left_join(matched.TCR,TCR.test,by=c('SampleID','clonetype','CloneSize','Ref.CloneSize'))
saveRDS(matched.TCR,paste0(opt$out,'/tcr.matched.tested.rds'))

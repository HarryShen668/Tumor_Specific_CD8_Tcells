library(Seurat)
library(tidyverse)

projection.plot <- function(ref, query = NULL, labels.col = "celltype", 
          cols = NULL, linesize = 1, pointsize = 1, ref.alpha = 0.3, 
          ref.size = NULL, title = "Projection of query on reference map") {
  labels <- ref[[labels.col]][, 1]
  states_all <- levels(factor(labels))
  nstates <- length(states_all)
  if (!is.null(cols)) {
    if (nstates > length(cols)) {
      warning("Not enough colors provided. Making an automatic palette")
      palette <- rainbow(n = nstates)
    }
    else {
      palette <- cols
    }
  }
  else {
    palette <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", 
                   "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", 
                   "#e812dd")
    if (nstates > length(palette)) {
        palette <- rainbow(n = nstates)
    }
  }
  cols_use <- scales::alpha(palette, ref.alpha)
  if (is.null(query)) {
    p <- DimPlot(ref, reduction = "umap", label = FALSE, 
                 group.by = labels.col, repel = TRUE, pt.size = ref.size, 
                 cols = cols_use) + ggtitle("Reference map") + theme(aspect.ratio = 1)
  }
  else {
    p <- DimPlot(ref, reduction = "umap", label = FALSE, 
                 group.by = labels.col, repel = TRUE, pt.size = ref.size, 
                 cols = cols_use) + geom_point(data.frame(query@reductions$umap@cell.embeddings), 
                                               mapping = aes(x = umap_1, y = umap_2), alpha = 0.6, 
                                               size = pointsize, shape = 17, color = "red") + 
      geom_density_2d(data = data.frame(query@reductions$umap@cell.embeddings),
                      mapping = aes(x = umap_1, y = umap_2), color = "black",
                      n = 200, h = 2, size = linesize) + ggtitle(title) +
      theme(aspect.ratio = 1)
  }
  return(p)
}
# color1 <- c('#d62e2d', '#e9787a', '#cd9c9b', '#684797', '#3377a9', '#96c3d8',
#   '#67a59b', '#70b866', '#6a9a52', '#a5d38f', '#ab9697', '#f19294',
#   '#f5b375', '#da8f6f', '#e0c880', '#f47d2f', '#ec8d63', '#e45d61',
#   '#a65a35', '#4b9d47')
# ref <- readRDS("/home/shenhaoyu/xinda/dataset/CD8TILreference/1.expr/seuList_v1_harmony_map.Rds")
# ref$celltype <- ref$meta.cluster
# ref <- RunHarmony(ref, group.by.vars = "cancerType", plot_convergence = TRUE, verbose = F,project.dim = FALSE)
# ref[['umap']] <- RunUMAP(Embeddings(ref, 'harmony')[, 1:30],
#                           assay='RNA', verbose=FALSE, umap.method='uwot', return.model=TRUE)
# DimPlot(ref, reduction = "umap", group.by = 'meta.cluster',cols=color1)
# 
# library(symphony)
# ref1 <- buildReferenceFromSeur(ref, assay = 'SCT',
#                                  verbose=TRUE, save_umap=TRUE, save_uwot_path='/home/liuzipei/project/cache_symphony_sct.uwot')
# 
# 
# 
# query <- readRDS("/home/shenhaoyu/xinda/dataset/Rosenberg_2022/2.mapping/query.map.cd8.Rds")
# ref$meta.cluster
# 
# ref$celltype
# ref[["celltype"]]
# ref[["umap"]]
# DimPlot(ref, group.by = "meta.cluster", cols = color1)
# DimPlot(query, group.by = "celltype", cols = color1)
# projection.plot(ref, query = query, cols = color1, linesize = 0.5, pointsize = 0.5)




library(Seurat)
library(SeuratData)
InstallData("pbmc3k")
pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")

levels(pbmc)

# Find differentially expressed features between CD14+ and FCGR3A+ Monocytes
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono")
# view results
head(monocyte.de.markers)
pbmc@assays$RNA@data

install.packages("tictoc")
library(tictoc)


tic()
pbmc.markers <- FindAllMarkers(pbmc)
toc()
#seurat_annotations

Seurat_top2 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

library(presto)
tic()
DE_list <- wilcoxauc(pbmc, 'seurat_annotations',assay = 'data')
toc()

presto_top2 <- top_markers(DE_list , n = 2)

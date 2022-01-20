all_samples_scATAC_Object <- readRDS("all_samples_scATAC_Object.rds")

library(tidyr)
library(wesanderson)
library(presto)
ATAC <- all_samples_scATAC_Object$ATAC
df <- as.data.frame(do.call(rbind, strsplit(row.names(ATAC@meta.data), '@')))
ATAC@meta.data[,c("sample", "hours", "cell")] <- df  %>%
  separate(V1,into = c("sample", "hours"),
           sep = "(?<=[A-Za-z])(?=[0-9])"
          )

DimPlot(ATAC, label = TRUE, reduction = "umap", group.by = "assign.ident", repel=T, pt.size = 0.5, label.size = 3) + labs(title = "Clustering of PBMCs from all samples")

DimPlot(ATAC, label = FALSE, reduction = "umap", group.by = "hours", repel=T, pt.size = 0.5, label.size = 2.5, cols = brewer.pal(3,"Set3"))

DimPlot(ATAC, label = FALSE, reduction = "umap", group.by = "hours", repel=T, pt.size = 0.5, label.size = 2.5, cols = brewer.pal(3,"Pastel2")) + labs(title = "Clustering of PBMCs in different time points")




X_matrix = GetAssayData(ATAC, slot = "data")
y <- Idents(ATAC) %>% unlist %>% as.character()
test.res = wilcoxauc(X_matrix, y)


remotes::install_github("immunogenomics/presto")

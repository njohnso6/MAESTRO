DimPlot(ATAC, label = TRUE, reduction = "umap", group.by = "assign.ident", repel=T, pt.size = 0.3, label.size = 3) + labs(title = "Umap showing annotation of 3 CLL and 3 Healthy Donors")

DimPlot(ATAC, label = FALSE, reduction = "umap", group.by = "tumor",repel=T, pt.size = 0.3, label.size = 2.5,cols = brewer.pal(2,"Dark2")) + labs(title = "Umap grouped by samples")


ATAC@meta.data <- ATAC@meta.data %>% mutate( tumor = ifelse(sample == "PBMC", "Healthy Donors", "CLL"))
ATAC@meta.data <- ATAC@meta.data %>% mutate( new_hours = ifelse(hours== "0", "0", "8+"))
ATAC@meta.data <- ATAC@meta.data %>% mutate( group = str_c(tumor , new_hours, sep="_"))


DimPlot(ATAC, label = FALSE, reduction = "umap", group.by = "group",repel=T, pt.size = 0.3, label.size = 2.5) + labs(title = "Umap grouped by samples")


DimPlot(object = ATAC, label = TRUE, pt.size = 0.15,
             group.by = "assign.ident", label.size = 3,
             cols = brewer.pal(8,"Set2")) +
  labs(title = "Umap showing annotation of PBMC samples from 3 CLL and 3 Healthy Donors") +
  theme_minimal() + NoLegend()

ggsave("output/UMAP_annot.png", width=6, height=5)


DimPlot(object = ATAC, label = FALSE, pt.size = 0.15,
        group.by = "group", label.size = 3,
        cols = brewer.pal(4,"Set2")) +
  labs(title = "Umap grouped by samples") +
  theme_minimal()

ggsave("output/UMAP_bysamples.png", width=9, height=5)

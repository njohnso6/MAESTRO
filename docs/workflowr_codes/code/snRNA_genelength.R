library(dplyr)

snRNA <- readRDS("./data/snrna/mm_LT_750_scRNA_Object.rds")
snRNA <- snRNA$RNA
genes <- rownames(snRNA@assays$RNA@counts)

library(EnsDb.Hsapiens.v86)

edb <- EnsDb.Hsapiens.v86
genes_ensemble <- genes(edb)
gene_length<- as.data.frame(genes_ensemble) %>% dplyr::select(gene_id, gene_name, width)


genes <- as.data.frame(genes)
colnames(genes) <- "gene_name"

genes_length <- inner_join(genes,gene_length)
genes <- genes_length$gene_name

raw_count <- snRNA@assays$RNA@data
count <- snRNA@assays$RNA@data[genes,]
dim(count)

cells <- colnames(count)
metadata <- snRNA@meta.data[cells,]
celltype <- metadata$seurat_clusters
names(celltype) <- rownames(metadata)
head(celltype)
table(celltype)


###On data slot
mat.sparse <- count
mat <- Matrix(count,sparse = FALSE)

##cell id
getPseudobulk <- function(mat, celltype) {
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
    cells <- names(celltype)[celltype==ct]
    pseudobulk <- rowSums(mat.sparse[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}

## test runtime
start_time1 <- Sys.time()
## call function
mat.summary <- getPseudobulk(mat, celltype)
end_time1 <- Sys.time()

## take a look
dim(mat.summary)
head(mat.summary)

sum <- as.data.frame(mat.summary)

sum$gene_name <- rownames(sum)
dim(sum)
head(sum)
pseudobulk_data <- inner_join(sum, gene_length)
head(gene_length)

#Plot
library(patchwork)
library(ggplot2)

p0 <- ggplot(data = pseudobulk_data, mapping = aes(x = pseudobulk_data$width, y = pseudobulk_data$`0`)) +
  geom_point(alpha = 0.1) +
  ylab("Cluster 0") +
  stat_cor(label.x = 3)



p2 <- ggplot(data = pseudobulk_data, mapping = aes(x = pseudobulk_data$width, y = pseudobulk_data$`2`)) +
  geom_point(alpha = 0.1) +
  ylab("Cluster 2") +
  stat_cor(label.x = 3)

p5 <- ggplot(data = pseudobulk_data, mapping = aes(x = pseudobulk_data$width, y = pseudobulk_data$`5`)) +
  geom_point(alpha = 0.1) +
  ylab("Cluster 5") +
  stat_cor(label.x = 3)


p0 | p2 | p5

ggsave("output/snRNA_genelengthVSLognormUMI.png", width = 12, height = 4)
boxplot(sum$`0`)

sum <- sum %>% filter(`0` < 40000)

####Two step normalization




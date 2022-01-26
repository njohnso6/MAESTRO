#1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#3. Divide the RPK values by the “per million” scaling factor. This gives you TPM.

#1
raw <- snRNA@assays$RNA@counts
raw_count <- snRNA@assays$RNA@counts[genes,]
dim(raw_count)
mat_count <- as.data.frame(Matrix(count,sparse = FALSE))
mat_count[1:5,1:5]
dim(mat_count)
max(mat_count)

genes_length <- genes_length %>% mutate(width_kb = width/1000)
dim(genes_length)
head(genes_length)

by_genelength <- function(x) {
  vector <- x/genes_length$width_kb
  return(vector)
}

dim(mat_count)
norm_genelength <- apply(mat_count, 2, by_genelength)
dim(norm_genelength)
max(norm_genelength)

#2
colsum <- apply(norm_genelength, 2, sum)
scaling_factor <- colsum/1000000

#3
by_scalefact <- function(x) {
  vector <- x/scaling_factor
  return(vector)
}

length(scaling_factor)
dim(norm_genelength)

TPM <- apply(norm_genelength, 1, by_scalefact)
max(TPM)
min(TPM)

#####
#The same as previous
#####

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
TPM_pseudo <- getPseudobulk(TPM, celltype)
dim(TPM_pseudo)
end_time1 <- Sys.time()

TPM_pseudo <- as.data.frame(TPM_pseudo)
head(TPM_pseudo)
dim(TPM_pseudo)
TPM_pseudo$gene_name <- rownames(TPM_pseudo)
dim(TPM_pseudo)
TPM_pseudo %>% head()

head(genes_length)
TPM_pseudo$gene_length <- genes_length$width
#Plot
library(patchwork)
library(ggplot2)

Cluster0 <- ggplot(data = TPM_pseudo, mapping = aes(x = TPM_pseudo$gene_length, y = TPM_pseudo$`0`)) +
  geom_point(alpha = 0.1) +
  ylab("Cluster 0") +
  stat_cor(label.x = 3)

Cluster2 <- ggplot(data = TPM_pseudo, mapping = aes(x = TPM_pseudo$gene_length, y = TPM_pseudo$`2`)) +
  geom_point(alpha = 0.1) +
  ylab("Cluster 2") +
  stat_cor(label.x = 3)

Cluster5 <- ggplot(data = TPM_pseudo, mapping = aes(x = TPM_pseudo$gene_length, y = TPM_pseudo$`5`)) +
  geom_point(alpha = 0.1) +
  ylab("Cluster 5") +
  stat_cor(label.x = 3)

Cluster0 | Cluster2 | Cluster5

ggsave("output/snRNA_TPM_VS_genelength.png", height = 6, width = 15)



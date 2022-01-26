library(ExperimentHub)
eh <- ExperimentHub()
query(eh, "Kang")

sce <- eh[["EH2259"]]

sce <- sce[rowSums(counts(sce) > 0) > 0, ]
counts(sce)
dim(sce)

library(scater)

qc <- perCellQCMetrics(sce)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)


sce$id <- paste0(sce$stim, sce$ind)
(sce <- prepSCE(sce,
                kid = "cell", # subpopulation assignments
                gid = "stim",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

sce <- runUMAP(sce, pca = 20)
.plot_dr <- function(sce, dr, col)
  plotReducedDim(sce, dimred = dr, colour_by = col) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_minimal() + theme(aspect.ratio = 1)

cs_by_k <- split(colnames(sce), sce$cluster_id)
cs100 <- unlist(sapply(cs_by_k, function(u)
  sample(u, min(length(u), 100))))

# plot t-SNE & UMAP colored by cluster & group ID
for (dr in c("TSNE", "UMAP"))
  for (col in c("cluster_id", "group_id"))
    .plot_dr(sce[, cs100], dr, col)


###DS
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)
t(head(assay(pb)))


# run DS analysis
res <- pbDS(pb, verbose = FALSE)
# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)
## [1] "B cells


k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))

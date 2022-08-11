suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(ggplot2))
require(reshape2)
require(dplyr)
library(doParallel)
library(foreach)

data_dir <- "/project2/xinhe/kevinluo/GSFA/data/"
res_dir <- "/project2/xinhe/kevinluo/GSFA/twostep_clustering/Stimulated_T_Cells/stimulated"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
setwd(res_dir)

n_cores <- 8

# 0. Load input data ####
cat("Load input data ... \n")
combined_obj <- readRDS(file.path(data_dir,"TCells_cropseq_data_seurat.rds"))
feature.names <- data.frame(fread(file.path(data_dir, "Stimulated_T_Cells_GSE119450_RAW_D1N_genes.tsv.gz"),
                                  header = FALSE), stringsAsFactors = FALSE)
metadata <- combined_obj@meta.data
table(metadata$orig.ident)

# Extract data for stimulated cells
combined_obj@meta.data$condition <- "unstimulated"
combined_obj@meta.data$condition[which(endsWith(combined_obj@meta.data$orig.ident, "S"))] <- "stimulated"

# combined_obj.list <- SplitObject(combined_obj, split.by = "condition")
# combined_obj <- combined_obj.list$stimulated
combined_obj <- subset(combined_obj, subset = condition == "stimulated")
combined_obj
table(combined_obj@meta.data$orig.ident)

# 1. Data preprocessing ####
cat("Data preprocessing ... \n")

cat("QC ... \n")
# The number of unique genes detected in each cell.
range(combined_obj$nFeature_RNA)

# The total number of molecules detected within a cell
range(combined_obj$nCount_RNA)

# The percentage of reads that map to the mitochondrial genome
range(combined_obj$percent_mt)

cat("Data filtering ... \n")
# We filter cells that have more than 500 genes identified.
combined_obj <- subset(combined_obj, subset = nFeature_RNA > 500)

# Normalizing the data
combined_obj <- NormalizeData(combined_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
combined_obj <- FindVariableFeatures(combined_obj, selection.method = "vst", nfeatures = 1000)
combined_obj

# Regress out total UMI counts per cell and percent of mitochondrial genes detected per cell and scaled to obtain gene level z-scores.
combined_obj <- ScaleData(combined_obj, vars.to.regress = c("nCount_RNA", "percent_mt"))

saveRDS(combined_obj, file = file.path(res_dir, "TCells_stimulated_seurat_processed_data.rds"))

# 2. Perform linear dimensional reduction
combined_obj <- readRDS(file.path(res_dir, "TCells_stimulated_seurat_processed_data.rds"))
combined_obj <- RunPCA(combined_obj, features = VariableFeatures(object = combined_obj))
pdf(file.path(res_dir, "TCells_stimulated_seurat_pca.pdf"), width = 5, height = 5)
ElbowPlot(combined_obj, ndims = 50)
dev.off()

# 3. Run non-linear dimensional reduction (UMAP/tSNE) ####
combined_obj <- RunUMAP(combined_obj, dims = 1:30)

# 4. Cluster the cells ####
combined_obj <- FindNeighbors(combined_obj, dims = 1:30)
combined_obj <- FindClusters(combined_obj)
cluster_labels <- Idents(combined_obj)
levels(combined_obj)

# Look at cluster IDs of the first 5 cells
head(cluster_labels, 5)
saveRDS(combined_obj, file = file.path(res_dir, "TCells_stimulated_seurat_clustered.rds"))

# 5. Finding differentially expressed features (cluster biomarkers) ####
combined_obj <- readRDS(file.path(res_dir, "TCells_stimulated_seurat_clustered.rds"))
cat("Run DE test using MAST...\n")
cat(length(levels(combined_obj)), "clusters.\n")
registerDoParallel(cores=n_cores)
ptm <- proc.time()
de.markers <- foreach(i=levels(combined_obj), .packages="Seurat") %dopar% {
  FindMarkers(combined_obj, ident.1 = i, test.use = "MAST")
}
proc.time() - ptm
stopImplicitCluster()
saveRDS(de.markers, file = file.path(res_dir, "TCells_stimulated_seurat_MAST_DEGs.rds"))

# cat("Run DE test using MAST...\n")
# de.markers <- vector("list", length = length(levels(combined_obj)))
# names(de.markers) <- levels(combined_obj)
# for(i in levels(combined_obj)){
#   cat("cluster", i, "\n")
#   system.time(
#     de.markers[[i]] <- FindMarkers(combined_obj, ident.1 = i, test.use = "MAST"))
# }
# saveRDS(de.markers, file = file.path(res_dir, "TCells_stimulated_seurat_MAST_DEGs.rds"))

# cat("Run DE test using DESeq2...\n")
# registerDoParallel(cores=n_cores)
# ptm <- proc.time()
# de.markers <- foreach(i=levels(combined_obj), .packages="Seurat") %dopar% {
#   FindMarkers(combined_obj, ident.1 = i, test.use = "DESeq2")
# }
# names(de.markers) <- levels(combined_obj)
# proc.time() - ptm
# stopImplicitCluster()
# saveRDS(de.markers, file = file.path(res_dir, "TCells_stimulated_seurat_DESeq2_DEGs.rds"))


# cat("Run DE test using DESeq2...\n")
# de.markers <- vector("list", length = length(levels(combined_obj)))
# names(de.markers) <- levels(combined_obj)
# for(i in levels(combined_obj)){
#   cat("cluster", i, "\n")
#   system.time(
#     de.markers[[i]] <- FindMarkers(combined_obj, ident.1 = i, test.use = "DESeq2"))
# }
# saveRDS(de.markers, file = file.path(res_dir, "TCells_stimulated_seurat_DESeq2_DEGs.rds"))


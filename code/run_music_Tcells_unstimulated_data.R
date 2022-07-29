suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(MUSIC))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5),
                             axis.title = element_text(size = 14),
                             axis.text = element_text(size = 13),
                             legend.title = element_text(size = 13),
                             legend.text = element_text(size = 12),
                             panel.grid.minor = element_blank())
)

data_dir <- "/project2/xinhe/yifan/Factor_analysis/Stimulated_T_Cells/"
dir.create("/scratch/midway2/kaixuan/GSFA/music_analysis/Stimulated_T_Cells/unstimulated", recursive = TRUE, showWarnings = FALSE)
setwd("/scratch/midway2/kaixuan/GSFA/music_analysis/Stimulated_T_Cells/unstimulated")
dir.create("./music_output", recursive = TRUE, showWarnings = FALSE)
n_cpus <- 5

# 0. Load input data ####
cat("Load input data ... \n")

combined_obj <- readRDS('/project2/xinhe/yifan/Factor_analysis/shared_data/TCells_cropseq_data_seurat.rds')
metadata <- combined_obj@meta.data
metadata[1:5, ]

table(metadata$orig.ident)

# stimulated_cells <- rownames(metadata)[which(endsWith(metadata$orig.ident, "S"))]
# cat(length(stimulated_cells), "stimulated cells. \n")
unstimulated_cells <- rownames(metadata)[which(endsWith(metadata$orig.ident, "N"))]
cat(length(unstimulated_cells), "unstimulated cells. \n")

feature.names <- data.frame(fread(paste0(data_dir, "GSE119450_RAW/D1N/genes.tsv"),
                                  header = FALSE), stringsAsFactors = FALSE)

expression_profile <- combined_obj@assays$RNA@counts[, unstimulated_cells]
rownames(expression_profile) <- feature.names$V2[match(rownames(expression_profile),
                                                       feature.names$V1)]

cat("Dimension of expression profile for unstimulated cells: \n")
dim(expression_profile)

targets <- names(combined_obj@meta.data)[4:24]
targets[targets == "NonTarget"] <- "CTRL"
perturb_information <- apply(combined_obj@meta.data[unstimulated_cells, 4:24], 1,
                             function(x){ targets[which(x > 0)] })

# 1. Data preprocessing ####
cat("Data preprocessing ... \n")

cat("QC ... \n")
crop_seq_list <- Input_preprocess(expression_profile, perturb_information)
saveRDS(crop_seq_list, "music_output/music_crop_seq_list.rds")

x <- c(1:100)
saveRDS(x, "music_output/x.rds")

crop_seq_qc <- Cell_qc(crop_seq_list$expression_profile,
                       crop_seq_list$perturb_information,
                       species = "Hs", plot = F)
saveRDS(crop_seq_qc, "music_output/music_crop_seq_qc.rds")

cat("Data imputation ... \n")
crop_seq_imputation <- Data_imputation(crop_seq_qc$expression_profile,
                                       crop_seq_qc$perturb_information,
                                       cpu_num = n_cpus)
saveRDS(crop_seq_imputation, "music_output/music_imputation.merged.rds")

cat("Data filtering ... \n")
crop_seq_filtered <- Cell_filtering(crop_seq_imputation$expression_profile,
                                    crop_seq_imputation$perturb_information,
                                    cpu_num = n_cpus)
saveRDS(crop_seq_filtered, "music_output/music_filtered.merged.rds")

# 2. Model building ####
cat("Model building ... \n")

crop_seq_vargene <- Get_high_varGenes(crop_seq_filtered$expression_profile,
                                      crop_seq_filtered$perturb_information, plot = T)
saveRDS(crop_seq_vargene, "music_output/music_vargene.merged.rds")

## Get_topics() can take up to a few hours to finish,
## depending on the size of data
cat("Fitting model with 5 topics ... \n")

system.time(
  topic_1 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 5))
saveRDS(topic_1, "music_output/music_merged_5_topics.rds")

cat("Fitting model with 10 topics ... \n")
system.time(
  topic_2 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 10))
saveRDS(topic_2, "music_output/music_merged_10_topics.rds")

cat("Fitting model with 15 topics ... \n")
system.time(
  topic_3 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 15))
saveRDS(topic_3, "music_output/music_merged_15_topics.rds")

cat("Fitting model with 20 topics ... \n")
system.time(
  topic_4 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 20))
saveRDS(topic_4, "music_output/music_merged_20_topics.rds")

## try fewer numbers of topics
system.time(
  topic_5 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 4))
saveRDS(topic_5, "music_output/music_merged_4_topics.rds")

system.time(
  topic_6 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 6))
saveRDS(topic_6, "music_output/music_merged_6_topics.rds")

system.time(
  topic_7 <- Get_topics(crop_seq_vargene$expression_profile,
                        crop_seq_vargene$perturb_information,
                        topic_number = 7))
saveRDS(topic_7, "music_output/music_merged_7_topics.rds")

# 3. Pick optimal number of topics ####
cat("Select optimal number of topics... \n")

topic_model_list <- list()
topic_model_list$models <- list()
topic_model_list$perturb_information <- topic_1$perturb_information
topic_model_list$models[[1]] <- topic_1$models[[1]]
topic_model_list$models[[2]] <- topic_2$models[[1]]
topic_model_list$models[[3]] <- topic_3$models[[1]]
topic_model_list$models[[4]] <- topic_4$models[[1]]
topic_model_list$models[[5]] <- topic_5$models[[1]]
topic_model_list$models[[6]] <- topic_6$models[[1]]
topic_model_list$models[[7]] <- topic_7$models[[1]]
saveRDS(topic_model_list, "music_output/topic_model_list.rds")

optimalModel <- Select_topic_number(topic_model_list$models,
                                    plot = T,
                                    plot_path = "music_output/select_topic_number.pdf")
saveRDS(optimalModel, "music_output/optimalModel.rds")

library("AUCell")
library("CellChat")
library("Seurat")
library("sceasy")
library("reticulate")
library("dplyr")
library("SpatialDecon")
library("InSituType")
library("tidyverse")
#library("schard")
library("ggplot2")
library("ggpubr")
library("here")
library("dittoSeq")
library("spatstat")


###############################################################
#################### COMMAND LINE INPUT #######################
###############################################################

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) < 1) {
  stop("Please provide one or more file paths as command-line arguments.")
}

# Iterate over each file path and print it (or do something else with it)
for (file in args) {
  print(paste("Processing file:", file))
  # You can add your file processing code here
}

##################################################
#################### INPUT #######################
##################################################

sample_name = args[1]
print(paste("Input sample name: ", sample_name))

sample_rds_path = args[2]
print(paste("Input sample fully labeled rds path: ", sample_rds_path))

insitutype_reference_sig_csv_file_path = args[3]
print(paste("Input inSituType ref sig profile from Nanostring: ", insitutype_reference_sig_csv_file_path))

output_folder_path = args[4]
print(paste("Input output folder: ", output_folder_path))

#Section 1 - Prepare Data and Genesets for AUCell analysis
#sample_name = "PC429"
#load cosmx data
#seurat_obj_norm <- readRDS("/cristealab/xiwang/Outputs/CosMx_Spatial_Pipelines/output_PCA429_08162024/InSituType/PCA429_semi_sup_insitutype_fully_labeled.rds")
seurat_obj_norm <- readRDS(sample_rds_path)

#select only count data from tumor cells
malignant <- subset(seurat_obj_norm, subset = level_2_clusters == "Tumor")
malignant_counts <- GetAssayData(object=malignant, slot="counts")

#load basal and classical gene sets
#genesets <- read.csv("/cristealab/ajordan/projects/BTC_analysis/data/external/gene_sets/basal_classical_and_fibroblast_gene_sets.csv")
genesets <- read.csv(insitutype_reference_sig_csv_file_path)
#select gene sets that represent basal and classical PDAC cells
polarity_genesets <- genesets[genesets$gs_name == "Raghavan_scClassical"| genesets$gs_name=="Raghavan_scBasal" | genesets$gs_name == "Raghavan_scIC",]

#output_folder_path = "/cristealab/xiwang/Outputs/CosMx_Spatial_Pipelines/output_PCA429_08162024/polarity_and_fibroblast"


###################################################
#################### OUTPUT #######################
###################################################

plot_section_4_1_correlation_AUCell_scores_path = paste0(output_folder_path, "/4_1_plot_correlation_AUCell_scores.jpeg")
plot_section_4_2_correlation_classical_basal_scores_path = paste0(output_folder_path, "/4_2_plot_correlation_classical_basal_scores.jpeg")

plot_section_4_3_classical_scores_vs_subtype_path = paste0(output_folder_path, "/4_3_plot_classical_scores_vs_subtype.jpeg")
plot_section_4_4_basal_scores_vs_subtype_path = paste0(output_folder_path, "/4_4_plot_basal_scores_vs_subtype.jpeg")
plot_section_4_5_coexpressor_scores_vs_subtype_path = paste0(output_folder_path, "/4_5_plot_coexpressors_scores_vs_subtype.jpeg")


plot_section_6_1_boxplots_distribution_basal_classical = paste0(output_folder_path, "/6_1_boxplots_distribution_basal_classical.jpeg")
rds_path_full_cohort_semi_sup_insitutype_tumor_subtype_labels = paste0(output_folder_path, "/", sample_name, "_full_cohort_semi_sup_insitutype_tumor_subtype_labels.rds")


plot_section_8_1_boxplot_level_3_cluster_fibrablast_subtype_frequency = paste0(output_folder_path, "/8_1_boxplot_level_3_cluster_fibrablast_subtype_frequency.jpeg")
plot_section_8_2_boxplot_level_3B_cluster_fibrablast_subtype_frequency = paste0(output_folder_path, "/8_2_boxplot_level_3B_cluster_fibrablast_subtype_frequency.jpeg")

plot_section_8_3_boxplot_level_2_cluster_fibrablast_subtype_frequency_cafs = paste0(output_folder_path, "/8_3_boxplot_level_2_cluster_fibrablast_subtype_frequency_cafs.jpeg")
plot_section_8_4_boxplot_level_3B_cluster_fibrablast_subtype_frequency_cafs = paste0(output_folder_path, "/8_4_boxplot_level_3B_cluster_fibrablast_subtype_frequency_cafs.jpeg")

plot_section_8_5_fap_caf_vlnplot = paste0(output_folder_path, "/8_4_boxplot_level_3B_cluster_fibrablast_subtype_frequency_cafs.jpeg")

######################################################
#################### SECTION 2 #######################
######################################################
#Section 2 - run AUCell to score cells for Classical and Basal programs
set.seed(123)
malignant_rankings <- AUCell_buildRankings(malignant_counts, plotStats=TRUE)


#filter gene sets to include genes present in the rankings
polarity_genesets_cosmx <-  polarity_genesets[polarity_genesets$gene %in% rownames(malignant_rankings),]

#create list of gene sets
polarity_genesets_cosmx <- split(polarity_genesets_cosmx$gene, polarity_genesets_cosmx$name)

#run AUC
cells_AUC <- AUCell_calcAUC(polarity_genesets_cosmx, malignant_rankings)
#save(cells_AUC, file="/cristealab/xiwang/Outputs/CosMx_Spatial_Pipelines/output_PCA429_08162024/polarity_and_fibroblast/cells_AUC.RData")
save(cells_AUC, file=file.path(output_folder_path, "cells_AUC.RData"))


#explore AUC thresholds 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 


######################################################
#################### SECTION 3 #######################
######################################################
#Section 3 - Explore Classical and Basal thresholds
#plot classical AUC threshold
geneSetName <- rownames(cells_AUC)[grep("Classical", rownames(cells_AUC))]

aucell_histplot = AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.25)
abline(v=0.25)

# save(aucell_histplot, file = "/cristealab/xiwang/Outputs/CosMx_Spatial_Pipelines/output_PCA429_08162024/polarity_and_fibroblast/3_plot_classical_AUC_threshold.jpeg")

auc_scores_classical <- getAUC(cells_AUC)[geneSetName, ]
# Filter cells with AUCell scores between 0.09 and 0.16

weakClassical <- names(which(auc_scores_classical >=0.02 & auc_scores_classical < 0.097))
Classical <- names(which(auc_scores_classical >= 0.097 & auc_scores_classical <= 0.167))
strongClassical <- names(which(getAUC(cells_AUC)[geneSetName,]>0.16))


geneSetName <- rownames(cells_AUC)[grep("Basal", rownames(cells_AUC))]

#plot basal auc threshold
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.25)
abline(v=0.25)

auc_scores_basal <- getAUC(cells_AUC)[geneSetName, ]
# Filter cells with AUCell scores between 0.09 and 0.16
weakBasal <- names(which(auc_scores_basal >=0.045 & auc_scores_basal < 0.093))
Basal <- names(which(auc_scores_basal >= 0.093 & auc_scores_basal <= 0.18))
strongBasal <- names(which(getAUC(cells_AUC)[geneSetName,]>0.18))

geneSetName <- rownames(cells_AUC)[grep("scIC", rownames(cells_AUC))]
auc_scores_IC <- getAUC(cells_AUC)[geneSetName, ]


######################################################
#################### SECTION 4 #######################
######################################################
#Section 4 - Classify cells based on euclidean distance
#create merged data frame with both classical and basal scores
basal_score <- as.data.frame(auc_scores_basal)
classical_score <- as.data.frame(auc_scores_classical)
coexpressor_score <- as.data.frame(auc_scores_IC)
merged <- merge(basal_score, classical_score, by = 'row.names', all = TRUE) 
rownames(merged) <- merged$Row.names
#rownames(merged) <- NULL
merged <- merge(merged, coexpressor_score, by = 'row.names', all = TRUE)
rownames(merged) <- merged$Row.names

#calculate euclidean distance
merged$euclidean_dist <- apply(merged[, c("auc_scores_basal", "auc_scores_classical")], 1, function(x) {
  sqrt((x[1] - x[2])^2)})


mal_counts<-as.data.frame(summary(malignant_counts))
## Correlate Euclidean distance with each gene's expression
correlations <- apply(malignant_counts, 1, function(gene_expr) {
  cor(gene_expr, merged$euclidean_dist, method = "pearson")
})


# Classify cells based on Euclidean distance and expression scores
merged$cell_classification <- ifelse(merged$euclidean_dist < 0.01, "Co-Expressor", 
                                             ifelse(merged$auc_scores_basal > merged$auc_scores_classical, "Basal", "Classical"))

merged$cell_stratification <- ifelse(merged$cell_classification == "Basal" & merged$auc_scores_basal >= 0.25, "Strong_Basal",
                                     ifelse(merged$cell_classification == "Basal" & merged$auc_scores_basal < 0.10, "Weak_Basal",
                                     ifelse(merged$cell_classification == "Basal" & merged$auc_scores_basal < 0.25 & merged$auc_scores_basal >= 0.10, "Basal",
                                            ifelse(merged$cell_classification == "Classical" & merged$auc_scores_classical >= 0.25, "Strong_Classical",
                                                   ifelse(merged$cell_classification == "Classical" & merged$auc_scores_classical < 0.25 & merged$auc_scores_classical >= 0.10, "Classical",
                                                    ifelse(merged$cell_classification == "Classical" & merged$auc_scores_classical < 0.10, "Weak_Classical","Co_Expressor"))))))


# head(merged)
subtypetable1 <- table(merged$cell_classification)
subtypetable2 <- table(merged$cell_stratification)

# Add the Euclidean distance and IC scores to the Seurat object's metadata
malignant[['level_2B_clusters']] <- merged$cell_classification[match(rownames(malignant@meta.data), rownames(merged))] 

malignant@meta.data$classical_score <- merged$auc_scores_classical[match(rownames(malignant@meta.data), rownames(merged))]
malignant@meta.data$basal_score <- merged$auc_scores_basal[match(rownames(malignant@meta.data), rownames(merged))]
malignant@meta.data$coexpressor_score <- merged$auc_scores_IC[match(rownames(malignant@meta.data), rownames(merged))]

malignant@meta.data$polarity_score <- malignant@meta.data$classical_score - malignant@meta.data$basal_score

#unique(malignant@meta.data)
malignant <- SetIdent(malignant, value='level_2B_clusters')


#plot correlation of AUCell scores
plot <- FeatureScatter(
  malignant,
  "polarity_score",
  "coexpressor_score",
  smooth = FALSE,
  combine = FALSE,
  slot = "data",
  plot.cor = TRUE,
  cols = c("blue","green","orange")
) + 
  ggtitle("AUCell Polarity Score vs Coexpressor Score") +
  xlab("Polarity Score (Classical Score - Basal Score)") +
  ylab("Coexpressor Score")

ggsave(plot_section_4_1_correlation_AUCell_scores_path)

# Plot cell classical score vs basal score
plot <- FeatureScatter(
  malignant,
  "classical_score",
  "basal_score",
  smooth = FALSE,
  combine = FALSE,
  slot = "data",
  plot.cor = TRUE,
  cols = c("blue","green","orange")
) + 
  ggtitle("AUCell Classical Score vs Basal Score") +
  xlab("Classical Score") +
  ylab("Basal Score")

ggsave(plot_section_4_2_correlation_classical_basal_scores_path)


#create violin plot of marker scores stratified by collapsed cell type
plot <- VlnPlot(malignant, group.by = "level_2B_clusters", features = c("classical_score", "basal_score","coexpressor_score"), ncol = 2, pt.size = 0, cols = c("orange","blue","green"), combine = FALSE) 

my_comparisons <- list( c("Classical", "Basal"), c("Classical", "Co-Expressor"), c("Basal", "Co-Expressor"))

plot[[1]] +
  geom_boxplot(width = 0.2) +
  ggtitle("Classical Score vs Tumor Cell Subtype") +
  xlab("Tumor Cell Subtype") +
  ylab("Classical Score")
ggsave(plot_section_4_3_classical_scores_vs_subtype_path)

plot[[2]] +
  geom_boxplot(width = 0.25) +
  ggtitle("Basal Score vs Tumor Cell Subtype") +
  xlab("Tumor Cell Subtype") +
  ylab("Basal Score")
ggsave(plot_section_4_4_basal_scores_vs_subtype_path)

plot[[3]] +
  geom_boxplot(width = 0.25) +
  ggtitle("Co-Expressor Score vs Tumor Cell Subtype") +
  xlab("Tumor Cell Subtype") +
  ylab("Co-Expressor Score")
ggsave(plot_section_4_5_coexpressor_scores_vs_subtype_path)


######################################################
#################### SECTION 5 #######################
######################################################
#Section 5 - Identify Differentially-Expressed Genes Between Basal and Classical Cells

#SKIPED


######################################################
#################### SECTION 6 #######################
######################################################
#Section 6 - Label Basal and Classical Cells in Main Object
cell_ids <- rownames(seurat_obj_norm@meta.data)
seurat_obj_norm@meta.data$level_2B_clusters <- seurat_obj_norm@meta.data$level_2_clusters
seurat_obj_norm@meta.data$level_3B_clusters <- seurat_obj_norm@meta.data$level_3_clusters %>%
  as.character()
seurat_obj_norm@meta.data$level_3C_clusters <- seurat_obj_norm@meta.data$level_3B_clusters %>%
  as.character()

seurat_obj_norm@meta.data$level_2B_clusters[cell_ids %in% rownames(merged)] <- merged$cell_classification[match(cell_ids[cell_ids %in% rownames(merged)],rownames(merged))]

seurat_obj_norm@meta.data$level_3B_clusters[cell_ids %in% rownames(merged)] <- merged$cell_classification[match(cell_ids[cell_ids %in% rownames(merged)],rownames(merged))]

seurat_obj_norm@meta.data$level_3C_clusters[cell_ids %in% rownames(merged)] <- merged$cell_stratification[match(cell_ids[cell_ids %in% rownames(merged)],rownames(merged))]

malignant_labeled <- subset(seurat_obj_norm, subset = level_1_clusters == "Tumor")

#create boxplots of distribution of basal and classical cells
dittoBarPlot(
    object = malignant_labeled,
    var = "level_2B_clusters",
    group.by = "patient",
    main = "Tumor Subtype Frequency")
ggsave(plot_section_6_1_boxplots_distribution_basal_classical, width=10, height=8)

# Print the updated metadata to verify changes
table(seurat_obj_norm@meta.data$malignant)

#save object with basal and classical labels
saveRDS(seurat_obj_norm, rds_path_full_cohort_semi_sup_insitutype_tumor_subtype_labels)


######################################################
#################### SECTION 8 #######################
######################################################
#Section 8 - Label Fibroblast Cells
fibroblast <- subset(seurat_obj_norm, subset = level_2B_clusters %in% c("FAP_positive_mycaf", "FAP_negative_mycaf","iCAF"))

dittoBarPlot(
    object = fibroblast,
    var = "level_3_clusters",
    group.by = "patient",
    main = "Fibroblast Subtype Frequency")
ggsave(plot_section_8_1_boxplot_level_3_cluster_fibrablast_subtype_frequency, width=10, height=8)

dittoBarPlot(
    object = fibroblast,
    var = "level_3B_clusters",
    group.by = "patient",
    main = "Fibroblast Subtype Frequency")
ggsave(plot_section_8_2_boxplot_level_3B_cluster_fibrablast_subtype_frequency, width=10, height=8)


Idents(seurat_obj_norm) <- "level_2_clusters"

caf_cells <- WhichCells(seurat_obj_norm, idents = c("myCAF","iCAF","Mesenchymal.Cell"))

FAP_positive_caf <- WhichCells(seurat_obj_norm, cells = caf_cells, expression = FAP > 0)
FAP_negative_caf <- WhichCells(seurat_obj_norm, cells = caf_cells, expression = FAP == 0)

# Assign 'FAP+ve' and 'FAP-ve' statuses to mycaf cells
seurat_obj_norm@meta.data$level_2B_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_positive_caf] <- "FAP_positive_caf"
seurat_obj_norm@meta.data$level_2B_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_negative_caf] <- "FAP_negative_caf"

# Assign 'FAP+ve' and 'FAP-ve' statuses to mycaf cells
seurat_obj_norm@meta.data$level_3B_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_positive_caf] <- "FAP_positive_caf"
seurat_obj_norm@meta.data$level_3B_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_negative_caf] <- "FAP_negative_caf"

seurat_obj_norm@meta.data$level_3C_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_positive_caf] <- "FAP_positive_caf"
seurat_obj_norm@meta.data$level_3C_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_negative_caf] <- "FAP_negative_caf"

caf <- subset(seurat_obj_norm, cells = caf_cells)

dittoBarPlot(
    object = caf,
    var = "level_2_clusters",
    group.by = "patient",
    main = "Fibroblast Subtype Frequency")
ggsave(plot_section_8_3_boxplot_level_2_cluster_fibrablast_subtype_frequency_cafs, width=10, height=8)


dittoBarPlot(
    object = caf,
    var = "level_3B_clusters",
    group.by = "patient",
    main = "Fibroblast Subtype Frequency")
ggsave(plot_section_8_4_boxplot_level_3B_cluster_fibrablast_subtype_frequency_cafs, width=10, height=8)


VlnPlot(caf, features = "FAP")
ggsave(plot_section_8_5_fap_caf_vlnplot)


# Check the distribution of the feature
feature_data <- FetchData(caf, vars = "FAP")
summary(feature_data)

# Extract the 'FAP' column as a numeric vector
fap_values <- feature_data$FAP

# Ensure that 'fap_values' is numeric
if (!is.numeric(fap_values)) {
  stop("The 'FAP' variable is not numeric.")
}

# Create histogram
hist(fap_values, breaks = 6, main = "Histogram of FAP expression", xlab = "Expression Level")

#hist(feature_data, breaks = 6, main = "Histogram of FAP expression", xlab = "Expression Level")

fibroblasttable1 <- table(caf@meta.data$level_3B_clusters)


#split cells into fap+ mycaf and fap+icaf
Idents(seurat_obj_norm) <- "level_3_clusters"

mycaf_cells <- WhichCells(seurat_obj_norm, idents = "myCAF")

FAP_positive_mycaf <- WhichCells(seurat_obj_norm, cells = mycaf_cells, expression = FAP > 0)
FAP_negative_mycaf <- WhichCells(seurat_obj_norm, cells = mycaf_cells, expression = FAP == 0)

# Assign 'FAP+ve' and 'FAP-ve' statuses to mycaf cells
seurat_obj_norm@meta.data$level_3B_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_positive_mycaf] <- "FAP_positive_mycaf"
seurat_obj_norm@meta.data$level_3B_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_negative_mycaf] <- "FAP_negative_mycaf"

seurat_obj_norm@meta.data$level_3C_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_positive_mycaf] <- "FAP_positive_mycaf"
seurat_obj_norm@meta.data$level_3C_clusters[rownames(seurat_obj_norm@meta.data) %in% FAP_negative_mycaf] <- "FAP_negative_mycaf"

#save seurat object

#saveRDS(seurat_obj_norm, here("./data/interim/cosmx/human/seurat_objects/BTC_full_cohort_semi_sup_insitutype_fully_labeled_07152024.RDS"))
saveRDS(seurat_obj_norm, file=file.path(output_folder_path, "BTC_full_cohort_semi_sup_insitutype_fully_labeled.RDS"))

# Function to run AUCell to identify CAFs close to FAP+ CAFs and generate plots
identify_fap_positive_cafs <- function(caf, genesets, output_dir, sample_name, full_cohort_norm) {
    # Ensure output directory exists for plots
    plots_dir <- file.path(output_dir, "plots")
    if (!dir.exists(plots_dir)) {
        dir.create(plots_dir, recursive = TRUE)
    }

    # Get counts data from the 'caf' Seurat object
    fibroblast_counts <- GetAssayData(object = caf, slot = "counts")

    # Build AUCell rankings
    fibroblast_rankings <- AUCell_buildRankings(fibroblast_counts, plotStats = TRUE)

    # Identify fibroblast gene sets for FAP+ CAFs
    fibroblast_genesets <- genesets[genesets$gs_name == "FAP_pos_caf", ]

    # Filter gene sets to include genes present in the rankings
    fibroblast_genesets_cosmx <- fibroblast_genesets[fibroblast_genesets$gene %in% rownames(fibroblast_rankings), ]

    # Create list of gene sets
    fibroblast_genesets_cosmx <- split(fibroblast_genesets_cosmx$gene, fibroblast_genesets_cosmx$gs_name)

    # Run AUCell to calculate AUC scores
    cells_AUC <- AUCell_calcAUC(fibroblast_genesets_cosmx, fibroblast_rankings)

    # Explore AUC thresholds (plots histogram)
    # Save the histogram directly to a file
    histogram_filename <- file.path(plots_dir, paste0(sample_name, "_FAP_CAF_AUC_Histogram.png"))
    png(filename = histogram_filename, width = 800, height = 600)
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, assign = TRUE)
    dev.off()

    # Extract the gene set name
    geneSetName <- rownames(cells_AUC)[grep("FAP", rownames(cells_AUC))]

    # Get AUC scores for FAP
    auc_scores_FAP <- getAUC(cells_AUC)[geneSetName, ]

    # Merge AUC scores into a data frame
    merged <- data.frame(auc_scores_FAP)

    # Add AUC scores to 'caf' metadata
    caf@meta.data$FAP_CAF_score <- merged$auc_scores_FAP[match(rownames(caf@meta.data), rownames(merged))]

    # Classify cells based on quantiles of AUC scores
    caf@meta.data$FAP_CAF_score_level <- ifelse(
        caf@meta.data$FAP_CAF_score > quantile(caf@meta.data$FAP_CAF_score, probs = 0.75),
        "FAP_like",
        ifelse(
            caf@meta.data$FAP_CAF_score < quantile(caf@meta.data$FAP_CAF_score, probs = 0.25),
            "FAP_low",
            "FAP_intermediate"
        )
    )

    # Update 'full_cohort_norm' metadata with new classifications
    cell_ids <- rownames(full_cohort_norm@meta.data)
    matching_cells <- cell_ids %in% rownames(caf@meta.data)
    matched_indices <- match(cell_ids[matching_cells], rownames(caf@meta.data))

    full_cohort_norm@meta.data$level_2B_clusters[matching_cells] <- caf@meta.data$FAP_CAF_score_level[matched_indices]
    full_cohort_norm@meta.data$level_3B_clusters[matching_cells] <- caf@meta.data$FAP_CAF_score_level[matched_indices]
    full_cohort_norm@meta.data$level_3C_clusters[matching_cells] <- caf@meta.data$FAP_CAF_score_level[matched_indices]

    # Generate dittoBarPlot for Fibroblast Subtype Frequency
    fibroblast_plot <- dittoBarPlot(
        object = caf,
        var = "FAP_CAF_score_level",
        group.by = "patient",
        main = "Fibroblast Subtype Frequency"
    )

    # Save the plot
    ggsave(
        filename = file.path(plots_dir, paste0(sample_name, "_Fibroblast_Subtype_Frequency.png")),
        plot = fibroblast_plot,
        width = 8,
        height = 6
    )

    # Return the updated 'caf' and 'full_cohort_norm' objects
    return(list(caf = caf, full_cohort_norm = full_cohort_norm))
}

#LEN: debugging
#options(error = function() traceback(3))

#CALL identify_fap_positive_cafs
result <- identify_fap_positive_cafs(caf, genesets, output_folder_path, sample_name, seurat_obj_norm)
caf <- result$caf
full_cohort_norm <- result$full_cohort_norm

# Save the updated Seurat objects
saveRDS(caf, file = file.path(output_folder_path, paste0(sample_name, "_caf_updated.RDS")))
saveRDS(full_cohort_norm, file = file.path(output_folder_path, paste0(sample_name, "_full_cohort_norm_updated.RDS")))


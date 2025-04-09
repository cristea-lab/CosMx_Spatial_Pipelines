library("Seurat")
library("sceasy")
library("reticulate")
library("dplyr")
library("SpatialDecon")
library("InSituType")
library("tidyverse")
library("here")
library("dittoSeq")

###############################################################
#################### COMMAND LINE INPUT #######################
###############################################################
## Example command-line input: 
# Rscript insitutype_single_sample_processing_pipeline.R PCA429 /cristealab/ajordan/projects/BTC_analysis/data/processed/cosmx/human/resegmented_data/PCA429 /cristealab/ajordan/projects/BTC_analysis/data/external/insitutype_references/Pancreas.profiles.csv /cristealab/xiwang/Outputs/CosMx_Spatial_Pipelines/output_PCA429_08262024/InSituType 12 14


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) < 1) {
  stop("Please provide one or more file paths as command-line arguments.")
}
if (length(args) != 7) {
    stop("Please double check your input to make sure only 7 command-line inputs are given.")
}


##################################################
#################### INPUT #######################
##################################################
sample_name = args[1]
print(paste("Input sample name: ", sample_name))

sample_rds_path = args[2]
print(paste("Input sample rds path: ", sample_rds_path))


sample_folder_path = args[3]
print(paste("Input sample folder path from AtoMx: ", sample_folder_path))

output_folder_path = args[4]
print(paste("Input output folder: ", output_folder_path))

insitutype_reference_sig_csv_file_path = args[5]
print(paste("Input inSituType ref sig profile from Nanostring: ", insitutype_reference_sig_csv_file_path))


# num_cluster, min and max represents the lower and upper end of the range (both inclusive)
min_num_cluster = args[6]
print(paste("Input lower end of the num of clusters: ", min_num_cluster))

max_num_cluster = args[7]
print(paste("Input higher end of the num of clusters: ", max_num_cluster))

#LEN: 2025-01-14 These params moved to 1_QC_preprocessing_single_sample_pipeline.R
# User input QC metrice after viewing the plot
#min_nCount_RNA = args[8] # >=
#print(paste("Input QC metrice of min nCount RNA: ", min_nCount_RNA))

#max_nFeature_negprobes = args[9] # <
#print(paste("Input QC metrice of max nFeature negprobes: ", max_nFeature_negprobes))


##################################################
#################### OUTPUT ######################
##################################################

#output_Rds_file_path = paste0(output_folder_path, "/A_1_pre_QC_sample.rds")

#LEN: 2025-02-10 DEPRECATED
#post_QC_and_normalization_20_counts_per_cell_rds_file_path = paste0(output_folder_path, "/A_2_post_QC_and_normalization_20_counts_per_cell.rds")



#-----------------------------------------------------------------------------------------------------

##################################################
#################### START #######################
##################################################



#################### LOAD SAMPLE ####################
#seurat_obj = readRDS(output_Rds_file_path)
seurat_obj_norm = readRDS(sample_rds_path)

#######################################################
#################### PREPROCESSING ####################
#######################################################
meta_data_file_path = paste0(sample_folder_path, "/", sample_name, "_metadata_file.csv")

# InSituType reference profile from Nanostring 
pancreas_nanostring_ref_file_path = insitutype_reference_sig_csv_file_path


#LEN: 2025-01-14 Normalization now handled by UMAP script which is upstream from this
#normalize data
#seurat_obj_norm <- NormalizeData(seurat_obj)
#seurat_obj_norm <- JoinLayers(seurat_obj_norm)

#saveRDS(seurat_obj_norm, post_QC_and_normalization_20_counts_per_cell_rds_file_path)


####################################################
#################### InSituType ####################
####################################################

# Read reference file
pancreas_nanostring = read.csv(pancreas_nanostring_ref_file_path, row.names = 1) %>% as.matrix()

# Prepare cosmx data for insitutype run
# btain and store count data
counts <- LayerData(seurat_obj_norm, layer = "counts") %>% as.matrix() %>% t()

neg_count_vector <- seurat_obj_norm@meta.data$nCount_negprobes %>% as.data.frame() %>% rowMeans()

# Prepare cell anchors for insitutype run
# calculate anchor stats
anchor_stats <- get_anchor_stats(counts = counts,
                           neg = neg_count_vector,
                           profiles = pancreas_nanostring)

# estimate per-cell background as a fraction of total counts:
negmean.per.totcount <- mean(neg_count_vector) / mean(rowSums(counts))
per.cell.bg <- rowSums(counts) * negmean.per.totcount

# choose cell anchors
anchors <- choose_anchors_from_stats(counts = counts, 
                                     neg = neg_count_vector, 
                                     bg = per.cell.bg,
                                     anchorstats = anchor_stats, 
                                     n_cells = 5000, 
                                     min_cosine = 0.4, 
                                     min_scaled_llr = 0.02, 
                                     insufficient_anchors_thresh = 20)

                                     #use anchors to update reference profiles
updatedprofiles <- updateReferenceProfiles(reference_profiles = pancreas_nanostring, 
                                           counts = counts, 
                                           neg = neg_count_vector, 
                                           bg = per.cell.bg,
                                           anchors = anchors) 


#create cell cohorts based on mIF data
#mif_data <- seurat_obj_norm@meta.data[,c("Mean.PanCK","Mean.CD45","Mean.CD68_CK8_18")] %>% as.matrix()
mif_data <- seurat_obj_norm@meta.data[,c("Mean.PanCK","Mean.CD45")] %>% as.matrix()
mif_cohort <- fastCohorting(mif_data, gaussian_transform = TRUE)


#run semi-supervised insitutype run
#set seed for reproducibility
#set random seed
set.seed(123)

# Given the range of cluster, run through all the number of clusters
for (num_cluster in min_num_cluster:max_num_cluster) {
  print(paste0("Running inSituType clustering analysis: ", num_cluster, " number of clusters"))

  # Define the output format
  cluster_output_folder_path = paste0(output_folder_path, "/", num_cluster,"_clusters")
  # Make the output path
  dir.create(cluster_output_folder_path)


  flightpath_initial_cluster_file_path = paste0(cluster_output_folder_path, "/2_1_flightpath_initial_cluster.png")
  flightpath_new_cluster_file_path = paste0(cluster_output_folder_path, "/2_2_flightpath_new_cluster.png")

  semi_sup_insitutype_refined_rds_file_path = paste0(cluster_output_folder_path, "/B_1_semi_sup_insitutype_refined.rds")

  dotplot_marker_genes_file_path =  paste0(cluster_output_folder_path, "/2_3_dotplot_marker_genes.png")

  semi_sup_insitutype_fully_labeled_rds_file_path = paste0(cluster_output_folder_path, "/B_2_semi_sup_insitutype_fully_labeled.rds")

  barplot_level_1_clusters_file_path = paste0(cluster_output_folder_path, "/2_4_barplot_level_1_clusters.png")
  barplot_level_2_clusters_file_path = paste0(cluster_output_folder_path, "/2_5_barplot_level_2_clusters.png")
  barplot_level_3_clusters_file_path = paste0(cluster_output_folder_path, "/2_6_barplot_level_3_clusters.png")

  violin_plot_marker_scores_1_file_path = paste0(cluster_output_folder_path, "/3_1_violin_plot_of_marker_scores_stratified_by_collapsed_cell_type_1.png")
  violin_plot_marker_scores_2_file_path = paste0(cluster_output_folder_path, "/3_2_violin_plot_of_marker_scores_stratified_by_collapsed_cell_type_2.png")

  spatial_plot_file_path = paste0(cluster_output_folder_path, "/4_1_spatial_plot.png")
  spatial_plot_file_path_level_1_clusters = paste0(cluster_output_folder_path, "/4_2_spatial_plot_level_1_clusters.png")


  #run semi-supervised insitutype clustering analysis
  initial_clusters <- insitutype(x = counts, neg = neg_count_vector, n_clusts = seq(1, num_cluster, by = 1),
    reference_profiles = updatedprofiles$updated_profiles,
    update_reference_profiles = FALSE
  ) 

  #generate flightpath plot
  colors <- colorCellTypes(freqs = table(initial_clusters$clust), palette = "brewers")
  fp <- flightpath_plot(insitutype_result = initial_clusters, col = colors[initial_clusters$clust])
  ggsave(flightpath_initial_cluster_file_path, width=10, height=8)


  #identify marker genes after insitutype run
  #set cell identities as insitutype clusters
  seurat_obj_norm@meta.data$insitutype_clusters <- initial_clusters$clust
  Idents(seurat_obj_norm) <- "insitutype_clusters"

  #identify marker genes for insitutype clusters
  cosmx_markers <- FindAllMarkers(seurat_obj_norm, only.pos = TRUE)

  new_clusters <- refineClusters(
    logliks = initial_clusters$logliks,
    # merges = c("Fibroblast" = "iCAF", "i" = "myCAF"),
    # subclustering via refineClusters is not recommended for semi-supervised
    # results
    #to_delete = c("T.cell.CD4"),
    subcluster = NULL,
    counts = counts,
    neg = neg_count_vector
  ) 

  fp <- flightpath_plot(insitutype_result = new_clusters, col = colors[new_clusters$clust])
  ggsave(flightpath_new_cluster_file_path, width=10, height=8)

  seurat_obj_norm@meta.data$insitutype_clusters <- new_clusters$clust
  Idents(seurat_obj_norm) <- "insitutype_clusters"

  #identify marker genes for insitutype clusters
  cosmx_markers <- FindAllMarkers(seurat_obj_norm, only.pos = TRUE)


  #save object for further analysis
  saveRDS(seurat_obj_norm, semi_sup_insitutype_refined_rds_file_path)



  #rename cluster ids
  level_3_clusters <- c("CD8.T.Cell", "CEACAM6+PDAC", "Macrophage.1", "B.Cell", "NK", "myCAF", "Endothelial.Cell","iCAF", "AGR2+PDAC", "Mesenchymal.Cell", "DST+PDAC", "Macrophage.2","LYZ+PDAC","CXCL5+PDAC","Monocyte","Alpha.Cell","Islet.Cell","Plasma.Cell","Acinar.Cell","VEGFA+PDAC","Treg","Delta.Cell","Dendritic.Cell","Beta.Cell")

  names(level_3_clusters) <- levels(seurat_obj_norm)
  seurat_obj_norm <- RenameIdents(seurat_obj_norm, level_3_clusters)
  seurat_obj_norm@meta.data$level_3_clusters <- Idents(seurat_obj_norm)


  #create dotplot of marker genes
  feature_vector <-c("KRT6A/B/C","KRT7","AGR2","CEACAM6","SERPINA1","SERPINA3","CUZD1","CPB1","PRSS2","REG1A","FN1","BGN","COL1A1","ACTA2","FAP","LMNA","COL4A1","PECAM1","CD74","C1QC","CD68","PTPRC","CD3E","CD8A","IGHD","IGKC","IFITM3","IGHG1","IGHG2","IGHM","CSF1R","S100A8","INS","TTR","PGR","ACTG2","MYH11","MYL9")
  DotPlot(seurat_obj_norm, features = feature_vector) + 
      theme(axis.text.x = element_text(size = 6)) +
      RotatedAxis() +
      coord_flip()
  ggsave(dotplot_marker_genes_file_path, width=10, height=8)



  #collapse clusters into tumor, stroma, and immune for level 1 clusters
  seurat_obj_norm@meta.data <- seurat_obj_norm@meta.data %>%
  mutate(level_1_clusters = case_when(
    level_3_clusters %in% c("Monocyte", "Macrophage.1", "Macrophage.2", "B.Cell", "Plasma.Cell", "NK", "Dendritic.Cell", "CD8.T.Cell", "Treg") ~ "Immune",
    level_3_clusters %in% c("CEACAM6+PDAC", "AGR2+PDAC", "DST+PDAC","LYZ+PDAC","CXCL5+PDAC", "VEGFA+PDAC") ~ "Tumor",
    level_3_clusters %in% c("myCAF", "iCAF", "Mesenchymal.Cell", "Endothelial.Cell") ~ "Stroma",
    level_3_clusters %in% c("Islet.Cell","Alpha.Cell", "Beta.Cell", "Delta.Cell", "Acinar.Cell") ~ "Parenchyma",
    TRUE ~ as.character(level_3_clusters)
  ))

  seurat_obj_norm@meta.data <- seurat_obj_norm@meta.data %>%
  mutate(level_2_clusters = case_when(
    level_3_clusters %in% c("CD8.T.Cell", "Treg") ~ "T.Cell",
    level_3_clusters %in% c("B.Cell","Plasma.Cell") ~ "B/Plasma",
    level_3_clusters %in% c("Monocyte", "Macrophage.1", "Macrophage.2", "Dendritic.Cell") ~ "Myeloid",
    level_3_clusters %in% c("CEACAM6+PDAC", "AGR2+PDAC", "DST+PDAC","LYZ+PDAC","CXCL5+PDAC", "VEGFA+PDAC") ~ "Tumor",
    level_3_clusters %in% c("Islet.Cell","Alpha.Cell", "Beta.Cell", "Delta.Cell") ~ "Islets",
    TRUE ~ as.character(level_3_clusters)
  ))

  #set identities to collapsed clusters for downstream plotting
  Idents(seurat_obj_norm) <- "collapsed_clusters"

  #save object for further analysis
  saveRDS(seurat_obj_norm, semi_sup_insitutype_fully_labeled_rds_file_path)

  dittoBarPlot(
      object = seurat_obj_norm,
      var = "level_1_clusters",
      group.by = "patient")
  ggsave(barplot_level_1_clusters_file_path)

  dittoBarPlot(
      object = seurat_obj_norm,
      var = "level_2_clusters",
      group.by = "patient")
  ggsave(barplot_level_2_clusters_file_path)

  dittoBarPlot(
      object = seurat_obj_norm,
      var = "level_3_clusters",
      group.by = "patient")
  ggsave(barplot_level_3_clusters_file_path)


  #####################################
  ######### EXCEL SHEET TABLE #########
  #####################################





  #review mIF staining for each cluster
  #create violin plot of marker scores stratified by collapsed cell type
  #VlnPlot(seurat_obj_norm, group.by = "level_2_clusters", features = c("Mean.PanCK", "Mean.CD45", "Mean.CD68_CK8_18"), ncol = 2, pt.size = 0, y.max = 10000)
  VlnPlot(seurat_obj_norm, group.by = "level_2_clusters", features = c("Mean.PanCK", "Mean.CD45"), ncol = 2, pt.size = 0, y.max = 10000)
  ggsave(violin_plot_marker_scores_1_file_path)

  #create violin plot of marker scores stratified by collapsed cell type
  VlnPlot(seurat_obj_norm, group.by = "level_2_clusters", features = c("CD8A", "CD4", "PTPRC"), ncol = 2, pt.size = 0)
  ggsave(violin_plot_marker_scores_2_file_path)



  #plot clusters in a spatially-resolved manner
  ImageDimPlot(seurat_obj_norm, fov = "fov", axes = TRUE, nmols = 10000)
  ggsave(spatial_plot_file_path)


  #color cases by patient to ensure the tumors are labeled correctly
  cluster_colors <- c(
    "Tumor" = "#0AF5E5",
    "Stroma" = "#8CFF00",
    "Parenchymal" = "#CC66FF",
    "Immune" = "#FF7777"
  )

  Idents(seurat_obj_norm) <- "level_1_clusters"
  #plot clusters in a spatially-resolved manner
  ImageDimPlot(seurat_obj_norm, fov = "fov", split.by = "patient", axes = TRUE, nmols = 10000)
  ggsave(spatial_plot_file_path_level_1_clusters)
}

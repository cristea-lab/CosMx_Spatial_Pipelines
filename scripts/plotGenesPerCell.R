# Len Taing 2025 (TGBTG)
# Given a path to RAW cosmx data OR a seurat object representing cosmx data
#
# PLOTS the distribution of Genes Per Cell
# and IF available, colorizes the tissue image with genes per cell

library("Seurat")
library("ggplot2")

loadCosMx <- function(data_path) {
  seurat_obj <- LoadNanostring(data.dir=data_path, fov="fov", assay="RNA")
  
  #load metadata
  sample_name <- basename(data_path)
  #NOTE: LoadNanostring does NOT read in all of the metadata so we have to do this
  meta_data_file_path = paste0(data_path, "/", sample_name, "_metadata_file.csv")
  #LEN: 2025-05-21 ADD check in case the metafile is .csv.gz
  if (!file.exists(meta_data_file_path)) {
    meta_data_file_path = paste0(meta_data_file_path, ".gz")
  }

  metadata <- read.csv(meta_data_file_path, header = TRUE)
  original_metadata <- seurat_obj@meta.data
  metadata <- cbind(original_metadata, metadata)
  seurat_obj@meta.data <- metadata
  
  #For some reason, it thinks that nCount_RNA and nFeature_RNA column names are
  #HYP: it's because we load the metadata twice
  #duplicated--deduplicate
  colnames(seurat_obj@meta.data) <- gsub("^\\s+|\\s+$", "", colnames(seurat_obj@meta.data))  # Trim leading/trailing spaces
  colnames(seurat_obj@meta.data) <- make.names(colnames(seurat_obj@meta.data), unique = TRUE)  # Ensure uniqueness

  return(seurat_obj)
}

#Plot the distrubtion of genes per cell
# input: seurat object, a output png path, and an optional int on where to
# draw a vertical line in the plot, e.g. the QC threshold
plotGenesPerCellHistogram <- function(seurat_obj, distrib_plot_out, x_vert_line=NULL){
  # 1) Compute mean and median
  mean_val   <- mean(seurat_obj@meta.data$nFeature_RNA, na.rm = TRUE)
  median_val <- median(seurat_obj@meta.data$nFeature_RNA, na.rm = TRUE)
  
  # 2) Create a data frame for vertical lines:
  #    We embed the numeric values directly in the factor labels so
  #    they show up in the legend.
  stats_lines <- data.frame(
    # If the user line is NULL, we omit it
    label = c(
      sprintf("Mean=%.2f", mean_val),
      sprintf("Median=%.2f", median_val),
      if(!is.null(x_vert_line)) sprintf("User=%.2f", x_vert_line) else character(0)
    ),
    x_val = c(
      mean_val,
      median_val,
      if(!is.null(x_vert_line)) x_vert_line else numeric(0)
    )
  )
  
  # 3) Build the base plot: disable legend for the histogram
  p <- ggplot(seurat_obj@meta.data, aes(x = nFeature_RNA)) +
    geom_histogram(show.legend = FALSE) +
    theme_classic() +
    xlab("No. of Genes Detected") +
    ylab("No. of Cells")
  
  # 4) Add vertical lines for mean, median, and (optionally) user line.
  #    We map color = label so each line is a different color and appears in the legend.
  p <- p +
    geom_vline(
      data = stats_lines,
      aes(xintercept = x_val, color = label),
      linetype = "dashed",
      size = 1
    ) +
    # Remove the legend title (you can keep or rename if you like)
    scale_color_discrete(name = NULL)
  
  # 5) Save the resulting plot
  ggsave(distrib_plot_out, plot = p, width = 6, height = 4)
}

plotGenesPerCellTissue <- function(seurat_obj, tissue_plot_out) {
  my_breaks <- c(seq(0, 300, by=30), Inf)
  bin_labels <- c(
    paste(seq(0, 270, by=30), seq(30, 300, by=30), sep="-"),
    ">300"
  )

  my_palette <- colorRampPalette(c("gray85", "lightblue", "blue", "red"))(length(bin_labels))
  seurat_obj$nFeature_RNA_bin <- cut(
    seurat_obj$nFeature_RNA,
    breaks = my_breaks,
    labels = bin_labels,
    include.lowest = TRUE,
    right = FALSE
  )
  
  ImageDimPlot(
      object = seurat_obj,
      fov = "fov",            # adjust to your dataset
      group.by = "nFeature_RNA_bin",
      cols= my_palette
    )
  ggsave(
    filename = tissue_plot_out,
    width = 6,    # adjust width (in inches) as needed
    height = 5,   # adjust height (in inches) as needed
    dpi = 300     # resolution in dots-per-inch
    #plot = p
  )
}

main <- function(input_data_path, distrib_plot_out, tissue_plot_out) {
  #LOAD the data
  seurat_obj <- loadCosMx(input_data_path)
  plotGenesPerCellHistogram(seurat_obj, distrib_plot_out, 20)
  #should check for valid FOV
  plotGenesPerCellTissue(seurat_obj, tissue_plot_out)
}

  
args <- commandArgs( trailingOnly = TRUE )
arg_cosmx_input_path = args[1] #or seurat object path
arg_distrib_plot_out = args[2]
arg_tissue_plot_out = args[3]

main(arg_cosmx_input_path, arg_distrib_plot_out, arg_tissue_plot_out)

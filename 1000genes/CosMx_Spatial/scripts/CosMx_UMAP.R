library(Seurat)
library(ggplot2)
#script to run UMAP analysis after initial filtering (should come right after normalization)
#input - QC'ed seurat object (after low-count cells, negative probes, and background probes have been removed)

#Runs normalization - umap cmds returns normalized/scaled/umapped object
runUMAP <- function(seurat_obj_QC_filtered) {
   seurat_obj_norm <- NormalizeData(seurat_obj_QC_filtered)
   seurat_obj_norm <- JoinLayers(seurat_obj_norm)
   seurat_obj_norm <- ScaleData(seurat_obj_norm, features = Features(seurat_obj_norm))
   seurat_obj_norm <- RunPCA(seurat_obj_norm, features = Features(seurat_obj_norm), seed.use = 1)
   seurat_obj_norm <- RunUMAP(object = seurat_obj_norm, dims = 1:50, spread = 2, min.dist = 0.05)
   return(seurat_obj_norm)
}

# plotUMAP
#Given a seurat object that has been normalized and whose dimensions have been reduced
#with UMAP (see runUMAP), a vector of genes/mif markers to plot in the umap, and plot
#output path, outputs 1) the umap plot, 2) plot of the genes/mif features in the umap,
#and 3) plot of cluster/Gene sets to the plot output path
#
# the umap will simple be called umap.png
# the features in umap pngs will be called {feature}_umap.png
# the cluster/feature pngs will be called {cluster}_{feature}_violin.png

#TODO: customize cluster/features vln plots
plotUMAP <- function(seurat_obj_norm, features, umap_plot_dir) {

   #Idents(seurat_obj_norm) <- "level_4_clusters" #identities should be set to the highest level of cluster annotation

   #LEN: BUG: we're not seeing clusters in the UMAP png
   #CAUSE: clustering happens downstream of this file!
   #HACK: for now we'll do seurat clustering here
   #NOTE: these clusters will not persist! b/c plotUMAP is called AFTER the
   #seurat obj is saved to file!
   #NOTE: resolution and dims are fixed!!!
   seurat_obj_norm <- FindNeighbors(seurat_obj_norm, reduction="pca", dims=1:50)
   seurat_obj_norm <- FindClusters(seurat_obj_norm, resolution=0.3)

   DimPlot(seurat_obj_norm, reduction = "umap") #first, a plot of the umap with the cluster annotations overlaid upon it should be generated
   png_path <- paste0(umap_plot_dir, "/umap.png")
   ggsave(png_path, width=10, height=8)

   #Next, the user should be able to create UMAPs with various genes of interest depicted or the mIF data. Len - would it be possible for
   # the user to type in what metadata they wish to display? This will change based on the clusters of interest (KRT7 for a basal tumor, etc. etc.
   # These feature plots below are examples and these represent the output
   #FeaturePlot(seurat_obj_norm, features = "INS", reduction = "umap")
   #FeaturePlot(seurat_obj_norm, features = "Max.PanCK", reduction = "umap")
   #FeaturePlot(seurat_obj_norm, features = "Max.CD45", reduction = "umap")
   #FeaturePlot(seurat_obj_norm, features = "Max.CD68_CK8_18", reduction = "umap")
   for (feature in features) {
     png_path <- paste0(umap_plot_dir, "/", feature, "_umap.png")
     FeaturePlot(seurat_obj_norm, features = feature, reduction = "umap")
     ggsave(png_path, width=10, height=8)
   }


   #Next, the user can create violin plots of the desired cluster annotations and then plot one or more features of interest (either Mean.PanCK, Mean.CD45, or Mean.CD68_CK8_18)
   # Len - could the user somehow select the cluster annotations they wish to plot via a drop-down menue?
   # png_path <- paste0(umap_plot_dir, "/", "level_1_Mean.PanCK_Mean.CD45_volin.png")
   # VlnPlot(seurat_obj_norm, group.by = "level_1_clusters", features = c("Mean.PanCK", "Mean.CD45"), ncol = 2, pt.size = 0, log = TRUE)
   # ggsave(png_path, width=10, height=8)
   
   # png_path <- paste0(umap_plot_dir, "/", "level_2_Mean.PanCK_Mean.CD45_volin.png")
   # VlnPlot(seurat_obj_norm, group.by = "level_2_clusters", features = c("Mean.PanCK", "Mean.CD45"), ncol = 2, pt.size = 0, log = TRUE)
   # ggsave(png_path, width=10, height=8)

   # png_path <- paste0(umap_plot_dir, "/", "level_3_Mean.PanCK_Mean.CD45_volin.png")
   # VlnPlot(seurat_obj_norm, group.by = "level_3_clusters", features = c("Mean.PanCK", "Mean.CD45"), ncol = 2, pt.size = 0, log = TRUE)
   # ggsave(png_path, width=10, height=8)
   
   # png_path <- paste0(umap_plot_dir, "/", "level_4_Mean.PanCK_Mean.CD45_volin.png")
   # VlnPlot(seurat_obj_norm, group.by = "level_4_clusters", features = c("Mean.PanCK", "Mean.CD45"), ncol = 2, pt.size = 0, log = TRUE)
   # ggsave(png_path, width=10, height=8)
#
}

#output - seurat object with scaled data, and PCA/UMAP dimensions

args <- commandArgs( trailingOnly = TRUE )
#input - QC'ed seurat object (after low-count cells, negative probes, and background probes have been removed) - .rds file
arg_cosmx_input_path = args[1]
args_features = strsplit(args[2], ",")[[1]] #comma separated list of features to plot in umap
#FOR testing
#args_features = c("INS", "Max.PanCK", "Max.CD45", "Max.CD68_CK8_18")
arg_seurat_normScaleUMAP_out = args[3]
args_seurat_plot_out_dir = args[4]

#load seurat obj:
seurat_obj <- readRDS(arg_cosmx_input_path)
seurat_obj_norm <- runUMAP(seurat_obj)
#SAVE seurat_obj_norm
saveRDS(seurat_obj_norm, arg_seurat_normScaleUMAP_out)

#plots
plotUMAP(seurat_obj_norm, args_features, args_seurat_plot_out_dir)

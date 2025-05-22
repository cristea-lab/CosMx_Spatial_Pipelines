library("Seurat")
library("sceasy")
library("reticulate")
library("dplyr")
library("SpatialDecon")
library("InSituType")
library("tidyverse")
library("here")
library("dittoSeq")


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) < 1) {
  stop("Please provide one or more file paths as command-line arguments.")
}

if (length(args) != 5) {
    stop("Please double check your input to make sure only 5 command-line inputs are given.")
}

##################################################
#################### INPUT #######################
##################################################

sample_name = args[1]
print(paste("Input sample name: ", sample_name))


sample_folder_path = args[2]
print(paste("Input sample folder path from AtoMx: ", sample_folder_path))

min_nCount_RNA = args[3]
max_proportion_negprobes = args[4]

output_folder_path = args[5]
print(paste("Input output folder: ", output_folder_path))

##################################################
#################### OUTPUT ######################
##################################################
meta_data_file_path = paste0(sample_folder_path, "/", sample_name, "_metadata_file.csv")

#LEN: 2025-05-21 ADD check in case the metafile is .csv.gz
if (!file.exists(meta_data_file_path)) {
  meta_data_file_path = paste0(meta_data_file_path, ".gz")
}

output_Rds_file_path = paste0(output_folder_path, "/A_1_pre_QC_and_filtered_sample.rds")

# QC
probe_count_histogram_file_path = paste0(output_folder_path, "/1_1_QC_probe_count_histogram.png")
feature_count_histogram_file_path = paste0(output_folder_path, "/1_2_QC_feature_count_histogram.png")
negative_probe_proportion_histogram_file_path = paste0(output_folder_path, "/1_3_QC_negative_probe_proportion_histogram.png")
cell_area_histogram_file_path = paste0(output_folder_path, "/1_4_QC_cell_area_histogram.png")

#-----------------------------------------------------------------------------------------------------

##################################################
#################### START #######################
##################################################



#################### LOAD SAMPLE ####################
#load sample
seurat_obj <- LoadNanostring(
  data.dir = sample_folder_path,
  fov = "fov",
  assay = "Nanostring"
)

#read in metadata to add to seurat object
metadata <- read.csv(meta_data_file_path, header = TRUE)

#withdraw original count and gene metadata from the seurat object
original_metadata <- seurat_obj@meta.data
  
#create fully unified metadata prior to adding to seurat object
metadata <- cbind(original_metadata, metadata)

#add metadata to seurat object
seurat_obj@meta.data <- metadata


#####################################################
#################### SAVE SEURAT ####################
#####################################################
#create function to edit seurat objects
edit_fov_column <- function(seurat_obj) {
    seurat_obj@meta.data$global_fov <- paste(seurat_obj@meta.data$Run_Tissue_name, seurat_obj@meta.data$fov, sep = "_")
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
        mutate(disease = "PDAC")
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
        relocate(disease, .after=orig.ident)
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
        rename(FOV = fov)

    return(seurat_obj)
}

seurat_obj = edit_fov_column(seurat_obj)

# Add column for the corresponding tissue name
seurat_obj@meta.data$patient = sample_name

#ALEX: 2025-03-12 - Added FOV QC steps to bring our pipeline in line with nanostring recommendations to remove bad FOVs
###################### FOV QC ####################
##################################################
#load Nanostring FOV QC functions
source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/Main/_code/FOV%20QC/FOV%20QC%20utils.R")

#load probe barcodes
allbarcodes <- readRDS(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))

#select only barcodes from 1000-plex human dataset
barcodemap <- allbarcodes$Hs_UCC

#create count matrix, set rownames and colnames to be the cell ids and gene names
counts <- seurat_obj@assays$Nanostring@layers$counts %>% t()
rownames(counts) <- rownames(seurat_obj@meta.data)
colnames(counts) <- rownames(seurat_obj@assays[["Nanostring"]]@features)

#obtain xy coordinates
xy <- seurat_obj@meta.data[, c("CenterX_local_px", "CenterY_local_px")] * 120.280945/1000000
colnames(xy) <- c("x_slide_mm", "y_slide_mm")

#withdraw FOV labels
fov <- seurat_obj@meta.data$FOV

#run FOV QC
res <- runFOVQC(counts = counts, xy = xy, fov = fov, barcodemap = barcodemap)

#plot flagged FOVs
mapFlaggedFOVs(res)

# list FOVs flagged for any reason, for loss of signal, for bias:
flaggedfov <- res$flaggedfovs

#list FOVs flagged for total counts
res$flaggedfovs_fortotalcounts

#list FOVs flagged for bias
res$flaggedfovs_forbias

# identify genes are involved in the flagged bits in those FOVs:
head(res$flagged_fov_x_gene)
head(res$flagged_fov_x_gene[, "gene"])

# count how many genes were impacted in one or more flagged FOVs:
length(unique(res$flagged_fov_x_gene[, "gene"]))

#ALEX: 2025-03-12 - added step to calculate proportion of negative probes per cell as this is recommended by NanoString
#calculate proportion of negative counts per cell
seurat_obj@meta.data$prop_negCount <- seurat_obj@meta.data$nCount_negprobes/seurat_obj@meta.data$nCount_Nanostring

#LEN: 2025-01-14- moving this cell QC filtering from 2_insitutype_single_sample_processing_pipeline.R to here
#ALEX: 2025-03-12 - added extra filter of removing cells with greater than 10% of probes representing negative probes (or whatever target is set), area > 3 standard devs, and bad FOVs with loss of signal or biased expression

# remove cells with less than target number of counts and 1 or more negative probe detected
seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= min_nCount_RNA &  prop_negCount <= max_proportion_negprobes & Area < 3*sd(seurat_obj@meta.data$Area))
`%notin%` <- Negate(`%in%`) #create not in operator

#remove the flagged FOVs
seurat_obj <- subset(seurat_obj, subset = FOV %notin% flaggedfov)

# remove all negative and falsecode probes
seurat_obj <- subset(seurat_obj, features = rownames(seurat_obj@assays[["Nanostring"]])[!grepl("SystemControl|Negative", rownames(seurat_obj@assays[["Nanostring"]]))])


saveRDS(seurat_obj, output_Rds_file_path)

############################################
#################### QC ####################
############################################
ggplot(seurat_obj@meta.data, aes(x=nCount_RNA)) + 
  geom_histogram(binwidth = 100, colour="black", fill = "grey") +
  geom_density(aes(y=..count..*100), color = "#000000", fill = "#F85700", alpha=.2) +
  geom_vline(aes(xintercept = 20, fill = "black")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme_classic() + 
  xlab("No. of Probes Detected") +
  ylab("No. of Cells") +
  xlim(min = 0, max = 3000)
ggsave(probe_count_histogram_file_path)

ggplot(seurat_obj@meta.data, aes(x=nFeature_RNA)) + 
  geom_histogram() +
  theme_classic() +
  xlab("No. of Genes Detected") +
  ylab("No. of Cells")
ggsave(feature_count_histogram_file_path)

ggplot(seurat_obj@meta.data, aes(x=prop_negCount)) + 
  geom_histogram(binwidth = 0.1) +
  theme(axis.text.y = element_text(size = 10)) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_classic() +
  xlab("Proportion of Negative Probes") +
  ylab("No. of Cells")
ggsave(negative_probe_proportion_histogram_file_path)

ggplot(seurat_obj@meta.data, aes(x=Area, y = ..count..)) +
  geom_histogram(binwidth = 1000, colour="black", fill = "grey") +
  geom_density(aes(y=..count..*1000), color = "#000000", fill = "#F85700", alpha=.2) +
  geom_vline(aes(xintercept = 4*sd(seurat_obj@meta.data$Area), fill = "black")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme_classic() + 
  xlab("Cell Area (pixels)") +
  ylab("No. of Cells") +
  xlim(min = 0, max = 50000)
ggsave(cell_area_histogram_file_path)

#plot FOV signal loss plot
FOVSignalLossSpatialPlot(res, outdir = output_folder_path, shownames = TRUE, plotwidth=6, plotheight=5)
FOVEffectsSpatialPlots(res = res, outdir = output_folder_path, bits = "all", plotwidth=6, plotheight=5)



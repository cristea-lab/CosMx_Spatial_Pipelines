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

if (length(args) != 3) {
    stop("Please double check your input to make sure only 3 command-line inputs are given.")
}

##################################################
#################### INPUT #######################
##################################################

sample_name = args[1]
print(paste("Input sample name: ", sample_name))


sample_folder_path = args[2]
print(paste("Input sample folder path from AtoMx: ", sample_folder_path))

output_folder_path = args[3]
print(paste("Input output folder: ", output_folder_path))

##################################################
#################### OUTPUT ######################
##################################################
meta_data_file_path = paste0(sample_folder_path, "/", sample_name, "_metadata_file.csv")

output_Rds_file_path = paste0(output_folder_path, "/A_1_pre_QC_sample.rds")

# QC
probe_count_histogram_file_path = paste0(output_folder_path, "/1_1_QC_probe_count_histogram.png")
feature_count_histogram_file_path = paste0(output_folder_path, "/1_2_QC_feature_count_histogram.png")
negative_probe_count_histogram_file_path = paste0(output_folder_path, "/1_3_QC_negative_probe_count_histogram.png")

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

saveRDS(seurat_obj, output_Rds_file_path)


############################################
#################### QC ####################
############################################
ggplot(seurat_obj@meta.data, aes(x=nCount_Nanostring)) + 
  geom_histogram() +
  theme_classic() +
  xlab("No. of Probes Detected") +
  ylab("No. of Cells") +
  xlim(min = 0, max = 4500)
ggsave(probe_count_histogram_file_path)

ggplot(seurat_obj@meta.data, aes(x=nFeature_Nanostring)) + 
  geom_histogram() +
  theme_classic() +
  xlab("No. of Genes Detected") +
  ylab("No. of Cells")
ggsave(feature_count_histogram_file_path)

ggplot(seurat_obj@meta.data, aes(x=nCount_negprobes)) + 
  geom_histogram() +
  theme_classic() +
  xlab("No. of Negative Probes Detected") +
  ylab("No. of Cells")
ggsave(negative_probe_count_histogram_file_path)
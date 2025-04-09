library("Seurat")
library("sceasy")
library("reticulate")
library("dplyr")
library("SpatialDecon")
library("InSituType")
library("tidyverse")
library("here")
library("dittoSeq")


#calculate fractional abundance and densities for level 1 annotations and create excel spreadsheet

###############################################################
#################### COMMAND LINE INPUT #######################
###############################################################
## Example command-line input: 
# Rscript 3_meta_data_excel_sheet.R PCA429 /cristealab/xiwang/Outputs/CosMx_Spatial_Pipelines/output_PCA429_08282024/InSituType/12_clusters/B_2_semi_sup_insitutype_fully_labeled.rds /cristealab/xiwang/Outputs/CosMx_Spatial_Pipelines/output_PCA429_08282024/InSituType/metrics_summary.csv


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

sample_rds_path = args[2]
print(paste("Input sample rds path: ", sample_rds_path))

output_csv = args[3]
print(paste("Output csv path: ", output_csv))



#-----------------------------------------------------------------------------------------------------

##################################################
#################### START #######################
##################################################



#################### LOAD SAMPLE ####################
seurat_obj = readRDS(sample_rds_path)


#calculate fractional abundance and densities for level 1 annotations and create excel spreadsheet
level_1_table <- table(seurat_obj$level_1_clusters, seurat_obj@meta.data$patient) %>%
  as.data.frame()
colnames(level_1_table) <- c("Cell_Type","Patient","Counts") #specify colnames

#create column to store fractional abundances and densities for each cluster
level_1_table[ ,c("Fraction")] <- NA 
level_1_table[ ,c("Density")] <- NA
#level_1_table[, c("Total_Cell_Num")] <-NA


summed_area <- seurat_obj@meta.data %>%
  group_by(patient, level_1_clusters) %>%
  summarise(
    #treatment_type = first(treatment_type),
   #treatment_name = first(treatment_name),
    area = sum(Area.um2),
    .groups = 'drop'  # This argument is used to ungroup after summarising
  )

#get number of cells for each slide 
num_cells <- nrow(seurat_obj@meta.data[seurat_obj@meta.data$patient == sample_name,])


#fill in fractional abundances for each cluster type
level_1_table[level_1_table$Patient == sample_name,]$Fraction <- level_1_table[level_1_table$Patient == sample_name,]$Counts/num_cells


#sort dataframes before binding
level_1_table <- level_1_table[order(level_1_table$Patient, level_1_table$Cell_Type),]
summed_area <- summed_area[order(summed_area$patient, summed_area$level_1_clusters),]

#add area column
level_1_table <- cbind(level_1_table, summed_area$area)

#calculate densities
level_1_table$Density <- level_1_table$Counts/level_1_table$`summed_area$area`*1000000

colnames(level_1_table) <- c("Cell_Type","Patient","Counts","Fraction","Density","Area") #specify colnames




###################################################################################################
#calculate fractional abundance and densities for level 2 annotations and create excel spreadsheet
level_2_table <- table(seurat_obj$level_2_clusters, seurat_obj@meta.data$patient) %>%
  as.data.frame()
colnames(level_2_table) <- c("Cell_Type","Patient","Counts") #specify colnames

#create column to store fractional abundances and densities for each cluster
level_2_table[ ,c("Fraction")] <- NA 
level_2_table[ ,c("Density")] <- NA
#level_2_table[, c("Total_Cell_Num")] <-NA


summed_area <- seurat_obj@meta.data %>%
  group_by(patient, level_2_clusters) %>%
  summarise(
    #treatment_type = first(treatment_type),
   #treatment_name = first(treatment_name),
    area = sum(Area.um2),
    .groups = 'drop'  # This argument is used to ungroup after summarising
  )

#get number of cells for each slide 
num_cells <- nrow(seurat_obj@meta.data[seurat_obj@meta.data$patient == sample_name,])


#fill in fractional abundances for each cluster type
level_2_table[level_2_table$Patient == sample_name,]$Fraction <- level_2_table[level_2_table$Patient == sample_name,]$Counts/num_cells

#sort dataframes before binding
level_2_table <- level_2_table[order(level_2_table$Patient, level_2_table$Cell_Type),]
summed_area <- summed_area[order(summed_area$patient, summed_area$level_2_clusters),]

#remove rows with zero counts
level_2_table <- level_2_table[level_2_table$Counts != 0,]

#add area column
level_2_table <- cbind(level_2_table, summed_area$area)

#calculate densities
level_2_table$Density <- level_2_table$Counts/level_2_table$`summed_area$area`*1000000

colnames(level_2_table) <- c("Cell_Type","Patient","Counts","Fraction","Density","Area") #specify colnames


###################################################################################################
#calculate fractional abundance and densities for level 3 annotations and create excel spreadsheet
level_3_table <- table(seurat_obj$level_3_clusters, seurat_obj@meta.data$patient) %>%
  as.data.frame()
colnames(level_3_table) <- c("Cell_Type","Patient","Counts") #specify colnames

#create column to store fractional abundances and densities for each cluster
level_3_table[ ,c("Fraction")] <- NA 
level_3_table[ ,c("Density")] <- NA
#level_3_table[, c("Total_Cell_Num")] <-NA

#get

summed_area <- seurat_obj@meta.data %>%
  group_by(patient, level_3_clusters) %>%
  summarise(
    #treatment_type = first(treatment_type),
   #treatment_name = first(treatment_name),
    area = sum(Area.um2),
    .groups = 'drop'  # This argument is used to ungroup after summarising
  )

#get number of cells for each slide 
num_cells <- nrow(seurat_obj@meta.data[seurat_obj@meta.data$patient == sample_name,])


#fill in fractional abundances for each cluster type
level_3_table[level_3_table$Patient == sample_name,]$Fraction <- level_3_table[level_3_table$Patient == sample_name,]$Counts/num_cells

#sort dataframes before binding
level_3_table <- level_3_table[order(level_3_table$Patient, level_3_table$Cell_Type),]
summed_area <- summed_area[order(summed_area$patient, summed_area$level_3_clusters),]

#remove rows with zero counts
level_3_table <- level_3_table[level_3_table$Counts != 0,]

#add area column
level_3_table <- cbind(level_3_table, summed_area$area)

#calculate densities
level_3_table$Density <- level_3_table$Counts/level_3_table$`summed_area$area`*1000000

colnames(level_3_table) <- c("Cell_Type","Patient","Counts","Fraction","Density","Area") #specify colnames




###################################################################################################
#calculate fractional abundance and densities for tumor cells and create excel spreadsheet
tumor_cells <- subset(seurat_obj, , subset = level_1_clusters %in% c("Tumor"))

level_3_table_tumor <- table(tumor_cells$level_3_clusters, tumor_cells@meta.data$patient) %>%
  as.data.frame()
colnames(level_3_table_tumor) <- c("Cell_Type","Patient","Counts") #specify colnames

#create column to store fractional abundances and densities for each cluster
level_3_table_tumor[ ,c("Fraction")] <- NA 
level_3_table_tumor[ ,c("Density")] <- NA
#level_3_table[, c("Total_Cell_Num")] <-NA

#get

summed_area <- tumor_cells@meta.data %>%
  group_by(patient, level_3_clusters) %>%
  summarise(
    #treatment_type = first(treatment_type),
   #treatment_name = first(treatment_name),
    area = sum(Area.um2),
    .groups = 'drop'  # This argument is used to ungroup after summarising
  )

#get number of cells for each slide 
num_cells <- nrow(tumor_cells@meta.data[tumor_cells@meta.data$patient == sample_name,])


#fill in fractional abundances for each cluster type
level_3_table_tumor[level_3_table_tumor$Patient == sample_name,]$Fraction <- level_3_table_tumor[level_3_table_tumor$Patient == sample_name,]$Counts/num_cells


#sort dataframes before binding
level_3_table_tumor <- level_3_table_tumor[order(level_3_table_tumor$Patient, level_3_table_tumor$Cell_Type),]
summed_area <- summed_area[order(summed_area$patient, summed_area$level_3_clusters),]

#remove rows with zero counts
level_3_table_tumor <- level_3_table_tumor[level_3_table_tumor$Counts != 0,]

#add area column
level_3_table_tumor <- cbind(level_3_table_tumor, summed_area$area)

#calculate densities
level_3_table_tumor$Density <- level_3_table_tumor$Counts/level_3_table_tumor$`summed_area$area`*1000000

colnames(level_3_table_tumor) <- c("Cell_Type","Patient","Counts","Fraction","Density","Area") #specify colnames





###################################################################################################
#calculate fractional abundance and densities for fibroblasts and create excel spreadsheet
fibroblasts <- subset(seurat_obj, subset = level_2_clusters %in% c("myCAF","iCAF","Mesenchymal.Cell"))

level_2_table_fibroblasts <- table(fibroblasts$level_2_clusters, fibroblasts@meta.data$patient) %>%
  as.data.frame()
colnames(level_2_table_fibroblasts) <- c("Cell_Type","Patient","Counts") #specify colnames

#create column to store fractional abundances and densities for each cluster
level_2_table_fibroblasts[ ,c("Fraction")] <- NA 
level_2_table_fibroblasts[ ,c("Density")] <- NA
#level_2_table_fibroblasts[, c("Total_Cell_Num")] <-NA

#get

summed_area <- fibroblasts@meta.data %>%
  group_by(patient, level_2_clusters) %>%
  summarise(
    #treatment_type = first(treatment_type),
   #treatment_name = first(treatment_name),
    area = sum(Area.um2),
    .groups = 'drop'  # This argument is used to ungroup after summarising
  )

#get number of cells for each slide 
num_cells <- nrow(fibroblasts@meta.data[fibroblasts@meta.data$patient == sample_name,])


#fill in fractional abundances for each cluster type
level_2_table_fibroblasts[level_2_table_fibroblasts$Patient == sample_name,]$Fraction <- level_2_table_fibroblasts[level_2_table_fibroblasts$Patient == sample_name,]$Counts/num_cells

#sort dataframes before binding
level_2_table_fibroblasts <- level_2_table_fibroblasts[order(level_2_table_fibroblasts$Patient, level_2_table_fibroblasts$Cell_Type),]
summed_area <- summed_area[order(summed_area$patient, summed_area$level_2_clusters),]

#remove rows with zero counts
level_2_table_fibroblasts <- level_2_table_fibroblasts[level_2_table_fibroblasts$Counts != 0,]

#add area column
level_2_table_fibroblasts <- cbind(level_2_table_fibroblasts, summed_area$area)

#calculate densities
level_2_table_fibroblasts$Density <- level_2_table_fibroblasts$Counts/level_2_table_fibroblasts$`summed_area$area`*1000000

colnames(level_2_table_fibroblasts) <- c("Cell_Type","Patient","Counts","Fraction","Density","Area") #specify colnames




### OUTPUT ALL DATA INTO REQUESTED FORMAT ###

# Initialize a list to store the metrics
metrics <- list(
  area.epithelial = NA,
  area.stroma = NA,
  area.overall = NA,
  cells.epithelial = NA,
  cells.stroma = NA,
  cells.overall = NA,
  fraction.cells.tumor = NA,
  fraction.cells.immune = NA,
  fraction.cells.fibroblast = NA,
  fraction.cells.other = NA,
  density.overall.cd3 = NA,
  density.overall.cd3cd4 = NA,
  density.overall.cd3cd8 = NA,
  density.overall.m1macrophage = NA,
  density.overall.m2macrophage = NA,
  density.overall.macrophage = NA,
  density.overall.cd14 = NA,
  density.overall.totalimmune = NA,
  density.overall.totalimmuneprolif = NA,
  fraction.tumor.basal = NA,
  fraction.tumor.classical = NA,
  fraction.tumor.coexpressor = NA,
  fraction.tumor.prolif = NA,
  fraction.fibroblast.mycaf = NA,
  fraction.fibroblast.icaf = NA
)

# Fill in the metrics from the tables
metrics$area.epithelial <- level_1_table %>% filter(Cell_Type == "Tumor") %>% pull(Area)
metrics$area.stroma <- level_1_table %>% filter(Cell_Type == "Stroma") %>% pull(Area)
metrics$area.overall <- sum(level_1_table$Area)

metrics$cells.epithelial <- level_1_table %>% filter(Cell_Type == "Tumor") %>% pull(Counts)
metrics$cells.stroma <- level_1_table %>% filter(Cell_Type == "Stroma") %>% pull(Counts)
metrics$cells.overall <- sum(level_1_table$Counts)

metrics$fraction.cells.tumor <- level_1_table %>% filter(Cell_Type == "Tumor") %>% pull(Fraction)
metrics$fraction.cells.immune <- level_1_table %>% filter(Cell_Type == "Immune") %>% pull(Fraction)
metrics$fraction.cells.fibroblast <- level_2_table %>% filter(Cell_Type %in% c("iCAF", "Mesenchymal.Cell", "myCAF")) %>% summarise(Fraction = sum(Fraction)) %>% pull(Fraction)
metrics$fraction.cells.other <- 1 - (metrics$fraction.cells.tumor + metrics$fraction.cells.immune + metrics$fraction.cells.fibroblast)

# Density metrics (assuming density is given in the tables)
metrics$density.overall.cd3 <- NA # Not provided in the tables
metrics$density.overall.cd3cd4 <- NA  # Not provided in the tables
metrics$density.overall.cd3cd8 <- NA  # Not provided in the tables
metrics$density.overall.m1macrophage <- level_3_table %>% filter(Cell_Type == "Macrophage.1") %>% pull(Density)
metrics$density.overall.m2macrophage <- level_3_table %>% filter(Cell_Type == "Macrophage.2") %>% pull(Density)
metrics$density.overall.macrophage <- level_3_table %>% filter(Cell_Type %in% c("Macrophage.1", "Marcophage.2")) %>% summarise(Density = sum(Density)) %>% pull(Density)
metrics$density.overall.cd14 <- NA  # Not provided in the tables
metrics$density.overall.totalimmune <- level_1_table %>% filter(Cell_Type == "Immune") %>% pull(Density)
metrics$density.overall.totalimmuneprolif <- NA  # Not provided in the tables

# Tumor fractions (assuming these are subtypes of tumor)
metrics$fraction.tumor.basal <- NA  # Not provided in the tables
metrics$fraction.tumor.classical <- NA  # Not provided in the tables
metrics$fraction.tumor.coexpressor <- NA  # Not provided in the tables
metrics$fraction.tumor.prolif <- NA  # Not provided in the tables

# Fibroblast fractions
metrics$fraction.fibroblast.mycaf <- level_2_table_fibroblasts %>% filter(Cell_Type == "myCAF") %>% pull(Fraction)
metrics$fraction.fibroblast.icaf <- level_2_table_fibroblasts %>% filter(Cell_Type == "iCAF") %>% pull(Fraction)

# Convert the list to a data frame
metrics_df <- data.frame(
  Metric = names(metrics),
  Value = unlist(metrics)
)

# # Write the data frame to a CSV file
write.csv(metrics_df, output_csv, row.names = FALSE)

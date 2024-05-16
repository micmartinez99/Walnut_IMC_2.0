# The purpose of this script is to try and get marker information for certain ROIs to relate the levels of free urolithin A to them


# Load libraries
library(imcRtools)
library(cytomapper)
library(RColorBrewer)
library(dittoSeq)
library(viridis)
library(gridGraphics)
library(cowplot)
library(scuttle)
library(dplyr)
library(pheatmap)
library(lisaClust)
library(CATALYST)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(scater)
library(scales)
library(BiocParallel)

spe <- readRDS("DATA/Revised_Rds_Objects/FULLY_REVISED_Walnut_Spe_Object_4.16.24.Rds")

# Calculate the mean counts for each ROI
image_mean <- aggregateAcrossCells(spe[rowData(spe)$use_channel], 
                                   ids = spe$roi_identifier_from_dr_moussa, 
                                   statistics="mean",
                                   use.assay.type = "counts")

# Transform the counts to asinh expression
assay(image_mean, "exprs") <- asinh(counts(image_mean))

# Get the expression information as a dataframe
expr <- as.data.frame(t(assay(image_mean, "exprs")))
write.csv(expr, file = "Walnut_IMC_Manuscript_Revisions_Outputs/Data_Files/asinh_expression_levels_Per_ROI_Per_Marker.csv")

# Get the expression levels for each individual cell so we can plot violins
cells <- as.data.frame(t(assay(spe, "exprs")))
ROIS <- spe$roi_identifier_from_dr_moussa
cells$ROIs <- ROIS
cells[,c(1,3:29)] <- NULL
cells$ROIs <- gsub("ROI ", "\\", cells$ROIs)
cells <- cells %>%
  mutate(Group = case_when(
    ROIs >= 1 & ROIs <= 10 ~ "Urolithin Low",
    ROIs >= 11 & ROIs <= 23 ~ "Urolithin High",
    TRUE ~ "Other"
  ))
cells$ROIs <- paste("ROI", cells$ROIs, sep = " ")
write.csv(cells, file = "Walnut_IMC_Manuscript_Revisions_Outputs/Data_Files/asinh_expression_levels_Per_ROI_Per_Marker.csv")






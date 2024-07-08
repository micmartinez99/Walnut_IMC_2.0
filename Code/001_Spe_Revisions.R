# The purpose of this script is to regenerate some of the figures for the walnut manuscript 


# Clear environment
rm(list = ls())

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

# Set working directory

# Read in the spe object
spe <- readRDS("DATA/Pratik/study_combined_spe_cluster_subset_by_WalnutIMCCombined_ordered_all_14ClustersLabels.csv_file.rds")

# Transform the counts
assay(spe, "exprs") <- asinh(counts(spe)/1)
rowData(spe)$use_channel <- !grepl("DNA1|DNA3", rownames(spe))

# Assign the labels to the clusters
cluster_celltype <- recode(spe$clusters_as_character_vector,
                           "C1" = "Slow (middle-bottom) epithelia",
                           "C2" = "TA-like/middle epithelia",
                           "C3" = "Antigen presenting cells",
                           "C4" = "Slow middle epithelia", 
                           "C5" = "Various immune",
                           "C6" = "Slow middle epithelia",
                           "C7" = "Surface epithelia",
                           "C8" = "Musculuar cells",
                           "C9" = "Bottom epithelia",
                           "C10" = "Antigen presenting cells",
                           "C11" = "Uncommon surface epithelia",
                           "C12" = "Non-specific stroma",
                           "C13" = "Myeloid lineage",
                           "C14" = "Fibrous stromal cells")

# Assign phenotype labels and assign Type information
spe$Phenotypes <- cluster_celltype
spe$Type <- ifelse(spe$Phenotypes == "Slow middle epithelia" | spe$Phenotypes == "Uncommon surface epithelia" |
                     spe$Phenotypes == "Surface epithelia" | spe$Phenotypes == "Bottom epithelia" | spe$Phenotypes == "TA-like/middle epithelia" |
                     spe$Phenotypes == "Slow (middle-bottom) epithelia", "Epithelial clusters", "Non-epithelial clusters")

# Calculate the mean expression of each Phenotype
celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$Phenotypes, 
                                      statistics = "mean",
                                      use.assay.type = "exprs")


# Assign cluster colors and type colors
clustersColors <- setNames(brewer.pal(length(unique(spe$Phenotypes)), name = "Paired"),
                           unique(spe$Phenotypes))
TypeColors <- setNames(c("cyan3", "pink"), c("Epithelial clusters", "Non-epithelial clusters"))

# Create a new empty list to hold the cluster color vecrors
colorings <- list()

# Append the color vectors to the colorings list
colorings$PhenotypeColors <- clustersColors
colorings$TypeColors <- TypeColors

# Append the colorings list to the metadata slot of spe
metadata(spe)$colors <- colorings

# Plot celltype mean heatmap, min max scaled
New <- dittoHeatmap(celltype_mean,
             assay = "exprs", 
             genes = rownames(spe)[rowData(spe)$use_channel],
             cluster_cols = TRUE, 
             cluster_rows = TRUE,
             scaled.to.max = TRUE,
             heatmap.colors.max.scaled = viridis(100),
             annot.by = c("Phenotypes", "Type"),
             annotation_colors = list(Phenotypes = metadata(spe)$colors$PhenotypeColors,
                                      Type = metadata(spe)$colors$TypeColors))
ggsave("New_Celltype_Mean_Heatmap.tiff", New, width = 8, height = 8, dpi = 300)

# Save rds file
saveRDS(spe, file = "DATA/Revised_Rds_Objects/April_2024_Edited_by_Mike_Walnut_Spe.Rds")


# Now we can visualize some of these clusters on the tissues
# Read in the images and masks
img <- readRDS("DATA/Pratik/study_combined_images_normalized_by___separateImages_=_FALSE___separateChannels_=_TRUE_2023_04_25_and_then_followed_by_normalize_0_and_0.2.rds")
masks <- readRDS("DATA/Pratik/study_combined_masks.rds")
saveRDS(img, file = "DATA/Revised_Rds_Objects/April_2024_Edited_by_Mike_asinh_Images.Rds")
saveRDS(masks, file = "DATA/Revised_Rds_Objects/April_2024_Edited_by_Mike_Walnut_Masks.Rds")

# Ensure that channelNames of the images match with the rownames of the spatial experiment object
channelNames(img) <- rownames(spe)

# Set image names across images/masks and spe
mcols(img) <- mcols(masks) <- DataFrame(sample_id = names(img),
                                        ROI = paste("ROI", unique(spe$roi_identifier_from_dr_moussa), sep = " "))

# Visualize clusters on the tissues
tissues <- names(img)

# Within a for loop, visualize each images with the ROI as the title
for(i in tissues) {
  cur_img <- img[i]
  cur_mask <- masks[i]
  
  projections <- plotCells(cur_mask,
            object = spe,
            cell_id = "ObjectNumber",
            img_id = "sample_id",
            image_title = list(text = mcols(cur_img)$ROI,
                               position = "topright",
                               font = 0),
            colour_by = "Type",
            colour = list(Phenotypes = metadata(spe)$colors$PhenotypeColors),
            save_plot = list(filename = paste(mcols(cur_img)$ROI, "cluster_projections.tiff", sep = "_")))
  
}

# I want to visualize just a single cluster on the images
for (i in tissues){
  
  # Define the current ROI and its associated mask
  cur_img <- img[i]
  cur_mask <- masks[i]
  
  for (j in unique(spe$Phenotypes)){
    
    # Subset the spe object to be just the cluster of the current iteration
    clust <- spe[,spe$Phenotypes == j]
    
    # Format "j" for filename
    j <- gsub("/", "_", j)
    
    # Plot the pixels of just this cluster for the current ROI
    clustCells <- plotCells(mask = cur_mask,
                              object = clust,
                              cell_id = "ObjectNumber",
                              img_id = "sample_id",
                              colour_by = "Phenotypes",
                              colour = list(Phenotypes = metadata(clust)$colors$PhenotypeColors),
                              thick = TRUE,
                              image_title = list(text = mcols(cur_img)$ROI,
                                                 position = "topright",
                                                 font = 0),
                              save_plot = list(filename = paste(mcols(cur_img)$ROI, j, "cluster_Projections.tiff", sep = "_")))
    rm(clust)
  }
}

# Read in the non-normalized images
nonNorm_img <- readRDS("DATA/Pratik/study_combined_images.rds")

# Set image names across images/masks and spe
mcols(nonNorm_img) <- mcols(masks) <- DataFrame(sample_id = names(nonNorm_img),
                                        ROI = paste("ROI", unique(spe$roi_identifier_from_dr_moussa), sep = " "))

# This is for visualizing composite images of marker stain and the cell masks
 for(i in tissues) {
   nncur_img <- nonNorm_img[i]
   cur_mask <- masks[i]
   
   projections <- plotPixels(nncur_img,
                            cur_mask,
                            object = spe,
                            cell_id = "ObjectNumber",
                            img_id = "sample_id",
                            image_title = list(text = mcols(cur_img)$ROI,
                                               position = "topright",
                                               font = 0),
                           colour_by = c("CD20", "CD8a"),
                            outline_by = "Phenotypes",
                            bcg = list(`CD20` = c(0,5,1),
                                       CD8a = c(0,5,1)),
                            colour = list(Phenotypes = metadata(clust)$colors$PhenotypeColors),
                           save_plot = list(filename = paste(mcols(nncur_img)$ROI, "cluster_projections.tiff", sep = "_")))
   
 }
 












# walnut_tumor_low_producer_01_01_2023_mcd1_roi1
# Build spatial graphs
spe <- readRDS("DATA/Revised_Rds_Objects/April_2024_Edited_by_Mike_Walnut_Spe.Rds")
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "knn", k = 20)
# Make to 6


# Spatial Community Analysis
set.seed(03061999)
spe <- aggregateNeighbors(spe, 
                          colPairName = "knn_interaction_graph", 
                          aggregate_by = "metadata", 
                          count_by = "Phenotypes")


set.seed(03061999)
cn_1 <- kmeans(spe$aggregatedNeighbors, centers = 6)
spe$cn_celltypes <- as.factor(cn_1$cluster)

# Examine the phenotype composition of the detected cellular neighborhoods. Look at the total number of cells per phenotype and the CN
for_plot1 <- table(as.character(spe$cn_celltypes), spe$Phenotypes)

# Plot total number of phenotypes in the CNs
pheatmap(for_plot1, 
         color = viridis(100), display_numbers = TRUE, 
         number_color = "white", number_format = "%.0f")


# Next, observe (per phenotype) the fraction of CN that they are distributed across
for_plot2 <- prop.table(table(as.character(spe$cn_celltypes), spe$Phenotypes), margin = 2)

pheatmap(for_plot2, 
         color = viridis(100), display_numbers = TRUE, 
         number_color = "white", number_format = "%.2f",
         scale = "column")

# Lastly, we can visualize the fraction of each CN made up of each cell type
for_plot3 <- prop.table(table(as.character(spe$cn_celltypes), spe$Phenotypes), margin = 1)

# Scale the columns to account for the relative phenotype abundances
pheatmap(for_plot3, 
         color = viridis(100), display_numbers = TRUE, 
         number_color = "white", number_format = "%.2f", 
         scale = "column")

# Visualize the enrichment of cell pehnotypes within cellular neighborhoods
regionMap(spe, 
          cellType = "Phenotypes",
          region = "cn_celltypes")

# Plot the spatial neighborhoods for wach tissue
for (i in tissues) {
  cur_img <- img[i]
  cur_mask <- masks[i]
  
  neighbors <- plotSpatial(spe[,spe$sample_id == i],
              img_id = "cn_celltypes", node_color_by = "Phenotypes", node_size_fix = 0.7) +
    scale_color_manual(values = metadata(spe)$colors$PhenotypeColors)
  ggsave(paste(mcols(cur_img)$ROI, "Neighbors.tiff", sep = "_"))
  
}


# We can essentially do the same thing by computing the mean expression across the 20-nearest neighbor prior to kmeans clustering (k = 6)
spe <- aggregateNeighbors(spe,
                          colPairName = "knn_interaction_graph",
                          aggregate_by = "expression",
                          assay_type = "exprs",
                          subset_row = rowData(spe)$use_channel)
set.seed(03061999)

# Cluster to 6 CNs
cn_2 <- kmeans(spe$mean_aggregatedExpression, centers = 6)
spe$cn_expression <- as.factor(cn_2$cluster)

for (i in tissues) {
  cur_img <- img[i]
  cur_mask <- masks[i]
  
  neighbors <- plotSpatial(spe[,spe$sample_id == i],
                           img_id = "sample_id", node_color_by = "cn_expression", node_size_fix = 0.7) 
  ggsave(paste(mcols(cur_img)$ROI, "Neighbors.tiff", sep = "_"))
  
}

# We can visualize the phenotype composition of each cellular neighborhood
for_plot4 <- prop.table(table(spe$cn_expression, spe$Phenotypes), 
                       margin = 1)

pheatmap(for_plot4, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")


# CATALYST-based visualiztion
library(CATALYST)

# Save SPE inCATALYST-compatible object with renamed colData entries and new metadata information
spe_cat <- spe
spe_cat$sample_id <- factor(paste("ROI", spe$roi_identifier_from_dr_moussa, sep = " "))
spe_cat$cluster_id <- factor(spe$Phenotypes)
spe_cat$condition <- factor(ifelse(spe$low_or_high_urolithin_producer == "low_producer", "Low", "High"), 
                            levels = c("Low", "High"))

# Add celltype information to metadata
metadata(spe_cat)$cluster_codes <- data.frame(celltype = factor(spe_cat$cluster_id))

# MDS pseudobulk by cell type
MDSplot <- pbMDS(spe_cat, 
      by = "cluster_id", 
      features = rownames(spe_cat)[rowData(spe_cat)$use_channel], 
      label_by = "cluster_id", 
      k = "celltype") +
  scale_color_manual(values = metadata(spe)$colors$PhenotypeColors) +
  theme_bw() +
  labs(size = "n cells",
       color = "Phenotypes",
       fill = "Phenotypes")
ggsave("Cluster_NMDS.tiff", MDSplot, width = 10, height = 10, dpi = 300)
  


# MDS pseudobulk by cell type and sample_id
pbMDS(spe_cat, 
      by = "both", 
      features = rownames(spe_cat)[rowData(spe_cat)$use_channel], 
      k = "celltype", 
      shape_by = "condition", 
      size_by = TRUE) +
  scale_color_manual(values = metadata(spe)$colors$PhenotypeColors) +
  theme_bw()
  
# CLR on cluster proportions across samples
PCA <- clrDR(spe_cat, 
      dr = "PCA", 
      by = "sample_id", 
      k = "celltype", 
      label_by = "sample_id", 
      arrow_col = "cluster_id") +
  scale_color_manual(values = metadata(spe)$colors$PhenotypeColors) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(color = "Phenotypes",
       size = "n cells",
       fill = "Producer")
ggsave("Low_vs_High_Clusters_PCA.png", PCA, width = 10, height = 10, dpi = 300)





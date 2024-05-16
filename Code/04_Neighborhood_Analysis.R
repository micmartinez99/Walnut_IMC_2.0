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
library(CATALYST)
library(ComplexHeatmap)
library(circlize)
library(scater)

setwd("Walnut_IMC_Manuscript_Revisions_Outputs/")

# Read in the modified spe object and normalized images
spe <- readRDS("/users/michaelmartinez/Desktop/Walnut_IMC_2.0/Walnut_IMC_Manuscript_Revisions_Outputs/Revised_Rds_Objects/April_2024_Edited_by_Mike_Walnut_Spe.Rds")

# Read in the images and masks
img <- readRDS("/users/michaelmartinez/Desktop/Walnut_IMC_2.0/Walnut_IMC_Manuscript_Revisions_Outputs/Revised_Rds_Objects/April_2024_Edited_by_Mike_asinh_Images.Rds")
masks <- readRDS("Walnut_IMC_Manuscript_Revisions_Outputs/Revised_Rds_Objects/April_2024_Edited_by_Mike_Walnut_Masks.Rds")

# Re-classify the clusters (4/11/24)
cluster_celltype <- recode(spe$Phenotypes,
                           "Slow (middle-bottom) epithelia" = "Non-proliferating epithelia",
                           "TA-like/middle epithelia" = "Proliferating epithelia",
                           "Antigen presenting cells" = "Antigen presenting cells",
                           "Slow middle epithelia" = "Non-proliferating epithelia", 
                           "Various immune" = "Non-specific immune",
                           "Slow middle epithelia" = "Non-proliferating epithelia",
                           "Surface epithelia" = "Non-proliferating epithelia",
                           "Musculuar cells" = "Muscle layer",
                           "Bottom epithelia" = "Non-proliferating epithelia",
                           "Antigen presenting cells" = "Antigen presenting cells",
                           "Uncommon surface epithelia" = "Non-proliferating epithelia",
                           "Non-specific stroma" = "Stroma",
                           "Myeloid lineage" = "Myeloid lineage",
                           "Fibrous stromal cells" = "Antigen presenting cells")

# Add cluster information to spe object
spe$Phenotypes <- cluster_celltype

# Find the number of cells in each cluster
sum(spe$Phenotypes == "Non-proliferating epithelia") #151581
sum(spe$Phenotypes == "Proliferating epithelia") #37895
sum(spe$Phenotypes == "Antigen presenting cells") #70295
sum(spe$Phenotypes == "Non-specific immune") #36145
sum(spe$Phenotypes == "Muscle layer") #25934
sum(spe$Phenotypes == "Stroma") #37755
sum(spe$Phenotypes == "Myeloid lineage") #19586

cluster_cellNums <- recode(spe$Phenotypes,
                           "Non-proliferating epithelia" = "Non-proliferating epithelia (151581)",
                           "Proliferating epithelia" = "Proliferating epithelia (37895)",
                           "Antigen presenting cells" = "Antigen presenting cells (70295)",
                           "Non-specific immune" = "Non-specific immune (36145)",
                           "Muscle layer" = "Muscle layer (25934)",
                           "Stroma" = "Stroma (37755)",
                           "Myeloid lineage" = "Myeloid lineage (19586)")
spe$Phenotype <- cluster_cellNums
spe$Phenotype <- factor(spe$Phenotype, levels = c("Non-proliferating epithelia (151581)",
                                                  "Proliferating epithelia (37895)",
                                                  "Non-specific immune (36145)",
                                                  "Muscle layer (25934)",
                                                  "Stroma (37755)",
                                                  "Antigen presenting cells (70295)",
                                                  "Myeloid lineage (19586)"))


# Assign Type
spe$Type <- ifelse(spe$Phenotypes == "Proliferating epithelia"|spe$Phenotypes == "Non-proliferating epithelia", "Epithelial clusters", "Non-epithelial clusters")
spe$Subgroup <- ifelse(spe$Phenotype == "Proliferating epithelia (37895)"|spe$Phenotype == "Non-proliferating epithelia (151581)", "Epithelial clusters", "Non-epithelial clusters")

# Calculate the mean expression of each Phenotype
celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$Phenotypes, 
                                      statistics = "mean",
                                      use.assay.type = "exprs")


# Assign cluster colors and type colors
phenotypeColors <- setNames(brewer.pal(length(unique(spe$Phenotype)), name = "Paired"),
                            unique(spe$Phenotype))
TypeColors <- setNames(c("cyan3", "pink"), c("Epithelial clusters", "Non-epithelial clusters"))
NumColors <- setNames(c("cyan3", "pink"), c("Epithelial clusters", "Non-epithelial clusters"))


# Create a new empty list to hold the cluster color vecrors
colorings <- list()

# Append the color vectors to the colorings list
colorings$PhenotypeColors <- phenotypeColors
colorings$TypeColors <- TypeColors
colorings$NumColors <- NumColors

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
                    annot.by = c("Phenotype", "Subgroup"),
                    annotation_colors = list(Phenotype = metadata(spe)$colors$PhenotypeColors,
                                             Subgroup = metadata(spe)$colors$NumColors
                                             ))
ggsave("Relabelled_Celltype_Mean_Heatmap.tiff", New, width = 8, height = 8, dpi = 300)







spe$roi_identifier_from_dr_moussa <- paste("ROI", spe$roi_identifier_from_dr_moussa, sep = " ")
spe$roi_identifier_from_dr_moussa <- factor(spe$roi_identifier_from_dr_moussa, levels = c(
  "ROI 1", "ROI 2", "ROI 3", "ROI 4", "ROI 5", "ROI 6", "ROI 7", "ROI 8", "ROI 9", "ROI 10",
  "ROI 11", "ROI 12", "ROI 13", "ROI 14", "ROI 15", "ROI 16", "ROI 17", "ROI 18", "ROI 19", "ROI 20",
  "ROI 21", "ROI 22", "ROI 23"
))

spe$groups <- ifelse(spe$groups == "normal_high_producer", "High Normal",
                     ifelse(spe$groups == "normal_low_producer", "Low Normal",
                            ifelse(spe$groups == "tumor_high_producer", "High Polyp", "Low Polyp")))

spe$groups <- factor(spe$groups, levels = c("Low Normal", "Low Polyp", "High Normal", "High Polyp"))

image_mean <- aggregateAcrossCells(spe, 
                                   ids = spe$roi_identifier_from_dr_moussa, 
                                   statistics="mean",
                                   use.assay.type = "counts")
assay(image_mean, "exprs") <- asinh(counts(image_mean))

group_coloring <- c("Low Normal" = "#8DD3C7",
                  "Low Polyp" = "#BEBADA",
                  "High Normal" = "#FB8072",
                  "High Polyp" = "#80B1D3")
spe$Group <- spe$groups
# individiual ROI heatmap
ROIHeatmap <- dittoHeatmap(image_mean, genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", 
             scaled.to.max = TRUE,
             cluster_cols = FALSE,
             heatmap.colors.max.scaled = viridis(100),
             show_colnames = TRUE,
             annot.by = c("groups"),
             annotation_colors = list(groups = group_coloring),
             fontsize_col = 10,
             gaps_col = c(10))

ggsave("ROI_heatmap.tiff", ROIHeatmap, width = 8, height = 8, dpi = 300)





# Barplot visualization
barplot <- dittoBarPlot(spe,
             var = "Phenotype",
             group.by = c("roi_identifier_from_dr_moussa"),
             split.by = "groups",
             split.adjust = list(scales = "free"),
             #x.reorder = c(1,12,13,14,15,16,17,18,19,20,21,22,23,2,3,4,5,6,7,8,9,10,11),
             legend.show = FALSE) +
  scale_fill_manual(values = metadata(spe)$colors$PhenotypeColors) +
  labs(fill = "Phenotypes",
       title = "") +
  theme(strip.text = element_text(face = "bold", size = 14),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 0))
ggsave("Barplot_Composition_By_ROI.tiff", barplot, width = 10, height = 10, dpi = 300)


# Ensure that channelNames of the images match with the rownames of the spatial experiment object
channelNames(img) <- rownames(spe)

# Set image names across images/masks and spe
mcols(img) <- mcols(masks) <- DataFrame(sample_id = names(img),
                                        ROI = unique(spe$roi_identifier_from_dr_moussa))

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
                           colour_by = "Phenotype",
                           colour = list(Phenotype = metadata(spe)$colors$PhenotypeColors),
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


# ROI 11 and ROI 22
cur_img <- img[c(9,20)]
cur_mask <- masks[c(9,20)]

eleven1 <- plotCells(cur_mask,
                         object = spe,
                         cell_id = "ObjectNumber",
                         img_id = "sample_id",
                         image_title = list(text = mcols(cur_img)$ROI,
                                            position = "topright",
                                            font = 0),
                         colour_by = "Phenotype",
                         colour = list(Phenotype = metadata(spe)$colors$PhenotypeColors),
                     return_plot = TRUE)
eleven2 <- plotCells(cur_mask,
                     object = spe,
                     cell_id = "ObjectNumber",
                     img_id = "sample_id",
                     colour_by = c("CD8a", "CD3", "CD20", "Epcam+ECadherin"),
                     colour = list(CD8a = c("black", "red"),
                                   CD3 = c("black", "burlywood1"),
                                   CD20 = c("black", "cyan2"),
                                   `Epcam+ECadherin` = c("black", "white")),
                     exprs_values = "exprs",
                     image_title = list(text = mcols(cur_img)$ROI,
                                        position = "topright",
                                        font = 0),
                     legend = list(colour_by.title.cex = 1,
                                   margin = 0.1),
                     return_plot = TRUE)
library(cowplot)
library(gridGraphics)
p11.1 <- ggdraw(eleven1$plot, clip = "on")
p11.2 <- ggdraw(eleven2$plot, clip = "on")

ROI11 <- plot_grid(p11.1, p11.2)
ggsave("ROI11_Figure.tiff", ROI11, width = 10, height = 10, dpi = 300)
                     
                     
                         

# Visualize pixels within cell masks on the tissues
tissues <- names(img)

# Within a for loop, visualize each images with the ROI as the title
for(i in tissues) {
  cur_img <- img[i]
  cur_mask <- masks[i]
  
  plotCells(cur_mask,
            object = spe,
            cell_id = "ObjectNumber",
            img_id = "sample_id",
            colour_by = c("CD8a", "CD3", "CD20", "Epcam+ECadherin"),
            colour = list(CD8a = c("black", "red"),
                          CD3 = c("black", "burlywood1"),
                          CD20 = c("black", "cyan2"),
                          `Epcam+ECadherin` = c("black", "white")),
            exprs_values = "exprs",
            image_title = list(text = mcols(cur_img)$ROI,
                               position = "topright",
                               font = 0),
            legend = list(colour_by.title.cex = 5,
                          margin = 0.1),
            save_plot = list(filename = paste(mcols(cur_img)$ROI, ".tiff", sep = "_")))

}

# Save SPE inCATALYST-compatible object with renamed colData entries and new metadata information
spe_cat <- spe
spe_cat$sample_id <- spe$roi_identifier_from_dr_moussa
spe_cat$cluster_id <- factor(spe$Phenotype)
spe_cat$condition <- spe$groups

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
       color = "Phenotype",
       fill = "Phenotype")
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
  theme(legend.position = "bottom",
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  labs(color = "Phenotype",
       size = "n cells",
       fill = "Group")
ggsave("Group_Clusters_PCA.tiff", PCA, width = 10, height = 10, dpi = 300)

# Sample down to 1000 cells
cur_cells <- sample(seq_len(ncol(spe)), 1000)
violins <- plotExpression(spe[,cur_cells],
               features = rownames(spe)[rowData(spe)$use_channel],
               x = "Phenotypes",
               exprs_values = "exprs",
               colour_by = "Phenotypes",
               one_facet = FALSE) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 0),
        strip.text = element_text(face = "bold", size = 14)) +
  scale_color_manual(values = metadata(spe)$colors$PhenotypeColors) +
  labs(color = "Phenotypes",
       x = "",
       y = "asinh Expression")
ggsave("Each_Cluster_Violin_Supplemental.tiff", violins, width = 10, height = 10, dpi = 300)

# Show heatmap by group
groupMean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),
                                  ids = spe$groups,
                                  statistics = "mean",
                                  use.assay.type = "exprs",
                                  subset.row = rownames(spe)[rowData(spe)$use_channel])

# Max expresison scaled heatmap
dittoHeatmap(groupMean,
             assay = "exprs",
             cluster_cols = TRUE,
             scaled.to.max = TRUE,
             heatmap.colors.max.scaled = viridis(100),
             annot.by = c("groups", "ncells"),
             annotation_colors = list(ncells = viridis(100)))

# Display sample frequency per cell type 
dittoBarPlot(spe, 
             var = "roi_identifier_from_dr_moussa", 
             group.by = "Phenotypes")


cell_density <- colData(spe) %>%
  as.data.frame() %>%
  group_by(groups, roi_identifier_from_dr_moussa) %>%
  # Compute the number of pixels covered by cells and 
  # the total number of pixels
  summarize(cell_area = sum(area),
            no_pixels = mean(width_px) * mean(height_px)) %>%
  # Divide the total number of pixels 
  # by the number of pixels covered by cells
  mutate(covered_area = 100*(cell_area / no_pixels))


group_colors <- brewer.pal(10, "Set3")
group_colors <- c("Low Normal" = "#8DD3C7",
                  "Low Polyp" = "#BEBADA",
                  "High Normal" = "#FB8072",
                  "High Polyp" = "#80B1D3")



# Visualize the image area covered by cells per image
test <- ggplot(cell_density, aes(x = groups, y = covered_area, fill = groups)) +
  geom_boxplot(position = "dodge2", outliers = FALSE, varwidth = FALSE, ) + 
  geom_point(position = "jitter") +
  geom_text_repel(aes(label = roi_identifier_from_dr_moussa)) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15)) +
  ylab("% covered area") + xlab("") +
  scale_fill_manual(values = c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")) +
  facet_grid(~groups, scales = "free") +
  theme(legend.position = "none")
ggsave("test.tiff", test, width = 10, height = 10, dpi = 300)

# Read in the non-normalized images
nonNorm_img <- readRDS("../Pratik/study_combined_images.rds")

# Set image names across images/masks and spe
mcols(nonNorm_img) <- mcols(masks) <- DataFrame(sample_id = names(nonNorm_img),
                                                ROI = paste("ROI", unique(spe$roi_identifier_from_dr_moussa), sep = " "))



# Sankey Diagram
sankDat <- read.csv("Walnut_ROI_Sankey.csv", header = TRUE, sep = ",")
sankDat$Group <- gsub("Urolithin High", "Urolithin\nHigh", sankDat$Group)
sankDat$Group <- gsub("Urolithin Low", "Urolithin\nLow", sankDat$Group)
sankDat$Pathology <- gsub("Adj normal", "Adj\nnormal", sankDat$Pathology)


library(alluvial)
tiff("Sankey.tiff", width = 6, height = 6, units = "in", res = 300)
x <- alluvial(sankDat[,1:3], freq = sankDat$IDtarget, col = ifelse(sankDat$Group == "Urolithin\nLow", "#8DD3C7", "#FB8072"), cex = 0.7)
dev.off()





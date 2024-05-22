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
library(patchwork)
library(scater)
library(scales)
library(BiocParallel)

# Read in the modified spe object and normalized images
spe <- readRDS("DATA/Revised_Rds_Objects/April_2024_Edited_by_Mike_Walnut_Spe.Rds")

#----------ADDITIONAL MODIFICATIONS----------#


# Re-classify the clusters (4/11/24)
cluster_celltype <- recode(spe$Phenotypes,
                           "Slow (middle-bottom) epithelia" = "Non-proliferating epithelia",
                           "TA-like/middle epithelia" = "Proliferating epithelia",
                           "Antigen presenting cells" = "Antigen presenting cells",
                           "Slow middle epithelia" = "Non-proliferating epithelia", 
                           "Various immune" = "Lymphocytes including T and B-cells",
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
sum(spe$Phenotypes == "Lymphocytes including T and B-cells") #36145
sum(spe$Phenotypes == "Muscle layer") #25934
sum(spe$Phenotypes == "Stroma") #37755
sum(spe$Phenotypes == "Myeloid lineage") #19586

# Create new rowData entry that is the cluster phenotype with the count of that cluster
cluster_cellNums <- recode(spe$Phenotypes,
                           "Non-proliferating epithelia" = "Non-proliferating epithelia (151581)",
                           "Proliferating epithelia" = "Proliferating epithelia (37895)",
                           "Antigen presenting cells" = "Antigen presenting cells (70295)",
                           "Lymphocytes including T and B-cells" = "Lymphocytes including T and B-cells (36145)",
                           "Muscle layer" = "Muscle layer (25934)",
                           "Stroma" = "Stroma (37755)",
                           "Myeloid lineage" = "Myeloid lineage (19586)")
spe$Phenotype <- cluster_cellNums

# Factor the Phenotype
spe$Phenotype <- factor(spe$Phenotype, levels = c("Non-proliferating epithelia (151581)",
                                                  "Proliferating epithelia (37895)",
                                                  "Lymphocytes including T and B-cells (36145)",
                                                  "Muscle layer (25934)",
                                                  "Stroma (37755)",
                                                  "Antigen presenting cells (70295)",
                                                  "Myeloid lineage (19586)"))

# Clean up the ROI information and factor
spe$roi_identifier_from_dr_moussa <- paste("ROI", spe$roi_identifier_from_dr_moussa, sep = " ")
spe$roi_identifier_from_dr_moussa <- factor(spe$roi_identifier_from_dr_moussa, levels = c(
  "ROI 1", "ROI 2", "ROI 3", "ROI 4", "ROI 5", "ROI 6", "ROI 7", "ROI 8", "ROI 9", "ROI 10",
  "ROI 11", "ROI 12", "ROI 13", "ROI 14", "ROI 15", "ROI 16", "ROI 17", "ROI 18", "ROI 19", "ROI 20",
  "ROI 21", "ROI 22", "ROI 23"
))

# Create a grouping variable
spe$groups <- ifelse(spe$groups == "normal_high_producer", "High Normal",
                     ifelse(spe$groups == "normal_low_producer", "Low Normal",
                            ifelse(spe$groups == "tumor_high_producer", "High Polyp", "Low Polyp")))
spe$groups <- ifelse(spe$roi_identifier_from_dr_moussa == "ROI 7", "Low Normal",
                     ifelse(spe$groups == "Low Normal", "Low Normal", 
                            ifelse(spe$groups == "Low Polyp", "Low Polyp",
                                   ifelse(spe$groups == "High Normal", "High Normal", "High Polyp"))))

# Factor the grouping variable
spe$groups <- factor(spe$groups, levels = c("Low Normal", "Low Polyp", "High Normal", "High Polyp"))


# Assign Type
spe$Type <- ifelse(spe$Phenotypes == "Proliferating epithelia"|spe$Phenotypes == "Non-proliferating epithelia", "Epithelial clusters", "Non-epithelial clusters")
spe$Subgroup <- ifelse(spe$Phenotype == "Proliferating epithelia (37895)"|spe$Phenotype == "Non-proliferating epithelia (151581)", "Epithelial clusters", "Non-epithelial clusters")

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


#---------- Figure 6a

# Calculate the mean expression of each Phenotype
celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$Phenotypes, 
                                      statistics = "mean",
                                      use.assay.type = "exprs")

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
ggsave("FINALIZED_FIGURES/FIGURE_6/Figure_6a.tiff", New, width = 8, height = 8, dpi = 300)

#########################
# INDIVIDUAL ROI HEATMAP
#########################

# Calculate the mean counts for each ROI
image_mean <- aggregateAcrossCells(spe[rowData(spe)$use_channel], 
                                   ids = spe$roi_identifier_from_dr_moussa, 
                                   statistics="mean",
                                   use.assay.type = "counts")

# Transform the counts to asinh expression
assay(image_mean, "exprs") <- asinh(counts(image_mean))

group_coloring <- c("Low Normal" = "#8DD3C7",
                    "Low Polyp" = "#BEBADA",
                    "High Normal" = "#FB8072",
                    "High Polyp" = "#80B1D3")


# Plot heatmap
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

#---------- Figure 6c

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
       title = "",
       x = "") +
  theme(strip.text = element_text(face = "bold", size = 22),
        axis.title.y = element_text(size = 34, face = "bold"),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 20))
ggsave("FINALIZED_FIGURES/FIGURE_6/Figure_6c.tiff", barplot, width = 10, height = 10, dpi = 300)

#----------Figure 6b

# Save SPE inCATALYST-compatible object with renamed colData entries and new metadata information
spe_cat <- spe
spe_cat$sample_id <- spe$roi_identifier_from_dr_moussa
spe_cat$cluster_id <- factor(spe$Phenotype)
spe_cat$condition <- spe$groups

# Add celltype information to metadata
metadata(spe_cat)$cluster_codes <- data.frame(celltype = factor(spe_cat$cluster_id))

# CLR on cluster proportions across samples
PCA <- clrDR(spe_cat, 
             dr = "PCA", 
             by = "sample_id", 
             k = "celltype", 
             label_by = "sample_id", 
             arrow_col = "cluster_id",
             arrow_opa = c(1)) +
  scale_color_manual(values = metadata(spe)$colors$PhenotypeColors) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        axis.title.y = element_text(size = 26, face = "bold"),
        axis.title.x = element_text(size = 26, face = "bold"),
        legend.text = element_text(size = 12)) +
  labs(color = "Phenotype",
       size = "n cells",
       fill = "Group") +
  scale_size(range = c(2, 10))
ggsave("FINALIZED_FIGURES/FIGURE_6/Figure_6b.tiff", PCA, width = 12, height = 12, dpi = 300)

#########################
# SUPPLEMENTAL VIOLIN
#########################

# Sample down to 1000 cells
cur_cells <- sample(seq_len(ncol(spe)), 1000)

# Plot violin plots for each marker and ROI
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

#########################
# SANKEY FLOW DIAGRAM
#########################

# Sankey Diagram
sankDat <- read.csv("Walnut_ROI_Sankey.csv", header = TRUE, sep = ",")
sankDat$Group <- gsub("Urolithin High", "Urolithin\nHigh", sankDat$Group)
sankDat$Group <- gsub("Urolithin Low", "Urolithin\nLow", sankDat$Group)
sankDat$Pathology <- gsub("Adj normal", "Adj\nnormal", sankDat$Pathology)


library(alluvial)
tiff("New_Sankey.tiff", width = 6, height = 6, units = "in", res = 300)
x <- alluvial(sankDat[,1:3], freq = sankDat$IDtarget, col = ifelse(sankDat$Group == "Urolithin\nLow", "#8DD3C7", "#FB8072"), blocks = TRUE, cex = 0.8)
dev.off()

#########################
# PLOT CELLS AND PLOT PIXELS
#########################

# Read in the images and masks
img <- readRDS("/users/michaelmartinez/Desktop/Walnut_IMC_2.0/Walnut_IMC_Manuscript_Revisions_Outputs/Revised_Rds_Objects/April_2024_Edited_by_Mike_asinh_Images.Rds")
masks <- readRDS("/users/michaelmartinez/Desktop/Walnut_IMC_2.0/Walnut_IMC_Manuscript_Revisions_Outputs/Revised_Rds_Objects/April_2024_Edited_by_Mike_Walnut_Masks.Rds")


# Ensure that channelNames of the images match with the rownames of the spatial experiment object
channelNames(img) <- rownames(spe)

# Set image names across images/masks and spe
mcols(img) <- mcols(masks) <- DataFrame(sample_id = names(img),
                                        ROI = unique(spe$roi_identifier_from_dr_moussa))

# Visualize clusters on the tissues
tissues <- names(img)

#########################
# CLUSTER PROJECTIONS (ALL)
#########################

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

#########################
# CLUSTER PROJECTIONS (SINGLE)
#########################

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

#########################
# PLOT CELLS  IMC MARKERS
#########################

# ROI 11 and ROI 22
cur_img <- img[c(9,20)]
cur_mask <- masks[c(9,20)]

# Within a for loop, visualize each images with the ROI as the title
for(i in tissues) {
  cur_img <- img[i]
  cur_mask <- masks[i]
  
  plotCells(cur_mask,
            object = spe,
            cell_id = "ObjectNumber",
            img_id = "sample_id",
            colour_by = c("Vimentin", "CD90", "MHC_II", "Epcam+ECadherin"),
            colour = list(Vimentin = c("black", "red"),
                          CD90 = c("black", "burlywood1"),
                          MHC_II = c("black", "cyan2"),
                          `Epcam+ECadherin` = c("black", "white")),
            exprs_values = "exprs",
            image_title = list(text = mcols(cur_img)$ROI,
                               position = "topright",
                               font = 0),
            legend = list(colour_by.title.cex = 5,
                          margin = 0.1),
            save_plot = list(filename = paste(mcols(cur_img)$ROI, ".tiff", sep = "_")))
  
}

#########################
# SAVE SPE OBJECT
#########################

saveRDS(spe, file = "FULLY_REVISED_Walnut_Spe_Object_4.16.24.Rds")


# Correlation figures
low <- as.data.frame(assay(image_mean, "exprs"))[,1:10]

low_corr <- rcorr(t(as.matrix(low)), type = "spearman")

# Extract correlation coefficients and p-values
cor_matrix <- as.matrix(low_corr$r)
p_values <- as.matrix(low_corr$P)
diag(p_values) <- 1

tiff("Low_IMC_Correlation_Plot.tiff", width = 6, height = 6, units = "in", res = 300)
PostUroCorr <- corrplot(cor_matrix, method = "circle", is.corr = TRUE, tl.cex = 0.5, order = "AOE",
                        addCoefasPercent = FALSE, #addCoef.col = TRUE, number.cex = 0.5, pch.cex = 0.5,
                        type = c("upper"), col = colorRampPalette(c("blue", "white", "red"))(500),
                        diag = TRUE, tl.col = "black", tl.srt = 45)
dev.off()


# Correlation figures
high <- as.data.frame(assay(image_mean, "exprs"))[,11:23]

high_corr <- rcorr(t(as.matrix(high)), type = "spearman")

# Extract correlation coefficients and p-values
cor_matrix <- as.matrix(high_corr$r)
p_values <- as.matrix(high_corr$P)
diag(p_values) <- 1

tiff("High_IMC_Correlation_Plot.tiff", width = 6, height = 6, units = "in", res = 300)
PostUroCorr <- corrplot(cor_matrix, method = "circle", is.corr = TRUE, tl.cex = 0.5, order = "AOE",
                        addCoefasPercent = FALSE, #addCoef.col = TRUE, number.cex = 0.5, pch.cex = 0.5,
                        type = c("upper"), col = colorRampPalette(c("blue", "white", "red"))(500),
                        diag = TRUE, tl.col = "black", tl.srt = 45)
dev.off()

#########################
# SPATIAL ANALYSIS
#########################
spe <- readRDS("FULLY_REVISED_Walnut_Spe_Object_4.16.24.Rds")
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "knn", k = 20)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "delaunay", max_dist = 20)

# Group cells based on information contained in their direct neighborhood.
# Perform Delaunay triangulation-based graph construction, neighborhood aggregation and then cluster cells.
# OR, construct a 10-nearest neighbor graph before aggregating informationa cross neighboring cells.

# For each cell, the function computes the fraction of cells of a certain type (e.g. cell type) among its neighbors
# For each cell it aggregates (e.g., mean) the expression counts across all neighboring cells
# Based on these measures, cells can now be clustered into cellular neighborhoods.

#########################
# BY EXPRESSION
#########################
set.seed(220705)

# Create L-curves, then cluster
spe <- lisaClust(spe, 
                 k = 6,
                 Rs = c(10, 20, 50),
                 spatialCoords = c("Pos_X", "Pos_Y"),
                 cellType = "Phenotype",
                 imageID = "sample_id")

# Plot the region map
regionMap(spe, 
          cellType = "region",
          region = "Phenotype")

# Visualize on the tissue
plotSpatial(spe[,spe$sample_id == "walnut_tumor_high_producer_03_07_2023_mcd3_roi1"],
            img_id = "cn_expression",
            node_color_by = "region",
            node_size_fix = 0.5) +
  scale_color_brewer(palette = "Set3")

# Plot composition as a heatmap
for_plot <- prop.table(table(spe$region, spe$Phenotype), 
                       margin = 1)
pheatmap(for_plot, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")

#########################
# SPATIAL CONTEXT ANALYSIS
#########################

# While CNs can represent sites of unique local processes, SCs describe tissue regions in which
# distinct CNs may be interacting. Hence, SCs may be interesting regions of specialized biological events.

# CN fractions are sorted from high-to-low and the SC of each cell is assigned as the minimal
# combination of SCs that additively surpass a user-defined threshold. 
# The default threshold of 0.9 aims to represent the dominant CNs, hence the most prevalent signals
# in a given window.

# Construct a 40 nearest neighbor graph
spe <- buildSpatialGraph(spe,
                         img_id = "sample_id",
                         type = "knn",
                         name = "knn_spatialcontext_graph",
                         k = 40)

# Compute the fraction of cellular neighborhoods around each cell
spe <- aggregateNeighbors(spe, 
                          colPairName = "knn_spatialcontext_graph",
                          aggregate_by = "metadata",
                          count_by = "region",
                          name = "aggregatedNeighborhood")

# Detect spatial contexts
spe <- detectSpatialContext(spe, 
                            entry = "aggregatedNeighborhood",
                            threshold = 0.90,
                            name = "spatial_context")

# Define SC color scheme
n_SCs <- length(unique(spe$spatial_context))
col_SC <- setNames(colorRampPalette(brewer.pal(9, "Paired"))(n_SCs), 
                   sort(unique(spe$spatial_context)))


# Visualize spatial contexts on images
plotSpatial(spe[,spe$sample_id == "walnut_tumor_high_producer_03_07_2023_mcd3_roi1"], 
            node_color_by = "spatial_context", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_manual(values = col_SC)

# Compare CN and SC for one patient 
p1 <- plotSpatial(spe[,spe$sample_id == "walnut_tumor_high_producer_03_07_2023_mcd3_roi1"], 
                  node_color_by = "region", 
                  img_id = "sample_id", 
                  node_size_fix = 0.5) +
  scale_color_brewer(palette = "Set3")

p2 <- plotSpatial(spe[,spe$sample_id == "walnut_tumor_high_producer_03_07_2023_mcd3_roi1"], 
                  node_color_by = "spatial_context", 
                  img_id = "sample_id", 
                  node_size_fix = 0.5) +
  scale_color_manual(values = col_SC, limits = force)
p1 + p2

## Filter spatial contexts
# By number of group entries
spe <- filterSpatialContext(spe, 
                            entry = "spatial_context",
                            group_by = "sample_id", 
                            group_threshold = 3,
                            name = "spatial_context_filtered")

plotSpatial(spe, 
            node_color_by = "spatial_context_filtered", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_manual(values = col_SC, limits = force)


# Filter out small and infrequent spatial contexts
spe <- filterSpatialContext(spe, 
                            entry = "spatial_context",
                            group_by = "sample_id", 
                            group_threshold = 3,
                            cells_threshold = 100,
                            name = "spatial_context_filtered")

plotSpatial(spe, 
            node_color_by = "spatial_context_filtered", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_manual(values = col_SC, limits = force)


## Plot spatial context graph 

# Colored by name, size by n_cells
plotSpatialContext(spe, 
                   entry = "spatial_context_filtered",
                   group_by = "sample_id",
                   node_color_by = "name",
                   node_size_by = "n_cells",
                   node_label_color_by = "name")

plotSpatialContext(spe, 
                   entry = "spatial_context_filtered",
                   group_by = "sample_id",
                   node_color_by = "n_cells",
                   node_size_by = "n_group",
                   node_label_color_by = "n_cells") +
  scale_color_viridis()

#########################
# INTERACTION ANALYSIS
#########################

# This section focuses on statistically testing the pairwise interaction between all
# phenotypes of the dataset. The `testInteractions` function implements interaction testing.

# Per grouping level, the function computes the averaged cell type/cell type interaction count
# and compares this count against an empirical null distribution which is generated by 
# permuting all cell labels (while maintaining in the tissue structure).
library(scales)
out2 <- testInteractions(spe, 
                         group_by = "groups",
                         label = "Phenotype", 
                         colPairName = "neighborhood",
                         BPPARAM = SerialParam(RNGseed = 03061999))

head(out2)

# Test the interactions by summing the interaction/avoidance values across each tissue
out2 %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Read in the differential marker data
diff <- read.csv("DifferentialMarkers.csv", header = TRUE, sep = ",")
diff$Marker <- factor(diff$Marker, levels = diff$Marker)
enrichmentColor <- brewer.pal(3, "Set1")
enrichColors <- c("Decreased in High" = "#E41A1C",
                  "Increased in High" = "#4DAF4A")


diffExpression <- ggplot(diff, aes(x = Fold.Change, y = Marker, fill = Group, label = Pcode)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_text(vjust = ifelse(diff$Group == "Decreased in High", 1.5, 0), size = 9) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A")) +  # Red for negative, blue for positive
  coord_flip() +  # Flip coordinates to have horizontal bars
  theme_classic() +  # Use a minimal theme
  labs(y = "", 
       x = "Log2 Fold Change") +
  facet_grid(~Group, scales = "free") +
  theme(panel.grid.major = element_blank(),  # Remove grid lines for clarity
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 26),
        axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        strip.text = element_text(size = 26, face = "bold"),
        legend.position = "none")
ggsave("Walnut_Differential_Marker_Expression.tiff", diffExpression, width = 10, height = 10, dpi = 300)









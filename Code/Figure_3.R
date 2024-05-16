# Figure 3: Blood serum biomarkers and correlation to urolithin A levels

# Load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(ggpubr)


#----------Figure 3a
# Serum Inflamatory Marker Heatmap
# Read in the serum inflammatory data
serum <- read.csv("/users/michaelmartinez/Desktop/Walnut_Project/SERUM_FIGURES/Serum_Inflammatory_Markers_Raw.csv", header = TRUE, sep = ",")

# Get a vector of sample IDs and Urolithin Group data
Groups <- serum[,2:3]
serum[,2:3] <- NULL
rownames(serum) <- serum$Sample.ID
serum$Sample.ID <- NULL

# Log transform all the values
serum_log1p <- sapply(serum, log1p)
rownames(serum_log1p) <- rownames(serum)
serum_log1p <- as.data.frame(serum_log1p)

# Heatmap values
heatmapNames <- c("c_peptide", "sICAM_1", "sIL_6R", "Ghrelin", "Haptoglobin", 
                  "TRAIL", "sVEGFR2", "MCP_2", "PDGF_AB_BB")

# Subset the serum_log1p dataframe to only include these columns
serum_subset <- serum_log1p[,colnames(serum_log1p) %in% heatmapNames]

# Append the values of the urolihin group
serum_subset$Urolithin <- Groups$URoAClass

# Take just the lows, mediums, and highs as separate dataframes
serum_subset_low <- serum_subset[serum_subset$Urolithin == "1",]
serum_subset_low$Urolithin <- NULL
low_means <- as.data.frame(colMeans(serum_subset_low))
colnames(low_means) <- c("Low")

serum_subset_med <- serum_subset[serum_subset$Urolithin == "2",]
serum_subset_med$Urolithin <- NULL
med_means <- as.data.frame(colMeans(serum_subset_med))
colnames(med_means) <- c("Med")

serum_subset_high <- serum_subset[serum_subset$Urolithin == "3",]
serum_subset_high$Urolithin <- NULL
high_means <- as.data.frame(colMeans(serum_subset_high))
colnames(high_means) <- c("High")

# Combine low, med, high into one dataframe
log1p_means <- cbind(cbind(low_means, med_means),high_means)

# Row order
markerOrder <- heatmapNames
log1p_means <- log1p_means[markerOrder,]

# Complex Heatmap
# Create a named vector for the colors
sample_colors <- c("#FF0000BF", "#A6CEE3", "#FF8000BF")
names(sample_colors) <- c("Low", "Med", "High")

#Set heatmap splitting pattern
hmSplit <- rep(1:3, c(1,1,1))

#Define the number of slices in the heatmap
slices <- 1+1+1

#Create a heatmap annotation
anno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#FF0000BF", "#A6CEE3", "#FF8000BF")), 
                     labels = c("Low", "Med", "High"),
                     labels_gp = gpar(col = "black", fontsize = 26, fontface = 2)),
  col = list(Group = sample_colors),
  show_annotation_name = FALSE
)


zscores <- t(apply(log1p_means, 1, scale))
rownames(zscores) <- c("C-peptide", "sICAM-1", "sIL-6R",
                       "Ghrelin", "Haptoglobin", "TRAIL",
                       "sVEGFR2", "MCP-2", "PDGF-AB/BB")

# Plot the heatmap
scaled <- Heatmap(zscores,
                  column_labels = colnames(zscores), 
                  show_column_names = FALSE,
                  name = " ",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  top_annotation = anno,
                  column_split = hmSplit,
                  column_title = NULL,
                  row_names_side = "left",
                  row_title = NULL,
                  col = viridis(100),
                  row_names_gp = gpar(fontsize = 16))

tiff("FINALIZED_FIGURES/FIGURE_3/Figure_3A_Walnuts_Log1p_colMeans_Heatmap.tiff", width = 8, height = 8, units = "in", res = 300)
draw(scaled)
dev.off()

#----------Figure 3c
# Serum inflammatory correaltion plots
# Read in the urolithin vaues
uroVals <- read.csv("/users/michaelmartinez/Desktop/Walnut_Project/SERUM_FIGURES/UroA_Delta.csv", header = TRUE, sep = ",")
uroVals <- uroVals[,c(1,3,4,7:8,11:ncol(uroVals))]

# Keep these columns
keep <- c("Sample", "ClassDelta", "Delta", "Sample.ID", "c_peptide", "sICAM_1", "EGF", "Ghrelin", "GRO_alpha", 
          "TRAIL", "Insulin", "PYY", "PDGF_AB_BB")
uroVals <- uroVals[,colnames(uroVals) %in% keep]
uroVals$Sample.ID <- NULL


# Get metadata
meta <- uroVals[,c(1:2)]
uroVals[,1:2] <- NULL
colnames(uroVals) <- c("Delta", "EGF" , "PYY", "sICAM-1", "PDGF-AB/BB" , "GRO-alpha", "TRAIL", "C-peptide", "Ghrelin", "Insulin")

# Get rid of negatives
#uroVals$Delta <- uroVals$Delta + 144.54

# Log transform the values
uroValsLog <- apply(uroVals, 2, log1p)
uroValsLog[is.nan(uroValsLog)] <- 0

# Set rownames
rownames(uroValsLog) <- meta$Sample
uroValsLog <- as.data.frame(uroValsLog)
uroValsLog$Group <- meta$ClassDelta
uroValsLog$Group <- ifelse(uroValsLog$Group == "1", "Low",
                           ifelse(uroValsLog$Group == "2", "Med", "High"))

# Pivot longer
uroValsLogLong <- pivot_longer(uroValsLog, -c("Delta", "Group"), names_to = "marker", values_to = "values")
uroValsLogLong$Group <- factor(uroValsLogLong$Group, levels = c("Low", "Med", "High"))
uroValsLogLong <- uroValsLogLong[uroValsLogLong$marker != "sICAM-1",]
uroValsLogLong$marker <- factor(uroValsLogLong$marker, levels = c("C-peptide", "Ghrelin", "EGF", "GRO-alpha", 
                                                                  "PDGF-AB/BB", "TRAIL","Insulin", "PYY"))


# Plot
serumCorr <- ggplot(uroValsLogLong, aes(x = Delta, y = values, color = Group)) +
  geom_point(size = 4) +
  facet_wrap(~marker, scales = "free_y", nrow = 2) +
  stat_cor(inherit.aes = FALSE, aes(x = Delta, y = values), method = "spearman", size = 8) +
  geom_smooth(inherit.aes = FALSE, aes(x = Delta, y = values), method = "lm", se = TRUE) +
  scale_color_manual(values = c("Low" = "#FF0000BF",
                                "Med" = "#A6CEE3",
                                "High" = "#FF8000BF")) +
  theme_classic() +
  labs(y = "Log serum marker levels (pg/ml)",
       x = "Log creatinine-normalized delta urolithin A levels (ng/mg)") +
  theme(strip.text = element_text(size = 32, face = "bold"),
        legend.position = "none",
        axis.title.x = element_text(size = 46, face = "bold"),
        axis.title.y = element_text(size = 46, face = "bold"),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32))

ggsave("FINALIZED_FIGURES/FIGURE_3/Final_Figure_3C.tiff", serumCorr, width = 19, height = 14, dpi = 300)




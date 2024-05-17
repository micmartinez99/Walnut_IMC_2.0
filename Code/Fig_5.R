# Figure 5: IMC analysis of polyps

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)


#---------- Figure 5b
# Differential Expression Barplot

# Read in the differential marker data
diff <- read.csv("DATA/Differential_Expression/DifferentialMarkers.csv", header = TRUE, sep = ",")

# Factor the markers so they appear in the order they are in the file, not alphabetically
diff$Marker <- factor(diff$Marker, levels = diff$Marker)

# Set a color scheme
enrichColors <- c("Decreased in High" = "#E41A1C",
                  "Increased in High" = "#4DAF4A")

# Plot
diffExpression <- ggplot(diff, aes(x = Fold.Change, y = Marker, fill = Group, label = Pcode)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_text(vjust = ifelse(diff$Group == "Decreased in High", 1.5, 0), size = 9) +
  scale_fill_manual(values = c("#E41A1C", "#4DAF4A")) +  # Red for negative, blue for positive
  coord_flip() +  # Flip coordinates to have horizontal bars
  theme_classic() +  # Use a minimal theme
  labs(y = "", 
       x = "Log2 fold change") +
  facet_grid(~Group, scales = "free") +
  theme(panel.grid.major = element_blank(),  # Remove grid lines for clarity
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 32),
        axis.title.y = element_text(size = 36, face = "bold"),
        axis.text.x = element_text(size = 24, face = "bold"),
        axis.title.x = element_text(size = 24),
        strip.text = element_text(size = 26, face = "bold"),
        legend.position = "none")
ggsave("FINALIZED_FIGURES/FIGURE_5/Figure_5b.tiff", diffExpression, width = 10, height = 10, dpi = 300)

#----------Figure 5c
# Vimentin Expression Violin Plots
# Vimentin violin

# Read in the asinh-transformed vimentin levels for each cell in the dataset
vimentin <- read.csv("DATA/IMC_Expression_Values/asinh_expression_levels_Per_ROI_Per_Marker.csv")

# Factor the ROIs
vimentin$ROIs <- factor(vimentin$ROI, levels = c(
  "ROI 1", "ROI 2", "ROI 3", "ROI 4", "ROI 5", "ROI 6", "ROI 7", "ROI 8", "ROI 9", "ROI 10",
  "ROI 11", "ROI 12", "ROI 13", "ROI 14", "ROI 15", "ROI 16", "ROI 17", "ROI 18", "ROI 19", "ROI 20",
  "ROI 21", "ROI 22", "ROI 23"))

# Assign the ROIs their urolithin group
vimentin$Group <- ifelse(vimentin$ROIs == "ROI 1" | 
                           vimentin$ROIs == "ROI 2" |
                           vimentin$ROIs == "ROI 3" |
                           vimentin$ROIs == "ROI 4" |
                           vimentin$ROIs == "ROI 5" |
                           vimentin$ROIs == "ROI 6" |
                           vimentin$ROIs == "ROI 7" |
                           vimentin$ROIs == "ROI 8" |
                           vimentin$ROIs == "ROI 9" |
                           vimentin$ROIs == "ROI 10", "Urolithin Low", "Urolithin High")

# Factor the group
vimentin$Group <- factor(vimentin$Group,
                         levels = c("Urolithin Low", "Urolithin High"))

# Set group colors
vimColors <- c("Urolithin Low" = "#FF0000BF",
               "Urolithin High" = "#FF8000BF")

#Plot width scaled violin plot
vimExpr <- ggplot(vimentin, aes(x = ROIs, y = as.numeric(Vimentin), fill = Group)) +
  theme_classic() +
  labs(x = "",
       y = "Vimentin expression") +
  theme(axis.text.x = element_text(size = 30, angle = 90),
        axis.text.y = element_text(size = 32),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 34, face = "bold"),
        strip.text = element_text(size = 32, face = "bold"),
        legend.position = "none") +
  geom_violin(trim = FALSE, scale = "width", adjust = 1, na.rm = TRUE, draw_quantiles = TRUE) +
  scale_fill_manual(values = vimColors) +
  facet_wrap(~Group, scales = "free_x") +
  geom_boxplot(outlier.shape = 19, outlier.size = 0.1, na.rm = TRUE, varwidth = TRUE, width=0.05)
ggsave("FINALIZED_FIGURES/FIGURE_5/Figure_5c_ainh.tiff", vimExpr, width = 10, height = 10, dpi = 300)


#----------Figure 5d
# Fecal Urolithin A Correlations with IMC markers
# Read in the data
data <- read.csv("DATA/IMC_Expression_Values/Mean_Markers_Per_ROIs.csv")
data[,2:6] <- NULL
rownames(data) <- data$ROI
data$ROI <- NULL

# Get the fecal data in a separate dataframe
fecal <- data[,1:4]
data[,1:4] <- NULL

# Take the rows that correspond to vimentin and CD163
vimentinCD163 <- data[,c(2,4)]

# Add in the fecal data and uroGroup information
vimentinCD163$FecalD <- fecal$Fecal_D
vimentinCD163$UroGroup <- fecal$UroA_Class

# Rename groups to be consistent with the rest of the paper and factor so urolithin low comes before urolithin high
vimentinCD163$UroGroup <- ifelse(vimentinCD163$UroGroup == "Low", "Urolithin Low", "Urolithin High")
vimentinCD163$UroGroup <- factor(vimentinCD163$UroGroup, levels = c("Urolithin Low", "Urolithin High"))
vimentinCD163$ROI <- rownames(vimentinCD163)
fecal$ROI <- rownames(fecal)

# Pivot longer for plotting
vim <- vimentinCD163 %>%
  pivot_longer(-c(ROI, UroGroup, FecalD), names_to = c("Marker"), values_to = "value") 
fecal[,c(2,3)] <- NULL

# Factor so vimentin comes before CD163 when plotting
vim$Marker <- factor(vim$Marker, levels = c("Vimentin", "CD163"))

# Plot a facetted plot
corplot <- ggplot(vim, aes(x = log(FecalD), y = value, color = UroGroup, label = ROI)) +
  geom_point(size = 4) +
  stat_cor(inherit.aes = FALSE, aes(x = log(FecalD), y = value), method = "spearman", size = 8, vjust = -3) +
  geom_smooth(inherit.aes = FALSE, aes(x = log(FecalD), y = value), method = "lm", se = TRUE) +
  theme_classic() +
  facet_wrap(~Marker, scales = "free") +
  geom_text_repel(aes(label = ROI)) +
  scale_color_manual(values = c("Urolithin Low" = "#FF0000BF", "Urolithin High" = "#FF8000BF")) +
  labs(y = "Mean expression",
       x = "Log10 Delta Fecal Urolithin A Levels (ng/ml)",
       color = "") +
  theme(strip.text = element_text(size = 32, face = "bold"),
        axis.title.x = element_text(size = 34, face = "bold"),
        axis.title.y = element_text(size = 40, face = "bold"),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32),
        legend.position = "bottom",
        legend.text = element_text(size = 32))
ggsave("FINALIZED_FIGURES/FIGURE_5/Figure_5d.tiff", corplot, width = 12, height = 12, dpi = 300)





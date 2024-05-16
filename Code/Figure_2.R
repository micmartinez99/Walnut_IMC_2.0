###### Figure 2: Urolithin Levels Presented as model-based clustering and correlation analysis

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(corrplot)
library(Hmisc)



#---------- Figure 2b
# Urolithin metabolites as paired boxplots

# Read in the data
uro_data <- read.csv("/users/michaelmartinez/Desktop/Walnut_Project/Code_RawData_for_Figure_Generation/Urolithin_Levels/Raw_Urolithin_Data/Urolithin_LongFormat.csv", header = TRUE, sep = ",")

# Change timepoint to "Pre" and "Post" instead of "Before" and "After"
uro_data$Timepoint <- ifelse(uro_data$Timepoint == "Before", "Pre", "Post")

# Factor the timepoint order so "Pre" comes before "Post"
time_order <- c("Pre",
                "Post")
uro_data <- uro_data %>%
  mutate(Timepoint = factor(Timepoint, levels = time_order))

# Factor the urolithin stratification order so it goes "Low", "Med", "High"
group_order <- c("Low", "Med", "High")
uro_data<- uro_data %>%
  mutate(Group = factor(Group, levels =  group_order))

# Sample 2B was missing Urolithin M7...value is 0

# Replace periods with spaces
uro_data$Metabolite <- gsub("\\.", " ", uro_data$Metabolite)

# Plot the paired plots
urolithins_paired <- ggpaired(uro_data, x = "Timepoint", y = "ng.mg",
                              fill = "Timepoint", line.color = "black", line.size = 0.1) +
  labs(y = "Raw urolithin levels (ng/mg)",
       x = "") +
  theme(legend.position = "right") +
  scale_fill_manual(values = c("Pre" = "#B2DF8A", "Post" = "#FB9A99")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 34, face = "bold"),
        axis.title.y = element_text(size = 40, face = "bold"),
        axis.text.y = element_text(size = 32),
        strip.text = element_text(size = 28, face = "bold")) +
  facet_wrap(~Metabolite, scales = "free_y") +
  stat_compare_means(method = "wilcox", paired = TRUE, label = "p.signif", vjust = 0.5, size = 8.5)
urolithins_paired
ggsave("FINALIZED_FIGURES/FIGURE_2/Final_Figure_2b.tiff", urolithins_paired, width = 14, height = 14, dpi = 300)


#---------- Figure 2c
# Correlation plot
# Read in the table that has creatinine-normalized urolithin values
vals <- read.csv("/users/michaelmartinez/Desktop/Walnut_Project/Code_RawData_for_Figure_Generation/Urolithin_Levels/Raw_Urolithin_Data/Urolithin 4.13.22_results_v1.xlsx - Urolithin 4.13.22_prep.csv")

# Take just the post columns
after <- vals[,c(1,14:ncol(vals))]
rownames(after) <- after$Sample
after$Sample <- NULL
colnames(after) <- c("Post isourolithin A", "Post urolithin A", "Post urolithin B",
                     "Post urolithin C", "Post urolithin D", "Post urolithin E",
                     "Post urolithin M5", "Post urolithin M6", "Post Urolitin M7")


# Generate the correlation matrix using Spearman
correlationMatrix2 <- rcorr(as.matrix(after), type = "spearman")

# Extract correlation coefficients and p-values
cor_matrix <- as.matrix(correlationMatrix2$r)
p_values <- as.matrix(correlationMatrix2$P)
diag(p_values) <- 1

# Save a tiff file
tiff("FINALIZED_FIGURES/FIGURE_2/Final_Figure_2c.tiff", width = 6, height = 6, units = "in", res = 300)
PostUroCorr <- corrplot(cor_matrix, method = "color", is.corr = TRUE,
                        addCoefasPercent = FALSE, addCoef.col = TRUE, number.cex = 0.5, pch.cex = 0.5,
                        type = c("lower"), col = colorRampPalette(c("blue", "white", "red"))(500),
                        diag = TRUE, tl.col = "black", tl.srt = 45)
dev.off()

write.csv(cor_matrix, file = "Urolithin_post_Spearman_Correlations.csv")
write.csv(p_values, file = "Urolithin_post_Spearman_Pvals.csv")

#---------- Figure 2e
# Urolithin PCA, based on the full panel of pre and post, clustered by delta GMM classification
# Read in the urolithin data
uros <- read.csv("/users/michaelmartinez/Desktop/Walnut_Project/Code_RawData_for_Figure_Generation/Urolithin_Levels/Raw_Urolithin_Data/Urolithin 4.13.22_results_v1.xlsx - Urolithin 4.13.22_prep.csv")

# Take the full panel (pre and post) for all 9 metabolites
uro <- uros[,5:22]

# Run PCA and scale
test <- prcomp(uro, scale. = TRUE)

# Extract scores
pc_scores <- as.data.frame(test$x)

# Append the delta stratification and assign Low, Med, High
pc_scores$UroClassDelta <- uros$ClassDelta
pc_scores$UroClassDelta <- ifelse(pc_scores$UroClassDelta == "C1", "Low",
                                  ifelse(pc_scores$UroClassDelta == "C2", "High", "Med"))

# Factor stratification so it is in the correct order
pc_scores$UroClassDelta <- factor(pc_scores$UroClassDelta, levels = c("Low", "Med", "High"))

# Assign color scheme to levels
pc_scores$color <- ifelse(pc_scores$UroClassDelta == "Low", "#FF0000BF", 
                          ifelse(pc_scores$UroClassDelta == "Med", "#A6CEE3", "#FF8000BF"))


# Plot PCA 
PCA <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = UroClassDelta)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("Low" = "#FF0000BF",
                                "Med" = "#A6CEE3",
                                "High" = "#FF8000BF")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 34),
        axis.title.x = element_text(size = 34, face = "bold"),
        axis.text.y = element_text(size = 32),
        axis.title.y = element_text(size = 34, face = "bold")) +
  labs(color = "GMM-based Urolithin A group")
ggsave("FINALIZED_FIGURES/FIGURE_2/Final_Figure_2e.tiff", PCA, width = 8, height = 8, dpi = 300)



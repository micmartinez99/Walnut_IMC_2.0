# The purpose of this script is to analyze the data from the Aksenov lab

library("ggplot2")
library("ggh4x")
library("ggrepel")
library("tidyverse")
library("dplyr")
library("ggpubr")
library("corrplot")
library("mclust")
library("PCAtools")

data <- read.csv("C18_Neg_Fecal.csv", header = TRUE, sep = ",")
pca <- prcomp(data[,2:3])
pcaRes <- as.data.frame(pca[["x"]])
pcaRes$Sample <- data$Sample
pcaRes$Group <- data$Group

ggplot(pcaRes, aes(x = PC1, y = PC2, color = Group)) +
  geom_point()


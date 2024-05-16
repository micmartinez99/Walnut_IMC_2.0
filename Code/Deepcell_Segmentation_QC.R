# The purpose of this script is to generate 3 Rds objects
# 1: Spatial experiment object (spe) for a combined steinbock output of all Fluidigm samples (Slide2_5 and Slide1_6)
# 2: Normalized images (two step normalization across images and channels)
# 3: Masks
# All three of these objects will have a unified metadata scheme as specified in the code below

setwd("/users/michaelmartinez/Desktop/Walnut_2.0/Deepcell_Segmentation_QC_Outputs/")

library(imcRtools)
library(cytomapper)
library(RColorBrewer)
library(dittoSeq)
library(viridis)
library(gridGraphics)
library(cowplot)

# First, we need to read in the Steinbock folder
spe <- read_steinbock("../Testing/1-1/Steinbock/")

# Check everything looks good so far...should be 26 total
length(unique(spe$sample_id))
unique(spe$sample_id)

# Specify a subset of the channels (we don't care about ICSK1-3 or DNA channels)
rowData(spe)$use_channel <- !grepl("DNA1|DNA3", rownames(spe))

# Transform the counts
assay(spe, "exprs") <- asinh(counts(spe)/1)

# Set channel names
rownames(spe) <- c("aSMA", "Vimentin", "CD66b","CD163","CD31","PDPN", "CD56",
                   "CD20", "CD90", "CD45RA", "CD11c", "CD4", "Epcam_ECad", "MHC_II",
                   "CD8a", "CD45", "CD15", "Foxp3", "Tryptase", "GZMB", "Ki67", "Collagen_I",
                   "CD3", "MHC_I", "CD45RO", "Pan_CK", "DNA1", "DNA3")

# Now, read in the images and deepcell masks
images <- loadImages("../Testing/1-1/Steinbock/img/")
masks <- loadImages("../Testing/1-1/Steinbock/masks/", as.is = TRUE)

# Ensure that channel names for images match the spe object
channelNames(images) <- rownames(spe)

# Check everything looks good again...should be 26 total
length(names(images))
names(images)
length(names(masks))
names(masks)

# Set image names across images/masks and spe
mcols(images) <- mcols(masks) <- DataFrame(sample_id = names(images))

# Normalize images and channels
set.seed(03061999)

# Normalize images first, then channels
cur_images <- cytomapper::normalize(images, separateImages = TRUE, ft = c(0,1))
cur_images2 <- cytomapper::normalize(cur_images, inputRange = c(0,0.1))

# Iterate through the ROIs and get segmentation QC images
ROIs <- names(cur_images2)

for(i in 1:length(ROIs)) {
  print(ROIs[i])
  
  segmentQC <- plotPixels(image = cur_images2[i],
                          #mask = masks[i],
                          object = spe,
                          cell_id = "ObjectNumber",
                          img_id = "sample_id",
                          missing_colour = "white",
                          colour_by = c("DNA1"),
                          colour = list(
                                        DNA1 = c("black", "white")),
                          image_title = list(text = c(ROIs[i]),
                                             position = "topleft",
                                             colour = "white",
                                             margin = c(0,5),
                                             cex = 2),
                          legend = list(colour_by.title.cex = 1,
                                        colour_by.labels.cex = 1),
                          return_plot = TRUE)
  
  QC <- ggdraw(segmentQC$plot, clip = "on")
  ggsave(paste(ROIs[i],"Pan_Cytokeratin_Normalized.tiff", sep = "_"), QC, width = 12, height = 8, dpi = 300)
}

# Save RDS objects
saveRDS(spe, file = "Spe_All_Fluidigm_EOLO_Samples.Rds")
saveRDS(cur_images2, file = "Spe_All_Fluidigm_Normalized_Images.Rds")
saveRDS(masks, file = "Spe_All_Fluidigm_Deepcell_Masks.Rds")


plotPixels(image = images[1],
           colour_by = c("Vimentin"),
           colour = list(Vimentin = c("black", "green")))


















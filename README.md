 # Background & Aims
Diet affects cancer risk and plant-derived polyphenols exhibit antioxidant and cancer-preventive properties. Walnuts are an exceptional source of the polyphenolic ellagitannin, pedunculagin, that is converted into urolithins by gut microflora. The following study examines the impact of urolithin metabolism on inflammatory markers in blood and colon polyp tissue. 

# Methods
In this clinical trial, we have evaluated the effects of walnut consumption on urinary urolithins, serum inflammatory markers and immune cell markers in polyp tissues obtained from 39 non-obese and obese subjects, with and without polyps. Together with detailed food frequency data, we have performed in-depth and integrated computational analysis of metabolomics data combined with serum inflammatory markers and spatial imaging of polyp tissues using imaging mass cytometry (IMC). 

# Results
LC-MS/MS analyses of urine and fecal samples identifies a widely divergent capacity to form nine different urolithin metabolites in this patient population. Subjects with higher urolithin A formation exhibit significantly lower levels of several key serological markers of inflammation, including C-peptide, sICAM 1, sIL6R, Ghrelin, TRAIL, sVEGFR2, PDGF and MCP2, an effect that is more pronounced in obese individuals for several key markers. Furthermore, spatial analysis of IMC data of colon polyp tissues taken from patients shows a significant reduction of vimentin and CD163 expression in subjects with a higher capacity to form urolithin A. 

# Conclusions
These studies indicate that the ability to form urolithin A is linked to general anti-inflammatory responses, warranting further studies to better understand the role of urolithins in the inflammatory process. (Clinicaltrials.gov, Number: NCT04066816).

# !!READ!! Repository Notes
This repository contains the analyses for metabolomics (urine and fecal), serum inflammatory markers, and imaging mass cytometry. (Note: the walnut project also contains microbiome and bulk RNA-Seq components. For microbiome component, please contact Dr. Christian Jobin and Dr. Raad Gharaibeh (University of Florida). For bulk RNA-Seq, analysis was conducted by the Computational Biology Core (UConn / UCHC). These sub-projects are not included in this manuscript or repository. Additionally, Bulk RNA-Seq on PBMCs from patients belonging to this manuscript was also conducted by the Giardina Lab at UConn.)

# Analysis of Metabolomics Data
Urine samples from two timepoints (before and after walnut supplementation) were submitte to Dr. Anthony Provatas (CESE, UConn). Hydrolyzed and unhydrolyzed metabolite levels were quantified for 9 urolitins (isourolithin A, urolithin A, B, C, D, E, M5, M6, and M7.) For preliminary analyses, hydrolyzed values were used by normalizing metabolite levels to creatinine. Samples were stratified into low, medium, or high urolithin producers based on the delta value of urolithin A (after-walnut creatinine-normalized hydrolyzed urolithin A value - before-walnut creatinine-normalized hydrolyzed urolithin A value.) Urolithin A was found to be the metabolite driving the differences between timepoints (as determined by principal component analysis.) Three group stratification was performed using Gaussian mixture modelling with M=3. 
For the analysis of fecal urolithin A levels, samples from the top and bottom tertiles of urinary urolithin A producers were choosen (low and high producers). From these patients, fecal smaples were processed for unhydrolyzed urolitin A analysis in Dr. Alexander Aksenov's Lab (UConn.) These values were used in conjunction with the imaging mass cytometry data to correlate fecal urolithin A levels with protein abundance. An important note: fecal samples were NOT re-stratified. I.e., if a sample came from a urinary low producer, the fecal sample was categorized as a low producer as well. 

# Analysis of Serum Inflammatory Markers
Serum from patients were sent to Eve Technologies for a 73-plex inflammatory marker panel. Data was correlated to the delta levels of urolithin A. Serum inflammatory markers were transformed using the log1p function. For serum inflammatory heatmap, marker levels were averaged across all the samples belonging to each group and then Z-scored for visualization. 







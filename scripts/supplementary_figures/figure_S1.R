# Packages
require(vegan)
require("ggrepel")
require(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
require(ggpubr)
require(ggfortify)
library(stringr) 
require(stats)
library(factoextra)
library(ggsci)
library(hrbrthemes)
library(gridExtra)
library(cowplot)

# parameters

text_size = 12
subtext_size = 10

# Directories

code_dir <- "~/Wasney-Briscoe/scripts/figures/"
data_dir <-  "~/"

source(paste0(code_dir,"helper_scripts/humanized_mouse_utilities.R"))


### reload relab data
mouse_relabs <- read.csv(paste0(data_dir,"merged_data/species/relative_abundance.txt.bz2"),sep="\t",stringsAsFactors = TRUE,row.names = 1)
mouse_relabs <- as.data.frame(t(mouse_relabs))
rownames(mouse_relabs) <- chartr(".", "-", (rownames(mouse_relabs)))

### Helper functions
source(paste0(code_dir,"helper_scripts/mouse_annotation_functions.R"))

### Renaming dataframes as necessary
mouse_metadata <- metadata

# MDS 
# STEP 1. Remove columns with all zeros
mouse_relabs <- mouse_relabs %>%
  select_if(colSums(.) != 0)

# STEP 2. Calculate BC dissimilarity
dissimilarity <- as.matrix(vegdist(mouse_relabs, "bray", na.rm=TRUE))

# STEP 3. MDS scaling
mds_object <- cmdscale(dissimilarity, 2, eig = TRUE, add = TRUE)
mds <- as.data.frame(mds_object$points)
eig <- mds_object$eig
var_explained <- (eig/sum(eig))*100

# Labeling the MDS dataframe
mds$accession <- rownames(mds)
mds <- mds %>%
  rowwise() %>%
  mutate(mouse = extract_mouse(accession, mouse_metadata),
         Diet = extract_diet(accession, mouse_metadata),
         Tissue = extract_tissue_type(accession, mouse_metadata))
mds$mouse <- as.character(mds$mouse)
mds[mds$accession == "TL1gDNAshort", "mouse"] <- "Inoculum"
mds <- mds %>%
  rowwise() %>%
  mutate(Tissue = ifelse(accession == "TL1gDNAshort", "Inoculum", Tissue))

# Plotting
# Factor
mds$tissue <- factor(mds$Tissue, levels = c("Duodenum",  "Jejunum", "Ileum", "Cecum", "Colon", "Inoculum"))
mds$mouse <- factor(mds$mouse , levels = c("1", "2", "3", "4", "5", "6", "7", "8", "Inoculum"))


x_axis_label <- paste0("MDS 1 (", round(var_explained[1], 1), "% variance explained)")
y_axis_label <- paste0("MDS 2 (", round(var_explained[2], 1), "% variance explained)")

figure_S1 <- ggplot(mds, aes(V1, V2, color = Tissue, label = mouse)) +
  # geom_blank() +
  geom_point(size = 3) +
  scale_color_manual(values=c('Duodenum'='#72c8f1',
                              'Jejunum'='#006c37',
                              'Ileum'='#dcc573',
                              'Cecum'='#d05c6f',
                              'Colon'='#312c77',
                              'Inoculum' = '#644117')) +
  geom_text_repel(aes(V1, V2, color = Tissue, label = mouse), size = 4, max.overlaps = Inf) +
  theme_bw() +
  labs(x = x_axis_label, y = y_axis_label) +
  theme(text = element_text(family = "Helvetica", size = text_size),
        # axis.title.x = element_text(),
        # axis.title.y = element_text(),
        axis.text = element_text(size = subtext_size),
        legend.text = element_text(size = subtext_size),
        legend.title = element_blank(),
        legend.key.width = unit(0.01, "npc"),
        legend.key.size = unit(0.8, 'lines'),
  ) +
  guides(shape = guide_legend(override.aes = list(size = 3)), color = guide_legend(override.aes = list(size = 3)))

figure_S1_legend <- as_ggplot(get_legend(figure_S1))

figure_S1 <- figure_S1 + theme(legend.position = "none")

legend_adjustment <- 0.2

figure_S1_final <- ggdraw() +
  draw_plot(figure_S1, 0, 0,1-legend_adjustment,1) +
  draw_plot(figure_S1_legend, 1-legend_adjustment,0,legend_adjustment,1) 

out_path <- "~/figure_S1.png"
ggsave(out_path, figure_S1_final, dpi = 300, bg = "white", height = 5, width = 6, units = "in")






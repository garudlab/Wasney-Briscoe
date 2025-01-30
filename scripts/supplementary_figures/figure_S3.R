# Packages
library(dplyr)
library(lmtest)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(scales)
library(compositions)
library(cowplot)
library(ggpubr)


# Directories
code_dir <- "~/Wasney-Briscoe/scripts/figures/"
data_dir <-  "~/"

# Functions
source(paste0(code_dir,"helper_scripts/humanized_mouse_utilities.R"))
source(paste0(code_dir,"helper_scripts/strain_phasing_helper_scripts.R"))

focal_species <- c("Alistipes_shahii_62199",
                   "Anaerostipes_hadrus_55206",
                   "Bacteroides_ovatus_58035",
                   "Clostridiales_bacterium_61057",
                   "Coprococcus_comes_61587",
                   "Eubacterium_hallii_61477",
                   "Faecalibacterium_prausnitzii_61481",
                   "Ruminococcus_obeum_61472",
                   "Ruminococcus_sp_55468",
                   "Sutterella_wadsworthensis_56828" 
)
p1 <- plot_strains(focal_species[1])
p2 <- plot_strains(focal_species[2])
p3 <- plot_strains(focal_species[3])
p4 <- plot_strains(focal_species[4])
p5 <- plot_strains(focal_species[5])
p6 <- plot_strains(focal_species[6])
p7 <- plot_strains(focal_species[7])
p8 <- plot_strains(focal_species[8])
p9 <- plot_strains(focal_species[9])
p10 <- plot_strains(focal_species[10])


Fig_s2_strain_legend <- plot_strains("Sutterella_wadsworthensis_56828", output = "legend")

# FINAL GRID #
y_grid_interval <- 1/10
# strain_fig_offset <- 0.012
strain_fig_offset <- 0.015
# legend_offset <- 1/12
legend_offset <- 1/15

figure_s3 <- ggdraw() +
  draw_plot(p1, strain_fig_offset, y_grid_interval*9, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p2, strain_fig_offset, y_grid_interval*8, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p3, strain_fig_offset, y_grid_interval*7, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p4, strain_fig_offset, y_grid_interval*6, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p5, strain_fig_offset, y_grid_interval*5, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p6, strain_fig_offset, y_grid_interval*4, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p7, strain_fig_offset, y_grid_interval*3, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p8, strain_fig_offset, y_grid_interval*2, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p9, strain_fig_offset, y_grid_interval*1, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p10, strain_fig_offset, y_grid_interval*0, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*10, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*9, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*8, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*7, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*6, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*5, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*4, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*3, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*2, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*1, legend_offset, y_grid_interval) +
  draw_plot(Fig_s2_strain_legend, 1-legend_offset, y_grid_interval*0, legend_offset, y_grid_interval) +
  draw_label("Frequency", x = 0.009, y = 1/2, angle = 90, size = 12, fontfamily = "Helvetica")

out_file <- "~/figure_S3.png"
ggsave(out_file, figure_s3, bg = "white", dpi = 300, width = 7.5, height = 16 + 2/3, units = "in")

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


#############################################################################
# B. vulgatus strain plot ###################################################
#############################################################################
# parameters

text_size = 12
subtext_size = 10

# Helper scripts

code_dir <- "~/Wasney-Briscoe/scripts/figures/"
data_dir <-  "~/"
source(paste0(code_dir,"helper_scripts/strain_phasing_helper_scripts.R"))

# Strain plots #

B_vulgatus_strain <- plot_strains("Bacteroides_vulgatus_57955")

B_vulgatus_strain_legend <- plot_strains("Bacteroides_vulgatus_57955", output = "legend")

# FINAL GRID #
y_grid_interval <- 1/1
strain_fig_offset <- 0.015
legend_offset <- 1/15


B_vulgatus_strain_final <- ggdraw() +
  draw_plot(B_vulgatus_strain, strain_fig_offset, y_grid_interval*5, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(B_vulgatus_strain_legend, 1-legend_offset, y_grid_interval*5, legend_offset, y_grid_interval) +
  draw_label("Frequency", x = 0.009, y = y_grid_interval*2 + (1-y_grid_interval*2)/2, angle = 90, size = 12, fontfamily = "Helvetica")

# figure_4

out_file <- "~/B_vulgatus_strain.png"
ggsave(out_file, B_vulgatus_strain_final, bg = "white", dpi = 300, width = 7.5, height = 10, units = "in")


#############################################################################
# B. vulgatus SNV plot ######################################################
#############################################################################

B_vulgatus_snvs <- plot_snvs("Bacteroides_vulgatus_57955")

B_vulgatus_strain_legend <- plot_snvs("Bacteroides_vulgatus_57955", output = "legend")


y_grid_interval <- 1/1
# strain_fig_offset <- 0.012
strain_fig_offset <- 0.015
# legend_offset <- 1/12
legend_offset <- 1/15

B_vulgatus_snvs_final <- ggdraw() +
  draw_plot(B_vulgatus_snvs, strain_fig_offset, y_grid_interval*5, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(B_vulgatus_strain_legend, 1-legend_offset, y_grid_interval*5, legend_offset, y_grid_interval) +
  draw_label("Frequency", x = 0.009, y = y_grid_interval*2 + (1-y_grid_interval*2)/2, angle = 90, size = 12, fontfamily = "Helvetica")

# B_vulgatus_snvs_final

out_file <- "~/B_vulgatus_snvs.png"
ggsave(out_file, B_vulgatus_snvs_final, bg = "white", dpi = 300, width = 7.5, height = 2, units = "in")



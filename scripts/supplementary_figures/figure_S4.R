# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(scales)
library(ggpubr)

# Directories
code_dir <- "~/Wasney-Briscoe-2024/scripts/figures/"
data_dir <-  "~/"

# Functions
source(paste0(code_dir,"helper_scripts/SNV_plot_helper_scripts.R"))

#####  plot

min_coverage <- 20

focal_species <- c("Alistipes_shahii_62199",
                   "Bacteroides_massiliensis_44749",
                   "Bacteroides_thetaiotaomicron_56941",
                   "Blautia_producta_56315",
                   "Blautia_wexlerae_56130",
                   "Clostridiales_bacterium_52743",
                   "Enterococcus_faecium_56947",
                   "Eubacterium_hallii_61477",
                   "Lachnospiraceae_bacterium_56833",
                   "Ruminococcus_bicirculans_59300",
                   "Ruminococcus_gauvreauii_59033")

p1 <- plot_snvs(focal_species[1])
p2 <- plot_snvs(focal_species[2])
p3 <- plot_snvs(focal_species[3])
p4 <- plot_snvs(focal_species[4])
p5 <- plot_snvs(focal_species[5])
p6 <- plot_snvs(focal_species[6])
p7 <- plot_snvs(focal_species[7])
p8 <- plot_snvs(focal_species[8])
p9 <- plot_snvs(focal_species[9])
p10 <- plot_snvs(focal_species[10])
p11 <- plot_snvs(focal_species[11])

p1_legend <- plot_snvs(focal_species[1], output = "legend")
p2_legend <- plot_snvs(focal_species[2], output = "legend")
p3_legend <- plot_snvs(focal_species[3], output = "legend")
p4_legend <- plot_snvs(focal_species[4], output = "legend")
p5_legend <- plot_snvs(focal_species[5], output = "legend")
p6_legend <- plot_snvs(focal_species[6], output = "legend")
p7_legend <- plot_snvs(focal_species[7], output = "legend")
p8_legend <- plot_snvs(focal_species[8], output = "legend")
p9_legend <- plot_snvs(focal_species[9], output = "legend")
p10_legend <- plot_snvs(focal_species[10], output = "legend")
p11_legend <- plot_snvs(focal_species[11], output = "legend")


y_grid_interval <- 1/11
# strain_fig_offset <- 0.012
strain_fig_offset <- 0.015
# legend_offset <- 1/12
legend_offset <- 1/15

figure_s4 <- ggdraw() +
  draw_plot(p1, strain_fig_offset, y_grid_interval*10, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p2, strain_fig_offset, y_grid_interval*9, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p3, strain_fig_offset, y_grid_interval*8, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p4, strain_fig_offset, y_grid_interval*7, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p5, strain_fig_offset, y_grid_interval*6, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p6, strain_fig_offset, y_grid_interval*5, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p7, strain_fig_offset, y_grid_interval*4, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p8, strain_fig_offset, y_grid_interval*3, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p9, strain_fig_offset, y_grid_interval*2, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p10, strain_fig_offset, y_grid_interval*1, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p11, strain_fig_offset, y_grid_interval*0, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(p1_legend, 1-legend_offset, y_grid_interval*10, legend_offset, y_grid_interval) +
  draw_plot(p2_legend, 1-legend_offset, y_grid_interval*9, legend_offset, y_grid_interval) +
  draw_plot(p3_legend, 1-legend_offset, y_grid_interval*8, legend_offset, y_grid_interval) +
  draw_plot(p4_legend, 1-legend_offset, y_grid_interval*7, legend_offset, y_grid_interval) +
  draw_plot(p5_legend, 1-legend_offset, y_grid_interval*6, legend_offset, y_grid_interval) +
  draw_plot(p6_legend, 1-legend_offset, y_grid_interval*5, legend_offset, y_grid_interval) +
  draw_plot(p7_legend, 1-legend_offset, y_grid_interval*4, legend_offset, y_grid_interval) +
  draw_plot(p8_legend, 1-legend_offset, y_grid_interval*3, legend_offset, y_grid_interval) +
  draw_plot(p9_legend, 1-legend_offset, y_grid_interval*2, legend_offset, y_grid_interval) +
  draw_plot(p10_legend, 1-legend_offset, y_grid_interval*1, legend_offset, y_grid_interval) +
  draw_plot(p11_legend, 1-legend_offset, y_grid_interval*0, legend_offset, y_grid_interval) +
  draw_label("Frequency", x = 0.009, y = 1/2, angle = 90, size = 12, fontfamily = "Helvetica")

out_file <- "~/figure_S4.png"
ggsave(out_file, figure_s4, bg = "white", dpi = 300, width = 7.5, height = 18 + 1/3 , units = "in")

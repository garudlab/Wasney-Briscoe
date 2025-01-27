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


# parameters

text_size = 12
subtext_size = 10

# Helper scripts

code_dir <- "~/Wasney-Briscoe-2024/scripts/figures/"
data_dir <-  "~/"
source(paste0(code_dir,"helper_scripts/strain_phasing_helper_scripts.R"))

# ANOVA #

#Choosing good species

data_single <- read.csv(paste0(data_dir, "popgen_stats/SinglePi_SchloissnigPi_cov4_SiteDownsampled.csv")) #%>%


inoculum_species <- data_single %>%
  filter(sample == "TL1gDNAshort") %>%
  select(species) %>%
  unique() %>%
  pull()

prevalent_species <- data_single %>%
  filter(sample != "TL1gDNAshort", species %in% inoculum_species) %>%
  group_by(species) %>%
  summarise(
    sample_count = n(),
    unique_mice = n_distinct(mouse)
  ) %>%
  filter(sample_count >= 3, unique_mice >= 2) %>%
  select(species) %>%
  pull()

data_single <- data_single %>%
  filter(species %in% prevalent_species) %>%
  rowwise()

species_list <- data_single %>%
  filter(sample == "TL1gDNAshort", genomewide_pi > 1e-3) %>%
  select(species) %>%
  unique() %>%
  pull()

anova_results_full <- data.frame(species = NA, variable = NA, VarExpl = NA, ss = NA, df = NA, pval = NA)

for (species in species_list) {
  print(paste0("Processing ", species))
  # DATA
  ## Paths
  project_folder <- "/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/"
  strain_abundance_path <- paste0(project_folder, "strain_phasing/strain_clusters/", species, "/", species, "_strain_frequency.csv")
  ## Loading the data
  strain_abundance <- read.csv2(strain_abundance_path, sep = "\t") %>%
    mutate(freq = as.numeric(freq)) %>%
    mutate(freq_trans = as.vector(clr(freq)))
  # mutate(freq_trans = asin(sqrt(freq)))
  
  
  ## Select the major strain
  
  major_strain <- strain_abundance %>% 
    group_by(strain) %>%
    summarise(mean_freq = mean(freq)) %>%
    filter(mean_freq == max(mean_freq)) %>%
    pull(strain)
  
  strain_abundance <- strain_abundance %>%
    filter(strain == major_strain, mouse_number != "Inoculum") 
  
  # ANOVA
  
  good_cages <- strain_abundance %>%
    filter(mouse_number != "Inoculum") %>%
    group_by(cage) %>%
    summarise(no_of_mice = n_distinct(mouse_number)) %>%
    rowwise() %>%
    mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
    filter(above_1) %>%
    select(cage) %>%
    pull()
  
  good_mice <- strain_abundance %>%
    filter(mouse_number != "Inoculum") %>%
    group_by(mouse_number) %>%
    summarise(no_of_regions = n_distinct(region)) %>%
    rowwise() %>%
    mutate(above_1 = ifelse(no_of_regions > 0, TRUE, FALSE)) %>%
    # mutate(above_1 = ifelse(no_of_regions > 1, TRUE, FALSE)) %>%
    filter(above_1) %>%
    select(mouse_number) %>%
    pull()
  
  good_regions <- strain_abundance %>%
    filter(mouse_number != "Inoculum") %>%
    group_by(region) %>%
    summarise(no_of_mice = n_distinct(mouse_number)) %>%
    rowwise() %>%
    mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
    filter(above_1) %>%
    select(region) %>%
    pull()
  

  strain_abundance_filtered <- strain_abundance %>%
    filter(region %in% good_regions, mouse_number %in% good_mice, cage %in% good_cages)

  no_cages <- length(strain_abundance_filtered %>% select(cage) %>% unique() %>% pull())
  no_mice <- length(strain_abundance_filtered %>% select(mouse_number) %>% unique() %>% pull())
  no_regions <- length(strain_abundance_filtered %>% select(region) %>% unique() %>% pull())
  
  
  no_cages_old <- -1
  no_mice_old <- -1
  no_regions_old <- -1
  
  while ((no_cages != no_cages_old) & (no_mice != no_cages_old) & (no_regions != no_cages_old)) {
  
    no_cages_old <- no_cages
    no_mice_old <- no_mice
    no_regions_old <- no_regions
    
    good_cages <- strain_abundance_filtered %>%
      filter(mouse_number != "Inoculum") %>%
      group_by(cage) %>%
      summarise(no_of_mice = n_distinct(mouse_number)) %>%
      rowwise() %>%
      mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
      filter(above_1) %>%
      select(cage) %>%
      pull()
    
    good_mice <- strain_abundance_filtered %>%
      filter(mouse_number != "Inoculum") %>%
      group_by(mouse_number) %>%
      summarise(no_of_regions = n_distinct(region)) %>%
      rowwise() %>%
      # mutate(above_1 = ifelse(no_of_regions > 0, TRUE, FALSE)) %>%
      mutate(above_1 = ifelse(no_of_regions > 1, TRUE, FALSE)) %>%
      filter(above_1) %>%
      select(mouse_number) %>%
      pull()
    
    good_regions <- strain_abundance_filtered %>%
      filter(mouse_number != "Inoculum") %>%
      group_by(region) %>%
      summarise(no_of_mice = n_distinct(mouse_number)) %>%
      rowwise() %>%
      mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
      filter(above_1) %>%
      select(region) %>%
      pull()
    
    strain_abundance_filtered <- strain_abundance_filtered %>%
      filter(region %in% good_regions, mouse_number %in% good_mice, cage %in% good_cages)
    
    no_cages <- length(strain_abundance_filtered %>% select(cage) %>% unique() %>% pull())
    no_mice <- length(strain_abundance_filtered %>% select(mouse_number) %>% unique() %>% pull())
    no_regions <- length(strain_abundance_filtered %>% select(region) %>% unique() %>% pull())
    
    
  }
  
  
  if ((no_cages >= 2) & (no_mice >= 2) & (no_regions >= 2)) {
    model <- aov(freq_trans ~ cage + mouse_number + region, data = strain_abundance)
    
  } else if ((no_mice >= 2) & (no_regions >= 2)) {
    model <- aov(freq_trans ~  mouse_number + region, data = strain_abundance)
    
  } else if ((no_mice >= 2)) {
    model <- aov(freq_trans ~  mouse_number, data = strain_abundance)
    
  } else {
    print(paste("Not enough good samples/mice/cages in", species))
    next
  }
  
  variance_explained <- data.frame(summary(model)[[1]])[, 2, drop = FALSE][[1]]/sum(data.frame(summary(model)[[1]])[, 2, drop = FALSE])
  degrees_of_freedom <- data.frame(summary(model)[[1]])[, 1, drop = FALSE][[1]]
  pvals <- data.frame(summary(model)[[1]])[, 5, drop = FALSE][[1]]
  variables <- trimws(rownames(data.frame(summary(model)[[1]])))
  sum_of_squares <- data.frame(summary(model)[[1]])[, 2, drop = FALSE][[1]]
  
  anova_results <- data.frame(species = species, variable = variables, VarExpl = variance_explained, ss = sum_of_squares, df = degrees_of_freedom, pval = pvals)
  
  anova_results_full <- rbind(anova_results_full, anova_results)
  
}

anova_results_full <- anova_results_full %>%
  filter(!is.na(species))

# PLOTTING

species_subset <- anova_results_full %>%
  group_by(species) %>%
  summarize(no_of_rows = n()) %>%
  filter(no_of_rows == 4) %>%
  select(species) %>%
  pull() %>%
  sort()

anova_df <- anova_results_full %>%
  filter(species %in% species_subset) %>%
  rowwise() %>%
  mutate(variable = ifelse(variable == "mouse_number", "Mouse", variable),
         species = paste0(strsplit(species, "_")[[1]][1], " ", strsplit(species, "_")[[1]][2])) %>%
  mutate(variable = ifelse(variable == "region", "Gut region", variable)) %>%
  mutate(variable = ifelse(variable == "cage", "Cage", variable))

anova_df$variable = factor(anova_df$variable, levels = c("Cage","Mouse", "Gut region", "Residuals"))

A <- ggplot(anova_df, aes(x = (species), y = VarExpl*100, fill = variable, group = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Cage" = "#483248", # Eggplant
                               "Mouse" = "#AA98A9", #  Lilac
                               "Gut region" = "#BDB5D5", # Wisteria
                               "Residuals" = "lightgrey")) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = text_size),
        axis.text.x = element_text(angle = 45, hjust = 1, size = subtext_size),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = subtext_size),
        axis.title.y = element_text(size = text_size),
        title = element_text(size = text_size),
        legend.title = element_blank(),
        legend.text = element_text(size=subtext_size),
        plot.margin = unit(c(0.5, 0.5, 0.5, 2.5), "lines"),
        # legend.position = "none"
        # legend.position = c(0.2, 0.8)
  ) +
  labs(title = "Partitioning observed variance in strain abundance",x = "Species", y = "Variance explained (%)", fill = "Variable")

# Strain plots #

A_equolifaciens <- plot_strains("Adlercreutzia_equolifaciens_60310")
A_hadrus <- plot_strains("Anaerostipes_hadrus_55206")
B_uniformis <- plot_strains("Bacteroides_uniformis_57318")
B_vulgatus <- plot_strains("Bacteroides_vulgatus_57955")
B_wexlerae <- plot_strains("Blautia_wexlerae_56130")
C_bacterium <- plot_strains("Clostridiales_bacterium_61057")
E_hallii <- plot_strains("Eubacterium_hallii_61477")

Fig_4_strain_legend <- plot_strains("Adlercreutzia_equolifaciens_60310", output = "legend")


# FINAL GRID #
y_grid_interval <- 1/6
strain_fig_offset <- 0.015
legend_offset <- 1/15


figure_4 <- ggdraw() +
  draw_plot(A_equolifaciens, strain_fig_offset, y_grid_interval*5, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(B_uniformis, strain_fig_offset, y_grid_interval*4, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(B_vulgatus, strain_fig_offset, y_grid_interval*3, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(B_wexlerae, strain_fig_offset, y_grid_interval*2, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(Fig_4_strain_legend, 1-legend_offset, y_grid_interval*5, legend_offset, y_grid_interval) +
  draw_plot(Fig_4_strain_legend, 1-legend_offset, y_grid_interval*4, legend_offset, y_grid_interval) +
  draw_plot(Fig_4_strain_legend, 1-legend_offset, y_grid_interval*3, legend_offset, y_grid_interval) +
  draw_plot(Fig_4_strain_legend, 1-legend_offset, y_grid_interval*2, legend_offset, y_grid_interval) +
  draw_plot(A, 0, 0, 1, y_grid_interval*2) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1,  y_grid_interval*2), size = 16, family = "Helvetica") +
  draw_label("Frequency", x = 0.009, y = y_grid_interval*2 + (1-y_grid_interval*2)/2, angle = 90, size = 12, fontfamily = "Helvetica")

# figure_4

out_file <- "~/figure_4.png"
ggsave(out_file, figure_4, bg = "white", dpi = 300, width = 7.5, height = 10, units = "in")



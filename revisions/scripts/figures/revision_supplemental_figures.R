#########################################################################################################
# Conventional mouse strain figure                                                                    #
#########################################################################################################
# remove list
rm(list = ls())

# Packages
library(dplyr)
library(cowplot)

# Functions

source("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/figures/helper_scripts/strain_phasing_helper_scripts.R")
source("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/helper_scripts/family_plot_helpers.R")

##### Selecting species ##### 


# Load species list
strain_phasing_dir <- "/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/strain_phasing/strain_clusters/"
species_list <- list.dirs(strain_phasing_dir, full.names = FALSE)[-1]

# constructing ANOVA table
anova_results_full <- data.frame(species = NA, variable = NA, VarExpl = NA, ss = NA, df = NA, pval = NA)


for (species in species_list) {
  print(paste0("Processing ", species))
  # DATA
  ## Paths
  project_folder <- "/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/"
  strain_abundance_path <- paste0(project_folder, "strain_phasing/strain_clusters/", species, "/", species, "_strain_frequency.csv")
  ## Loading the data
  strain_abundance <- read.csv2(strain_abundance_path, sep = "\t") %>%
    mutate(freq = as.numeric(freq)) %>%
    mutate(freq_trans = as.vector(clr(freq)))
  # mutate(freq_trans = asin(sqrt(freq)))
  
  # rename mice
  strain_abundance <- strain_abundance %>%
    rowwise() %>%
    mutate(mouse_number = paste0("Mouse ",as.character(extract_mouse(`sample`))))
  
  
  ## Select the major strain
  
  major_strain <- strain_abundance %>% 
    group_by(strain) %>%
    summarise(mean_freq = mean(freq)) %>%
    filter(mean_freq == max(mean_freq)) %>%
    pull(strain)
  
  # ANOVA
  
  good_cages <- strain_abundance %>%
    group_by(cage) %>%
    summarise(no_of_mice = n_distinct(mouse_number)) %>%
    rowwise() %>%
    mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
    filter(above_1) %>%
    select(cage) %>%
    pull()
  
  good_mice <- strain_abundance %>%
    group_by(mouse_number) %>%
    summarise(no_of_regions = n_distinct(region)) %>%
    rowwise() %>%
    mutate(above_1 = ifelse(no_of_regions > 0, TRUE, FALSE)) %>%
    # mutate(above_1 = ifelse(no_of_regions > 1, TRUE, FALSE)) %>%
    filter(above_1) %>%
    select(mouse_number) %>%
    pull()
  
  good_regions <- strain_abundance %>%
    group_by(region) %>%
    summarise(no_of_mice = n_distinct(mouse_number)) %>%
    rowwise() %>%
    mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
    filter(above_1) %>%
    select(region) %>%
    pull()
  
  strain_abundance_filtered <- strain_abundance %>%
    filter(region %in% good_regions, mouse_number %in% good_mice, cage %in% good_cages)
  
  # no_cages <- length(strain_abundance %>% select(cage) %>% unique() %>% pull())
  # no_mice <- length(strain_abundance %>% select(mouse_number) %>% unique() %>% pull())
  # no_regions <- length(strain_abundance %>% select(region) %>% unique() %>% pull())
  
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

##### A. Strain phasing 100555 #####

A <- plot_strains("100555", output = "normal")

##### B. Strain phasing 100600 #####

B <- plot_strains("100600", output = "normal")

##### C. Strain phasing 213583 #####

C <- plot_strains("213583", output = "normal")

##### D. Strain phasing 214603 #####

D <- plot_strains("214603", output = "normal")

##### E. Strain phasing 253991 #####

E <- plot_strains("253991", output = "normal")

##### F. Strain phasing 261672 #####

F_plot <- plot_strains("261672", output = "normal")

##### Grid #####

strain_fig_offset <- 0.0333


rect1 <- rectGrob(
  x = 0.085, y = 5/6,            # center of panels 1-3
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect2 <- rectGrob(
  x = 0.359, y = 5/6,            # center of panels 4-6
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#AA98A9", col = NA),
  just = c("left", "bottom")
)

rect3 <- rectGrob(
  x = 0.6325, y = 5/6,            # center of panels 4-6
  width = 0.275, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect4 <- rectGrob(
  x = 0.085, y = 4/6,            # center of panels 1-3
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect5 <- rectGrob(
  x = 0.359, y = 4/6,            # center of panels 4-6
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#AA98A9", col = NA),
  just = c("left", "bottom")
)

rect6 <- rectGrob(
  x = 0.6325, y = 4/6,            # center of panels 4-6
  width = 0.275, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect7 <- rectGrob(
  x = 0.085, y = 3/6,            # center of panels 1-3
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect8 <- rectGrob(
  x = 0.359, y = 3/6,            # center of panels 4-6
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#AA98A9", col = NA),
  just = c("left", "bottom")
)

rect9 <- rectGrob(
  x = 0.6325, y = 3/6,            # center of panels 4-6
  width = 0.275, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect10 <- rectGrob(
  x = 0.085, y = 2/6,            # center of panels 1-3
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect11 <- rectGrob(
  x = 0.359, y = 2/6,            # center of panels 4-6
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#AA98A9", col = NA),
  just = c("left", "bottom")
)

rect12 <- rectGrob(
  x = 0.6325, y = 2/6,            # center of panels 4-6
  width = 0.275, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect13 <- rectGrob(
  x = 0.085, y = 1/6,            # center of panels 1-3
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect14 <- rectGrob(
  x = 0.359, y = 1/6,            # center of panels 4-6
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#AA98A9", col = NA),
  just = c("left", "bottom")
)

rect15 <- rectGrob(
  x = 0.6325, y = 1/6,            # center of panels 4-6
  width = 0.275, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect16 <- rectGrob(
  x = 0.085, y = 0/6,            # center of panels 1-3
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)

rect17 <- rectGrob(
  x = 0.359, y = 0/6,            # center of panels 4-6
  width = 0.3, height = 0.1425,
  gp = gpar(fill = "#AA98A9", col = NA),
  just = c("left", "bottom")
)

rect18 <- rectGrob(
  x = 0.6325, y = 0/6,            # center of panels 4-6
  width = 0.275, height = 0.1425,
  gp = gpar(fill = "#e7e1ef", col = NA),
  just = c("left", "bottom")
)


grid <- ggdraw() +
  draw_grob(rect1) +
  draw_grob(rect2) +
  draw_grob(rect3) +
  draw_grob(rect4) +
  draw_grob(rect5) +
  draw_grob(rect6) +
  draw_grob(rect7) +
  draw_grob(rect8) +
  draw_grob(rect9) +
  draw_grob(rect10) +
  draw_grob(rect11) +
  draw_grob(rect12) +
  draw_grob(rect13) +
  draw_grob(rect14) +
  draw_grob(rect15) +
  draw_grob(rect16) +
  draw_grob(rect17) +
  draw_grob(rect18) +
  draw_plot(A, strain_fig_offset,5/6,1-strain_fig_offset,1/6) +
  draw_plot(B, strain_fig_offset,4/6,1-strain_fig_offset,1/6) +
  draw_plot(C, strain_fig_offset,3/6,1-strain_fig_offset,1/6) +
  draw_plot(D, strain_fig_offset,2/6,1-strain_fig_offset,1/6) +
  draw_plot(E, strain_fig_offset,1/6,1-strain_fig_offset,1/6) +
  draw_plot(F_plot, strain_fig_offset,0/6,1-strain_fig_offset,1/6)
  # draw_plot_label(c("A", "B", "C", "D", "E","F"), c(0,0,0,0,0,0), c(1,5/6,4/6,3/6,2/6,1/6), size = 16, family = "Helvetica")


out_file <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/mouse_strain_phases.png"
ggsave(out_file, grid, bg = "white", dpi = 300, width = 7.5, height = 10, units = "in")


#########################################################################################################
# Shalon strain phases                                                                  #
#########################################################################################################

# remove list
rm(list = ls())

# Packages
library(dplyr)
library(cowplot)

source("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/figures/helper_scripts/strain_phasing_helper_scripts.R")

##### Strain phasing Acidaminococcus_intestini_54097, subject 9 ##### 

A <- plot_strains("Acidaminococcus_intestini_54097", subject = 9, brief_x_labels = TRUE) # 2 panels

##### Strain phasing Adlercreutzia_equolifaciens_60310, subject 12 #####

B <- plot_strains("Adlercreutzia_equolifaciens_60310", subject = 12, brief_x_labels  = TRUE)

##### Strain phasing Anaerostipes_hadrus_55206, subject 6 #####

C <- plot_strains("Anaerostipes_hadrus_55206", subject = 6, brief_x_labels = TRUE, timepoints_to_use = c(3,5,6,7))

##### Strain phasing Anaerostipes_hadrus_55206, subject 8 #####

D <- plot_strains("Anaerostipes_hadrus_55206", subject = 8, brief_x_labels = TRUE, timepoints_to_use = c(2,5,6,8))

##### Strain phasing Bacteroides_vulgatus_57955, subject 2 #####

E <- plot_strains("Bacteroides_vulgatus_57955", subject = 2, brief_x_labels = TRUE, timepoints_to_use = c(3,5,6,7)) #, 

##### Strain phasing Bacteroides_vulgatus_57955, subject 8 #####

F_plot <- plot_strains("Bacteroides_vulgatus_57955", subject = 8, brief_x_labels = TRUE, timepoints_to_use = c(3,5,7,9)) #, , 

##### Strain phasing Bacteroides_vulgatus_57955, subject 11 #####

G <- plot_strains("Bacteroides_vulgatus_57955", subject = 11, brief_x_labels = TRUE, timepoints_to_use = c(3,5,6,10)) #, , 

##### Strain phasing Bifidobacterium_adolescentis_56815, subject 1 #####

H <- plot_strains("Bifidobacterium_adolescentis_56815", subject = 1, brief_x_labels = TRUE, timepoints_to_use = c(1,2,3,4)) #, 

##### Strain phasing Bifidobacterium_longum_57796, subject 5 #####

I <- plot_strains("Bifidobacterium_longum_57796", subject = 5, brief_x_labels = TRUE, timepoints_to_use = c(1,4)) # 2 panels

##### Strain phasing Bilophila_wadsworthia_57364, subject 9 #####

J <- plot_strains("Bilophila_wadsworthia_57364", subject = 9, brief_x_labels = TRUE) 

##### Strain phasing Blautia_wexlerae_56130, subject 8 #####

K <- plot_strains("Blautia_wexlerae_56130", subject = 8, brief_x_labels = TRUE, timepoints_to_use = c(2,5,7,8)) 

##### Strain phasing Burkholderiales_bacterium_56577, subject 6 #####

L <- plot_strains("Burkholderiales_bacterium_56577", subject = 6, brief_x_labels = TRUE, timepoints_to_use = c(1,3,4,5)) 

##### Strain phasing Catenibacterium_mitsuokai_61547, subject 10 #####

M <- plot_strains("Catenibacterium_mitsuokai_61547", subject = 10, brief_x_labels = TRUE, timepoints_to_use = c(3,4,5,6)) 

##### Strain phasing Dorea_formicigenerans_56346, subject 12 #####

N <- plot_strains("Dorea_formicigenerans_56346", subject = 12, brief_x_labels = TRUE, timepoints_to_use = c(3,4,5))  # 3 panels

##### Strain phasing Enterobacter_cloacae_58148, subject 13 #####

O <- plot_strains("Enterobacter_cloacae_58148", subject = 13, brief_x_labels = TRUE, timepoints_to_use = c(4,5,6,7)) 

##### Strain phasing Eubacterium_rectale_56927, subject 1 #####

P <- plot_strains("Eubacterium_rectale_56927", subject = 9, brief_x_labels = TRUE, timepoints_to_use = c(5,6,7)) # 3 panels

##### Strain phasing Eubacterium_rectale_56927, subject 14 #####

Q <- plot_strains("Eubacterium_rectale_56927", subject = 14, brief_x_labels = TRUE, timepoints_to_use = c(6,8,10))  # 3 panels

##### Strain phasing Guyana_massiliensis_60772, subject 12 #####

R <- plot_strains("Guyana_massiliensis_60772", subject = 12, brief_x_labels = TRUE, timepoints_to_use = c(1,3,4,5))

##### Strain phasing Parabacteroides_distasonis_56985, subject 9 #####

S <- plot_strains("Parabacteroides_distasonis_56985", subject = 9, brief_x_labels = TRUE, timepoints_to_use = c(1,2,3,5))

##### Strain phasing Ruminococcus_obeum_61472, subject 6 #####

Tee <- plot_strains("Ruminococcus_obeum_61472", subject = 6, brief_x_labels = TRUE) # 3 pabels

##### Strain phasing Ruminococcus_obeum_61472, subject 11 #####

U <- plot_strains("Ruminococcus_obeum_61472", subject = 11, brief_x_labels = TRUE, timepoints_to_use = c(3,5,6,8)) #

##### Strain phasing Ruminococcus_obeum_61472, subject 13 #####

V <- plot_strains("Ruminococcus_obeum_61472", subject = 13, brief_x_labels = TRUE, timepoints_to_use = c(3,6,7,8)) 

##### Strain phasing Ruminococcus_sp_55468, subject 6 #####

W <- plot_strains("Ruminococcus_sp_55468", subject = 6, brief_x_labels = TRUE, timepoints_to_use = c(1,2,4,6))

##### Strain phasing Ruminococcus_torques_62045, subject 8 #####

X <- plot_strains("Ruminococcus_torques_62045", subject = 8, brief_x_labels = TRUE, timepoints_to_use = c(3,5,6,7))

##### Grid 1 #####

text_size = 12
subtext_size = 10
strain_fig_offset <- 0.0333
x_axis_title_offset <- 0.013
panel_size = 0.196

grid <- ggdraw() +
  draw_plot(A, strain_fig_offset,x_axis_title_offset+5*(1-x_axis_title_offset)/6, panel_size*3-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(B, strain_fig_offset,x_axis_title_offset+4*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(C, strain_fig_offset,x_axis_title_offset+3*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(D, strain_fig_offset,x_axis_title_offset+2*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(E, strain_fig_offset,x_axis_title_offset+1*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(F_plot, strain_fig_offset,x_axis_title_offset+0*(1-x_axis_title_offset)/6,1-strain_fig_offset,1/6) +
  draw_label("Device type", x = 0.4375 , y = 0.005,hjust = 0,vjust = 0.0,size = text_size,fontfamily = "Helvetica", color = "black")
# draw_plot_label(c("A", "B", "C", "D", "E","F"), c(0,0,0,0,0,0), c(1,5/6,4/6,3/6,2/6,1/6), size = 16, family = "Helvetica")


out_file <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/Shalon_strain_phase_1.png"
ggsave(out_file, grid, bg = "white", dpi = 300, width = 7.5, height = 10, units = "in")


##### Grid 2 #####

text_size = 12
subtext_size = 10
strain_fig_offset <- 0.0333
x_axis_title_offset_2 <- x_axis_title_offset*2
strain_fig_offset <- 0.0333
panel_size = 0.196


grid_2 <- ggdraw() +
  draw_plot(G, strain_fig_offset,x_axis_title_offset+5*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(H, strain_fig_offset,x_axis_title_offset+4*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(I, strain_fig_offset,x_axis_title_offset+3*(1-x_axis_title_offset)/6,panel_size*3-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(J, strain_fig_offset,x_axis_title_offset+2*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(K, strain_fig_offset,x_axis_title_offset+1*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(L, strain_fig_offset,x_axis_title_offset+0*(1-x_axis_title_offset)/6,1-strain_fig_offset,1/6) +
  draw_label("Device type", x = 0.4375 , y = 0.005,hjust = 0,vjust = 0.0,size = text_size,fontfamily = "Helvetica", color = "black")
# draw_plot_label(c("A", "B", "C", "D", "E","F"), c(0,0,0,0,0,0), c(1,5/6,4/6,3/6,2/6,1/6), size = 16, family = "Helvetica")

out_file <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/Shalon_strain_phase_2.png"
ggsave(out_file, grid_2, bg = "white", dpi = 300, width = 7.5, height = 10, units = "in")


##### Grid 3 #####

text_size = 12
subtext_size = 10
strain_fig_offset <- 0.0333
strain_fig_offset <- 0.0333
panel_size = 0.196


grid_3 <- ggdraw() +
  draw_plot(N, strain_fig_offset,x_axis_title_offset+5*(1-x_axis_title_offset)/6,panel_size*4.05-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(O, strain_fig_offset,x_axis_title_offset+4*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(P, strain_fig_offset,x_axis_title_offset+3*(1-x_axis_title_offset)/6,panel_size*4.05-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(Q, strain_fig_offset,x_axis_title_offset+2*(1-x_axis_title_offset)/6,panel_size*4.05-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(R, strain_fig_offset,x_axis_title_offset+1*(1-x_axis_title_offset)/6,1-strain_fig_offset,(1-x_axis_title_offset)/6) +
  draw_plot(S, strain_fig_offset,x_axis_title_offset+0*(1-x_axis_title_offset)/6,1-strain_fig_offset,1/6) +
  draw_label("Device type", x = 0.4375 , y = 0.005,hjust = 0,vjust = 0.0,size = text_size,fontfamily = "Helvetica", color = "black")
# draw_plot_label(c("A", "B", "C", "D", "E","F"), c(0,0,0,0,0,0), c(1,5/6,4/6,3/6,2/6,1/6), size = 16, family = "Helvetica")


out_file <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/Shalon_strain_phase_3.png"
ggsave(out_file, grid_3, bg = "white", dpi = 300, width = 7.5, height = 10, units = "in")

##### Grid 4 #####

text_size = 12
subtext_size = 10
strain_fig_offset <- 0.0333
x_axis_title_offset_2 <- x_axis_title_offset*1.25
strain_fig_offset <- 0.0333
panel_size = 0.196


grid_4 <- ggdraw() +
  draw_plot(Tee, strain_fig_offset,x_axis_title_offset_2+3*(1-x_axis_title_offset_2)/4,panel_size*4.05-strain_fig_offset,(1-x_axis_title_offset_2)/4) +
  draw_plot(V, strain_fig_offset,x_axis_title_offset_2+2*(1-x_axis_title_offset_2)/4,1-strain_fig_offset,(1-x_axis_title_offset_2)/4) +
  draw_plot(W, strain_fig_offset,x_axis_title_offset_2+1*(1-x_axis_title_offset_2)/4,1-strain_fig_offset,(1-x_axis_title_offset_2)/4) +
  draw_plot(X, strain_fig_offset,x_axis_title_offset_2+0*(1-x_axis_title_offset_2)/4,1-strain_fig_offset,(1-x_axis_title_offset_2)/4) +
  draw_label("Device type", x = 0.4375 , y = 0.0075,hjust = 0,vjust = 0.0,size = text_size,fontfamily = "Helvetica", color = "black")
# draw_plot_label(c("A", "B", "C", "D", "E","F"), c(0,0,0,0,0,0), c(1,5/6,4/6,3/6,2/6,1/6), size = 16, family = "Helvetica")


out_file <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/Shalon_strain_phase_4.png"
ggsave(out_file, grid_4, bg = "white", dpi = 300, width = 7.5, height = (10/6)*5, units = "in")


##### SCRATCH #####

p <- plot_strains("Anaerostipes_hadrus_55206", subject = 6, timepoints_to_use = c(2,3,5,6,7))
out_file <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/Ahadrus_Subject6_demo.png"
ggsave(out_file, p, bg = "white", dpi = 300, width = 7.5, height = 3, units = "in")


#########################################################################################################
# Supplementary figure X (correlation)                                                                  #
#########################################################################################################

# Clearing environment
rm(list=ls()) 

# packages
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
library(scales)
library(latex2exp)


# parameters

text_size = 12
subtext_size = 10

# Directories

study <- "HumanizedMouse_Batch2" 
code_dir <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/HumanizedMouse/figures/figures/"
data_dir <-  "/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/"
# good_species_path <- paste0(data_dir, "metadata/good_cvg_species_list.txt")

source(paste0(code_dir,"helper_scripts/humanized_mouse_utilities.R"))
source(paste0(code_dir, "helper_scripts/mouse_annotation_functions.R"))
source(paste0(code_dir,"helper_scripts/mouse_annotation_functions.R"))


# Loading the data

## Good species
# good_species <- read.csv(good_species_path, header = FALSE)$V1

## Data

data_single <- read.csv(paste0(data_dir, "popgen_stats/SinglePi_SchloissnigPi_cov4_SiteDownsampled.csv")) #%>%
# mutate(good_species = as.logical(good_species))

## Filtering the data
data_single <- data_single %>%
  mutate(total_loci = as.numeric(total_loci)) %>%
  # filter(good_species == TRUE) %>%
  filter(total_loci > 5e5) %>%
  rowwise() %>%
  mutate(species_original = species,
         species = if_else(
           str_detect(species, "Ruminococcus_sp") | str_detect(species, "Faecalibacterium_prausnitzii"),
           paste(str_split(gsub("_", " ", species), " ")[[1]][1:3], collapse = " "),  # Keep 3 elements for "Ruminococcus"
           paste(str_split(gsub("_", " ", species), " ")[[1]][1:2], collapse = " ")   # Keep 2 elements otherwise
         ))

## Choosing good species

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
  rowwise() %>%
  mutate(species = if_else(
    str_detect(species, "Faecalibacterium prausnitzii"),
    paste(str_split(gsub(" ", " ", species), " ")[[1]][1:2], collapse = " "),  # Keep 3 elements for "Ruminococcus"
    species 
  ))


## Calculate species order
species_order <- data_single %>%
  filter(sample == "TL1gDNAshort") %>%
  arrange(desc(genomewide_pi))  %>%
  select(species) %>%
  pull()

species_order_rev <- data_single %>%
  filter(sample == "TL1gDNAshort") %>%
  arrange(genomewide_pi)  %>%
  select(species) %>%
  pull()

## Merging the data
data_single <- data_single %>% 
  mutate(pi_type = ifelse(sample == "TL1gDNAshort","Inoculum pi","Intra-sample pi"))

all_pi_data <- data_single

## Summary stats
average_labels <- all_pi_data %>% group_by(species,species_original,pi_type) %>% 
  summarise(mean_pi = mean(genomewide_pi)) %>%
  filter(pi_type == "Inoculum pi") 
all_pi_data <- all_pi_data %>%
  rowwise() %>%
  mutate(pi_type = ifelse(pi_type == "Intra-sample pi", "Intra-sample pi (average)", pi_type))


## annotating/filtering data

all_pi_data <- all_pi_data %>%
  filter(pi_type == "Intra-sample pi (average)", sample != "TL1gDNAshort") %>% 
  rowwise() %>%
  mutate(gut_segment = gut_site)
# rowwise() %>%
# mutate(mouse = extract_mouse_number(sample),
#        gut_segment = extract_tissue_type(sample, metadata))



merged_pi_df <- all_pi_data %>%
  left_join(average_labels, by = "species") %>%
  rename(mouse_pi = genomewide_pi, inoculum_pi = mean_pi) %>%
  select(c(species,sample,mouse,cage,diet,gut_site,gut_region,mouse_pi,inoculum_pi))

avg_pi_df <- merged_pi_df %>%
  select(species,mouse_pi, inoculum_pi) %>%
  group_by(species) %>%
  summarise(mouse_pi = mean(mouse_pi),inoculum_pi = mean(inoculum_pi))

pi_corr <- cor(merged_pi_df$mouse_pi, merged_pi_df$inoculum_pi, method = "spearman")



avg_pi_corr <- cor(avg_pi_df$mouse_pi, avg_pi_df$inoculum_pi, method = "spearman")
pi_test <- cor.test(avg_pi_df$mouse_pi,
                    avg_pi_df$inoculum_pi,
                    method = "spearman")


mouse_summary <- merged_pi_df %>%
  group_by(species) %>%
  summarise(
    mean_mouse_pi = mean(mouse_pi, na.rm = TRUE),
    sd_mouse_pi = sd(mouse_pi, na.rm = TRUE),
    inoculum_pi = unique(inoculum_pi)
  ) %>%
  arrange(desc(inoculum_pi)) %>%
  mutate(species = factor(species, levels = species))

plot_df <- mouse_summary %>%
  pivot_longer(
    cols = c(mean_mouse_pi, inoculum_pi),
    names_to = "source",
    values_to = "pi_value"
  ) %>%
  mutate(source = factor(source,
                         levels = c("mean_mouse_pi", "inoculum_pi"),
                         labels = c("Mouse mean pi", "Inoculum pi")))

plot_df_wide <- plot_df %>%
  pivot_wider(id_cols = species, names_from = `source`, values_from = pi_value) %>%
  rename(inoculum_pi = `Inoculum pi`, mouse_mean_pi = `Mouse mean pi`)

ols_model <- lm(mouse_mean_pi ~ inoculum_pi, data = plot_df_wide)
summary(ols_model)


p <- ggplot(plot_df, aes(x = species)) +
  geom_errorbar(
    data = subset(mouse_summary),
    aes(
      ymin = mean_mouse_pi - sd_mouse_pi,
      ymax = mean_mouse_pi + sd_mouse_pi
    ),
    width = 0.2,
    color = "black"
  ) +
  geom_point(aes(y = pi_value, color = source), size = 3) +
  labs(
    x = "Species",
    y = expression(pi)
  ) +
  scale_color_manual(
    values = c("Mouse mean pi" = "black", "Inoculum pi" = "brown"),
    labels = c(expression("Average mouse " * pi), expression("Inoculum " * pi))  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = subtext_size),
    axis.text.y = element_text(size = subtext_size),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = text_size),
    text = element_text(family = "Helvetica"),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.title = element_blank(),
    legend.background = element_rect(fill = alpha("white", 0.7), color = "black"),
    legend.text.align = 0 ,
    plot.margin = margin(10, 10, 10, 20)
  )
p

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/pi_corr_fig_old.png"
ggsave(out_path, p, dpi = 300, height = 5, width = 7.5)


p_1 <- ggplot() +
  # Error bars for mouse mean ± SD
  geom_errorbar(
    data = mouse_summary,
    aes(
      x = inoculum_pi,
      ymin = mean_mouse_pi - sd_mouse_pi,
      ymax = mean_mouse_pi + sd_mouse_pi
    ),
    width = 0.0001,
    color = "black"
  ) +
  # Points for mean mouse π
  geom_point(
    data = plot_df_wide,
    aes(x = inoculum_pi, y = mouse_mean_pi),
    fill = "red",          
    shape = 21,            
    size = 3,              
    stroke = 0.5
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +   # y = x line
  labs(
    x = expression("Inoculum " * pi),
    y = expression("Average mouse " * pi)
  ) +
  annotate(
    "text",
    x = 0.0055,
    y = 0.006,  # a bit above the line
    label = "y = x",
    color = "gray40",
    size = 4,
    hjust = 1,
    family = "Helvetica",
    angle = 39
  ) +
  annotate(
    "text",
    x = max(plot_df_wide$inoculum_pi, na.rm = TRUE) *0.65,
    y = max(plot_df_wide$mouse_mean_pi, na.rm = TRUE) * 1.6,
    label = expression("Spearman's " * rho * " = 0.845"),
    color = "black",
    family = "Helvetica",
    size = 4,
    hjust = 0
  ) +
  annotate(
    "text",
    x = max(plot_df_wide$inoculum_pi, na.rm = TRUE) *0.65,
    y = max(plot_df_wide$mouse_mean_pi, na.rm = TRUE) * 1.5,
    label = expression("p-value = " ~ 7.549 %*% 10^-7),
    color = "black",
    family = "Helvetica",
    size = 4,
    hjust = 0
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica"),
    legend.title = element_blank(),
    legend.text.align = 0,
    axis.title = element_text(size = text_size),
    axis.text.x = element_text(size = subtext_size),
    axis.text.y = element_text(size = subtext_size)
  )

# p_1

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/pi_corr_fig.png"
ggsave(out_path, p_1, dpi = 300, height = 5, width = 7.5)



#########################################################################################################
# Supplementary figure 4 (variance)                                                                  #
#########################################################################################################

# remove list
rm(list = ls())

# Packages
library(dplyr)

# species subset

species_subset <- c("Adlercreutzia_equolifaciens_60310",
                    "Anaerostipes_hadrus_55206",
                    "Bacteroides_uniformis_57318",
                    "Bacteroides_vulgatus_57955",
                    "Blautia_wexlerae_56130",
                    "Clostridiales_bacterium_61057",
                    "Eubacterium_hallii_61477")



within_cage_variance_df <- data.frame(
  species = character(),
  cage = character(),
  region = character(),
  variance = numeric()
)

between_cage_variance_df <- data.frame(
  species = character(),
  region = character(),
  variance = numeric()
)

for (species in species_subset) {
  project_folder <- "/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/"
  strain_abundance_path <- paste0(project_folder, "strain_phasing/strain_clusters/", species, "/", species, "_strain_frequency.csv")
  
  strain_abundance <- read.csv2(strain_abundance_path, sep = "\t") %>%
    mutate(freq = as.numeric(freq)) 
  
  within_cage_variance <- strain_abundance %>%
    filter(cage != "Inoculum", strain == 1) %>%
    group_by(cage, region) %>%
    summarise(variance = stats::var(freq)) %>%
    arrange(region, desc(variance)) %>%
    mutate(`species` = species)
  
  between_cage_variance <- strain_abundance %>%
    filter(cage != "Inoculum", strain == 1) %>%
    group_by(region) %>%
    summarise(variance = stats::var(freq)) %>%
    arrange(region, desc(variance)) %>%
    mutate(`species` = species)
  
  within_cage_variance_df <- rbind(within_cage_variance_df, within_cage_variance)
  between_cage_variance_df <- rbind(between_cage_variance_df, between_cage_variance)
  
}

between_cage_avg_variance_df <- between_cage_variance_df %>%
  group_by(`species`) %>%
  summarise(variance = mean(variance, na.rm = TRUE)) %>%
  mutate(variance_type = "Between cage")

within_cage_avg_variance_df <- within_cage_variance_df %>%
  group_by(`species`) %>%
  summarise(variance = mean(variance, na.rm = TRUE)) %>%
  mutate(variance_type = "Within cage")

variance_df <- rbind(between_cage_avg_variance_df, within_cage_avg_variance_df) %>%
  rowwise() %>%
  mutate(species_label = paste(str_split(species,"_")[[1]][1], str_split(species,"_")[[1]][2]))

text_size = 12
subtext_size = 10

variance_plot <- ggplot(variance_df, aes(species_label, variance,fill = variance_type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = subtext_size), 
        axis.text.y = element_text(size = subtext_size), 
        axis.title.y = element_text(size = text_size), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 10, l = 15)) +
  scale_fill_brewer(palette = 7) +
  labs(y = "Average variance") +
  coord_cartesian(clip = "off")

variance_plot

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/strain_freq_avg_variance.png"
ggsave(out_path, variance_plot, dpi = 300, width = 7.5, height = 6)


#########################################################################################################
# Shalon Within Host SNVs plot                                                                          #
#########################################################################################################

# remove list
rm(list = ls())

# Packages
library(dplyr)
library(cowplot)

source("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/figures/helper_scripts/SNV_plot_helper_scripts_shalon.R")


A <- plot_snvs("Bacteroides_massiliensis_44749", subject = 2, timepoints_to_use = c(2,4,5,7)) #
B <- plot_snvs("Bacteroides_vulgatus_57955", subject = 8, timepoints_to_use = c(8,9,10,11)) # 
C <- plot_snvs("Bacteroides_vulgatus_57955", subject = 11, timepoints_to_use = c(6,12,13,15))
D <- plot_snvs("Desulfovibrio_piger_61475", subject = 9, timepoints_to_use = c(2,4,5,8), sample_to_polarize_by = "SRR18585211")

A_legend <- plot_snvs("Bacteroides_massiliensis_44749",subject = 2, output = "legend")
B_legend <- plot_snvs("Bacteroides_vulgatus_57955",subject = 8, output = "legend")
C_legend <- plot_snvs("Bacteroides_vulgatus_57955",subject = 11, output = "legend")
D_legend <- plot_snvs("Desulfovibrio_piger_61475",subject = 9, output = "legend")
 

# FINAL GRID #

text_size = 12
subtext_size = 10
strain_fig_offset <- 0.0333
x_axis_title_offset <- 0.013
strain_fig_offset <- 0.0333
panel_size = 0.196
legend_offset <- 1/12




grid <- ggdraw() +
  draw_plot(A, strain_fig_offset, x_axis_title_offset+3*(1-x_axis_title_offset)/4, 1-strain_fig_offset - legend_offset, (1-x_axis_title_offset)/4) +
  draw_plot(B, strain_fig_offset, x_axis_title_offset+2*(1-x_axis_title_offset)/4, 1-strain_fig_offset - legend_offset, (1-x_axis_title_offset)/4) +
  draw_plot(C, strain_fig_offset, x_axis_title_offset+1*(1-x_axis_title_offset)/4, 1-strain_fig_offset - legend_offset, (1-x_axis_title_offset)/4) +
  draw_plot(D, strain_fig_offset, x_axis_title_offset+0*(1-x_axis_title_offset)/4, 1-strain_fig_offset - legend_offset, (1-x_axis_title_offset)/4) +
  draw_plot(A_legend, 1-legend_offset, x_axis_title_offset+3*(1-x_axis_title_offset)/4, legend_offset, (1-x_axis_title_offset)/4) +
  draw_plot(B_legend, 1-legend_offset, x_axis_title_offset+2*(1-x_axis_title_offset)/4, legend_offset, (1-x_axis_title_offset)/4) +
  draw_plot(C_legend, 1-legend_offset, x_axis_title_offset+1*(1-x_axis_title_offset)/4, legend_offset, (1-x_axis_title_offset)/4) +
  draw_plot(D_legend, 1-legend_offset, x_axis_title_offset+0*(1-x_axis_title_offset)/4, legend_offset, (1-x_axis_title_offset)/4) +
  draw_label("Device type", x = 0.4375 , y = 0.005,hjust = 0,vjust = 0.0,size = text_size,fontfamily = "Helvetica", color = "black") +
  draw_label("Frequency", x = 0.009, y = 0.5, angle = 90, size = 12, fontfamily = "Helvetica")

out_file <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/Shalon_SNVS_WithinTimepoint.png"
ggsave(out_file, grid, bg = "white", dpi = 300, width = 7.5, height = (10/6)*4, units = "in")

  
#########################################################################################################
# Family FOLD CHANGE PLOT                                                                          #
#########################################################################################################
library(dplyr)
library(compositions)

code_dir <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/HumanizedMouse/figures/figures/"
data_dir <-  "/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/"

### Helper script
source( paste0(code_dir, "helper_scripts/family_plot_helpers.R"))

study = "HumanizedMouse_Batch2"

### Load the data
data <- read.csv(paste0(data_dir,"merged_data/species/relative_abundance.txt.bz2"),sep="\t",stringsAsFactors = TRUE)
metadata <- read.csv(paste0(data_dir,"metadata/metadata.csv"),sep="\t",stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(tissue_type = ifelse(accession == "TL1gDNAshort", "Inoculum", tissue_type))
taxonomy <- read.csv(paste0(data_dir,"metadata/genome_taxonomy.txt"),sep="\t",header=TRUE)
taxInfo <- read.csv(paste0(data_dir,"metadata/genome_info.txt"),sep="\t",header = TRUE)


### Clean the data
data_sub = data[,2:ncol(data)]
rownames(data_sub) = data$species_id
range(data_sub)
species_sufficient = names(which(rowSums(data_sub > 0.001) > 0))
data_sub = data_sub[species_sufficient,]
data_sub$species_id = row.names(data_sub)
data = data_sub

### processing taxonomy dataframes

taxInfo  = taxInfo %>% filter(rep_genome > 0)
taxInfo = taxInfo %>% distinct(genome_name, species_id)

to_plot = data %>% pivot_longer(!species_id,names_to = "sample", values_to = "abundance")
to_plot$accession = to_plot$sample
to_plot <- to_plot %>% 
  mutate(accession = chartr(".", "-", accession))
intersect(to_plot$accession, metadata$accession)
to_plot = merge(to_plot, metadata, by = "accession")
nrow(to_plot)
to_plot = merge(to_plot, taxInfo, by = "species_id")
nrow(to_plot)

taxonomy_distinct = taxonomy %>% distinct(genome_name, phylum, class, order, family, genus)
to_plot = merge(to_plot,taxonomy_distinct, by = "genome_name")
to_plot$tissue_type = factor(to_plot$tissue_type,levels = c("Inoculum", "Duodenum","Jejunum","Ileum","Cecum","Colon"))

## Make stacked barplots
to_plot_ = to_plot %>% group_by(subject_id,family, order, tissue_type) %>%
  summarize(TotalAbundance = sum(abundance,na.rm=TRUE))
to_plot_stats = to_plot_ %>% group_by(order,family) %>% 
  summarize(MeanAbundance = mean(TotalAbundance,na.rm=TRUE))

#length(unique(to_plotstats$family))
length(unique(to_plot$family))

# to_plot_$subject_id = ifelse(to_plot_$subject_id != 0 ,paste0("Mouse ", to_plot_$subject_id), "Inoculum")
to_plot_$subject_id = ifelse(to_plot_$subject_id != 0 ,paste0(to_plot_$subject_id), "Inoculum")
sortie = to_plot_stats %>% arrange(desc(MeanAbundance))
to_plot_ = to_plot_ %>% arrange(order, family)
to_plot_$order[to_plot_$order == ""] = "Unidentified"
to_plot_$family[to_plot_$family == ""] = "unidentified"
to_plot_$family_order = paste0(to_plot_$order, "; ",to_plot_$family)
to_plot_$family_order = factor(to_plot_$family_order,levels = unique(to_plot_$family_order ))

family_df <- to_plot_


logFC_paired_df <- family_df %>%
  filter(subject_id != "Inoculum") %>%
  rowwise() %>%
  mutate(coarse_position = ifelse(tissue_type %in% c("Cecum", "Colon"), "Lower gut", "Upper gut")) %>%
  group_by(subject_id, family, order, family_order, coarse_position) %>%
  summarise(AverageAbundance = mean(TotalAbundance)) %>%
  tidyr::pivot_wider(names_from = coarse_position, values_from = AverageAbundance) %>%
  mutate(log2FC = log2(`Lower gut` + 1) - log2(`Upper gut` + 1),
         FC = (`Lower gut` + 1) / (`Upper gut` + 1))

blacklist_species <- logFC_paired_df %>% 
  filter(FC == 1) %>% 
  group_by(family_order) %>% 
  summarise(count = n()) %>% 
  filter(count > 1) %>% select(family_order) %>% 
  mutate(family_order = as.character(family_order)) %>% 
  pull()

summary_paired_df <- logFC_paired_df %>%
  filter(!(family_order %in% blacklist_species)) %>%
  group_by(family, order, family_order) %>%
  summarise(
    median_log2FC = median(log2FC),
    median_FC = median(FC),
    
    paired_test = list(
      wilcox.test(`Lower gut`, `Upper gut`, paired = TRUE, conf.int = TRUE, exact = FALSE)
    ),
    
    logFC_test = list(
      wilcox.test(log2FC, mu = 0, conf.int = TRUE, exact = FALSE)
    ),
    
    FC_test = list(
      wilcox.test(FC, mu = 1, conf.int = TRUE, exact = FALSE)
    )
  ) %>%
  mutate(
    p_value = paired_test[[1]]$p.value,
    
    p_value_logFC = logFC_test[[1]]$p.value,
    lower_ci_FC  = logFC_test[[1]]$conf.int[1],
    upper_ci_FC  = logFC_test[[1]]$conf.int[2],
    estimate_FC  = logFC_test[[1]]$estimate,
    
    lower_ci = FC_test[[1]]$conf.int[1],
    upper_ci = FC_test[[1]]$conf.int[2],
    estimate = FC_test[[1]]$estimate
  )



family_order_order <- summary_paired_df %>%
  arrange(desc(estimate_FC)) %>%
  select(family_order) %>%
  pull()

summary_paired_df$family_order <- factor(summary_paired_df$family_order, levels = rev(family_order_order))

p1 <- ggplot(summary_paired_df, aes(x = family_order, y = estimate_FC)) +
  geom_pointrange(aes(ymin = lower_ci_FC, ymax = upper_ci_FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_bw() + 
  labs(x = "Family; order", y = "Median within-host log2 fold change (Large vs. small intestine)") + 
  coord_flip()
p1

### SAVING
out_path_1 <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/HumanizedMouse/figures/revisions/family_log2FC_test.png"
ggsave(out_path_1, p1, dpi = 300, width = 7, height = 8)

  


#########################################################################################################
# STRAIN ANOVA PLOTS                                                                         #
#########################################################################################################

######## Humanized mice ########

# Clearing environment
rm(list=ls()) 

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
source("/u/project/ngarud/michaelw/Diversity-Along-Gut/HumanizedMouse/figures/figures/helper_scripts/strain_phasing_helper_scripts.R")

# ANOVA #

species_list <- c("Alistipes_shahii_62199",
                  "Adlercreutzia_equolifaciens_60310",
                  "Anaerostipes_hadrus_55206", 
                  "Bacteroides_uniformis_57318", 
                  "Bacteroides_vulgatus_57955", 
                  "Blautia_wexlerae_56130", 
                  "Clostridiales_bacterium_61057",
                  "Eubacterium_hallii_61477", 
                  "Ruminococcus_obeum_61472",
                  "Ruminococcus_sp_55468",
                  "Sutterella_wadsworthensis_56828")


species_subset <- c("Adlercreutzia_equolifaciens_60310",
                    "Anaerostipes_hadrus_55206",
                    "Bacteroides_uniformis_57318",
                    "Bacteroides_vulgatus_57955",
                    "Blautia_wexlerae_56130", 
                    "Clostridiales_bacterium_61057",
                    "Eubacterium_hallii_61477")

anova_results_full <- data.frame(species = NA, variable = NA, VarExpl = NA, ss = NA, df = NA, pval = NA)
lrtest_results_full <- data.frame(species = NA, p_val_region = NA, p_val_mouse = NA)


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
  
  # strain_abundance <- strain_abundance %>%
  #   filter(region %in% good_regions, mouse_number %in% good_mice, cage %in% good_cages)
  
  strain_abundance_filtered <- strain_abundance %>%
    filter(region %in% good_regions, mouse_number %in% good_mice, cage %in% good_cages)
  
  # no_cages <- length(strain_abundance %>% select(cage) %>% unique() %>% pull())
  # no_mice <- length(strain_abundance %>% select(mouse_number) %>% unique() %>% pull())
  # no_regions <- length(strain_abundance %>% select(region) %>% unique() %>% pull())
  
  no_cages <- length(strain_abundance_filtered %>% select(cage) %>% unique() %>% pull())
  no_mice <- length(strain_abundance_filtered %>% select(mouse_number) %>% unique() %>% pull())
  no_regions <- length(strain_abundance_filtered %>% select(region) %>% unique() %>% pull())
  
  
  no_cages_old <- -1
  no_mice_old <- -1
  no_regions_old <- -1
  
  while ((no_cages != no_cages_old) & (no_mice != no_cages_old) & (no_regions != no_cages_old)) {
    # no_cages_old <- no_cages
    # no_mice_old <- no_mice
    # no_regions_old <- no_regions
    # 
    # good_cages <- strain_abundance %>%
    #   filter(mouse_number != "Inoculum") %>%
    #   group_by(cage) %>%
    #   summarise(no_of_mice = n_distinct(mouse_number)) %>%
    #   rowwise() %>%
    #   mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
    #   filter(above_1) %>%
    #   select(cage) %>%
    #   pull()
    # 
    # good_mice <- strain_abundance %>%
    #   filter(mouse_number != "Inoculum") %>%
    #   group_by(mouse_number) %>%
    #   summarise(no_of_regions = n_distinct(region)) %>%
    #   rowwise() %>%
    #   mutate(above_1 = ifelse(no_of_regions > 0, TRUE, FALSE)) %>%
    #   filter(above_1) %>%
    #   select(mouse_number) %>%
    #   pull()
    # 
    # good_regions <- strain_abundance %>%
    #   filter(mouse_number != "Inoculum") %>%
    #   group_by(region) %>%
    #   summarise(no_of_mice = n_distinct(mouse_number)) %>%
    #   rowwise() %>%
    #   mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
    #   filter(above_1) %>%
    #   select(region) %>%
    #   pull()
    # 
    # strain_abundance <- strain_abundance %>%
    #   filter(region %in% good_regions, mouse_number %in% good_mice, cage %in% good_cages)
    # 
    # no_cages <- length(strain_abundance %>% select(cage) %>% unique() %>% pull())
    # no_mice <- length(strain_abundance %>% select(mouse_number) %>% unique() %>% pull())
    # no_regions <- length(strain_abundance %>% select(region) %>% unique() %>% pull())
    
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

# anova_df$variable = factor(anova_df$variable, levels = c("cage/diet","mouse", "gut segment", "residuals"))
anova_df$variable = factor(anova_df$variable, levels = c("Cage","Mouse", "Gut region", "Residuals"))



A <- ggplot(anova_df, aes(x = (species), y = VarExpl*100, fill = variable, group = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  # scale_fill_manual(values = c("Cage" = "#483248", # Eggplant
  #                              "Mouse" = "#AA98A9", #  Lilac
  #                              "Gut region" = "#BDB5D5", # Wisteria
  #                              "Residuals" = "lightgrey")) +
  scale_fill_manual(values = c("Cage" = "#AA98A9", 
                               "Mouse" = "#BDB5D5", 
                               "Gut region" = "#CC0000", 
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
  labs(x = "Species", y = "Variance explained (%)", fill = "Variable")


A

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/hm_strain_anova.png"
ggsave(out_path, A, dpi = 300, height = 5, width = 7.5)

######## conventional mice ######## 

# Packages
library(dplyr)
library(lmtest)
library(tidyr)
library(ggplot2)

source("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/helper_scripts/family_plot_helpers.R")


text_size = 12
subtext_size = 10


# Load species list
strain_phasing_dir <- "/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/strain_phasing/strain_clusters/"
species_list <- list.dirs(strain_phasing_dir, full.names = FALSE)[-1]

# constructing ANOVA table
anova_results_full <- data.frame(species = NA, variable = NA, VarExpl = NA, ss = NA, df = NA, pval = NA)


for (species in species_list) {
  print(paste0("Processing ", species))
  # DATA
  ## Paths
  project_folder <- "/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/"
  strain_abundance_path <- paste0(project_folder, "strain_phasing/strain_clusters/", species, "/", species, "_strain_frequency.csv")
  ## Loading the data
  strain_abundance <- read.csv2(strain_abundance_path, sep = "\t") %>%
    mutate(freq = as.numeric(freq)) %>%
    mutate(freq_trans = as.vector(clr(freq)))
  # mutate(freq_trans = asin(sqrt(freq)))
  
  # rename mice
  strain_abundance <- strain_abundance %>%
    rowwise() %>%
    mutate(mouse_number = paste0("Mouse ",as.character(extract_mouse(`sample`))))
  
  
  ## Select the major strain
  
  major_strain <- strain_abundance %>% 
    group_by(strain) %>%
    summarise(mean_freq = mean(freq)) %>%
    filter(mean_freq == max(mean_freq)) %>%
    pull(strain)
  
  # ANOVA
  
  good_cages <- strain_abundance %>%
    group_by(cage) %>%
    summarise(no_of_mice = n_distinct(mouse_number)) %>%
    rowwise() %>%
    mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
    filter(above_1) %>%
    select(cage) %>%
    pull()
  
  good_mice <- strain_abundance %>%
    group_by(mouse_number) %>%
    summarise(no_of_regions = n_distinct(region)) %>%
    rowwise() %>%
    mutate(above_1 = ifelse(no_of_regions > 0, TRUE, FALSE)) %>%
    # mutate(above_1 = ifelse(no_of_regions > 1, TRUE, FALSE)) %>%
    filter(above_1) %>%
    select(mouse_number) %>%
    pull()
  
  good_regions <- strain_abundance %>%
    group_by(region) %>%
    summarise(no_of_mice = n_distinct(mouse_number)) %>%
    rowwise() %>%
    mutate(above_1 = ifelse(no_of_mice > 1, TRUE, FALSE)) %>%
    filter(above_1) %>%
    select(region) %>%
    pull()
  
  strain_abundance_filtered <- strain_abundance %>%
    filter(region %in% good_regions, mouse_number %in% good_mice, cage %in% good_cages)
  
  # no_cages <- length(strain_abundance %>% select(cage) %>% unique() %>% pull())
  # no_mice <- length(strain_abundance %>% select(mouse_number) %>% unique() %>% pull())
  # no_regions <- length(strain_abundance %>% select(region) %>% unique() %>% pull())
  
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
  
  strain_abundance <- strain_abundance %>%
    filter(strain == major_strain)
  
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
  mutate(variable = ifelse(variable == "mouse_number", "Mouse", variable)) %>%
  mutate(variable = ifelse(variable == "region", "Gut region", variable)) %>%
  mutate(variable = ifelse(variable == "cage", "Cage", variable))

anova_df$variable = factor(anova_df$variable, levels = c("Cage","Mouse", "Gut region", "Residuals"))

species_metadata_path <- "/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc/metadata.tsv"
species_metadata <- read.csv2(species_metadata_path, sep = "\t")

anova_df <- anova_df %>%
  rowwise() %>%
  mutate(species_name = extract_species(species_metadata[species_metadata["species_id"] == species, "species"]),
         species_name = ifelse(species == "213583", "Anaerotruncus sp.", species_name))


B <- ggplot(anova_df, aes(x = (species_name), y = VarExpl*100, fill = variable, group = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  # scale_fill_manual(values = c("Cage" = "#AA98A9", # Eggplant
  #                              "Mouse" = "#483248", #  Lilac
  #                              "Gut region" = "#BDB5D5", # Wisteria
  #                              "Residuals" = "lightgrey")) +
  # scale_fill_manual(values = c("Cage" = "#fdbf6f", 
  #                              "Mouse" = "#ff7f00", 
  #                              "Gut region" = "#e31a1c", #"#CC0000", # Wisteria
  #                              "Residuals" = "lightgrey")) +
  scale_fill_manual(values = c("Cage" = "#AA98A9", 
                               "Mouse" = "#BDB5D5", 
                               "Gut region" = "#CC0000", 
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
  labs(x = "Species", y = "Variance explained (%)", fill = "Variable")
# B


out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/cm_strain_anova.png"
ggsave(out_path, B, dpi = 300, height = 5, width = 7.5)



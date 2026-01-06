######################################################################################
# STRAIN PHASING #####################################################################
######################################################################################

# Clearing environment
rm(list=ls()) 

# Packages
library(dplyr)
library(ggplot2)
library(cowplot)

# Helper scripts
source("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/figures/helper_scripts/strain_phasing_helper_scripts.R")

# Plotting strains ####################################################################

species_to_plot <- c("100320", 
                     "100555", 
                     "100600", 
                     "207693",
                     "207693",
                     "213583",
                     "214603",
                     "217240",
                     "229547",
                     "253991",
                     "261672",
                     "263672")

for (species_name in species_to_plot){
  # Plotting
  strain_plot <- plot_strains(species_name)
  # Saving
  out_file <- paste0("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/figures/strain_phasing/", species_name, "_strain_plot.png")
  ggsave(out_file, strain_plot, bg = "white", dpi = 300, width = 7.5, height = 2.5, units = "in")
  
}

strain_plot_100320 <- plot_strains("100320")
strain_plot_100555 <- plot_strains("100555")
strain_plot_100600 <- plot_strains("100600")
strain_plot_207693 <- plot_strains("207693")
strain_plot_213583 <- plot_strains("213583")
strain_plot_214603 <- plot_strains("214603")
strain_plot_217240 <- plot_strains("217240")
strain_plot_253991 <- plot_strains("253991")
strain_plot_261672 <- plot_strains("261672")
strain_plot_263672 <- plot_strains("263672")

strain_plot_207693 <- plot_strains("207693", output = "normal")
strain_plot_207693

out_file <- paste0("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/figures/strain_phasing/", "207693", "_QCbio.png")
ggsave(out_file, strain_plot_207693, bg = "white", dpi = 300, width = 7.5, height = 2.5, units = "in")


# return to this species
strain_plot_229547 <- plot_strains("229547")


# ANOVA ####################################################################
# Clearing environment
rm(list=ls()) 
# Packages
library(dplyr)
library(ggplot2)
library(cowplot)

# parameters

text_size = 12
subtext_size = 10

# functions
source("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/helper_scripts/family_plot_helpers.R")



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

anova_df <- anova_results_full %>%
  filter(species %in% species_subset) %>%
  rowwise() %>%
  mutate(variable = ifelse(variable == "mouse_number", "Mouse", variable)) %>%
  mutate(variable = ifelse(variable == "region", "Gut region", variable)) %>%
  mutate(variable = ifelse(variable == "cage", "Cage", variable))

anova_df$variable = factor(anova_df$variable, levels = c("Cage","Mouse", "Gut region", "Residuals"))

# label with taxonomy

species_metadata_path <- "/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc/metadata.tsv"
species_metadata <- read.csv2(species_metadata_path, sep = "\t")

anova_df <- anova_df %>%
  rowwise() %>%
  mutate(species = extract_species(species_metadata[species_metadata["species_id"] == species, "species"]))



# plotting
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


A

out_file <- paste0("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/figures/strain_phasing/anova.png")
ggsave(out_file, A, bg = "white", dpi = 300, width = 7, height = 5, units = "in")


# SCRATCH #########################################################################

species = "100320"
output = "normal"
strain_colors = c("#654321", "#D2B48C")

print(paste0("Processing ", species, "."))

# Loading data

strain_trajectory_path <- paste0("/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/strain_phasing/strain_clusters/", 
                                 species, 
                                 "/",
                                 species, 
                                 "_strain_frequency.csv")
print(paste0("Processing ", species, "."))

# Setting parameters

font_size <- 12
subtext_size <- 10
x_axis_font_size <- 8

# Loading data

strain_trajectory_path <- paste0("/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/strain_phasing/strain_clusters/", 
                                 species, 
                                 "/",
                                 species, 
                                 "_strain_frequency.csv")


#Loading strain frequencies
strain_trajectory_df <- read.csv2(strain_trajectory_path, sep = "\t")
strain_trajectory_df$freq <- as.numeric(strain_trajectory_df$freq)
strain_trajectory_df$quantile_25 <- as.numeric(strain_trajectory_df$quantile_25)
strain_trajectory_df$quantile_75 <- as.numeric(strain_trajectory_df$quantile_75)
strain_trajectory_df$upper_ci <- as.numeric(strain_trajectory_df$upper_ci)
strain_trajectory_df$lower_ci <- as.numeric(strain_trajectory_df$lower_ci)
strain_trajectory_df$species <- as.character(strain_trajectory_df$species)
strain_trajectory_df$strain <- as.character(strain_trajectory_df$strain)
strain_trajectory_df <- strain_trajectory_df %>%
  mutate(mouse_id = gsub(" ", "_", paste(cage, mouse_number)))

# Getting species names
species_metadata_path <- "/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc/metadata.tsv"
species_metadata_df <- read.csv2(species_metadata_path, sep = "\t") 
species_metadata_df$species_id <- as.character(species_metadata_df$species_id) 

taxonomy <- str_split(species_metadata_df[species_metadata_df['species_id'] == species,"species"], "__")[[1]][2]


## Changing mouse number

change_mouse_number <- function(mouse_id) {
  cage <- as.numeric(str_split(mouse_id, "_")[[1]][2])
  mouse <- as.numeric(str_split(mouse_id, "_")[[1]][4])
  new_mouse_number = (cage - 1)*3 + mouse
  mouse_label = paste0("Mouse ", as.character(new_mouse_number))
  return(mouse_label)
}

strain_trajectory_df$mouse_label <- sapply(strain_trajectory_df$mouse_id, change_mouse_number)

## Identify strain order

mice_present <- as.vector(strain_trajectory_df %>% filter(!is.na(freq)) %>% select(mouse_id) %>% unique() %>% pull())

strain_order <- strain_trajectory_df %>% 
  group_by(mouse_id, strain) %>% 
  summarise(avg_freq = mean(freq, na.rm = TRUE)) %>% 
  group_by(strain) %>% 
  summarise(avg_strain_freq = mean(avg_freq)) %>% 
  arrange(desc(avg_strain_freq)) %>% 
  select(strain) %>% 
  pull() 

#Shortening gut labels
shorten_gut_label <- function(label) {
  if (is.na(label)) {
    return("Inoculum")
  } else if ((label == "Duodenum") | (label == "Jejunum") | (label == "Ileum")) {
    return(substr(label, 1, 1))
  } else if ((label == "Cecum") | (label == "Colon")) {
    return(substr(label, 1, 2))
  } else if (label == "Inoculum") {
    return(label)
  }
}

strain_trajectory_df$gut_label <- sapply(strain_trajectory_df$region, shorten_gut_label)

# Add empty rows into mice in which the strain is absent

if (!("Cage_1_Mouse_1" %in% strain_trajectory_df$mouse_id)) {
  strain_trajectory_df <- strain_trajectory_df %>%
    add_row(species = species, strain = "1", mouse_number = "Mouse 1", region = "Colon", cage = "Cage 1", sample = "Col_1_1", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_1_Mouse_1", mouse_label = "Mouse 1", gut_label = "Co")
}
if (!("Cage_1_Mouse_2" %in% strain_trajectory_df$mouse_number)) {
  strain_trajectory_df <- strain_trajectory_df %>%
    add_row(species = species, strain = "1", mouse_number = "Mouse 2", region = "Colon", cage = "Cage 1", sample = "Col_1_2", freq = 0, quantile_25 = NA, quantile_75 = NA,upper_ci = NA, lower_ci = NA, mouse_id = "Cage_1_Mouse_2", mouse_label = "Mouse 2", gut_label = "Co")
}
if (!("Cage_1_Mouse_3" %in% strain_trajectory_df$mouse_number)) {
  strain_trajectory_df <- strain_trajectory_df %>%
    add_row(species = species, strain = "1", mouse_number = "Mouse 3", region = "Colon", cage = "Cage 1", sample = "Col_1_3", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_1_Mouse_3", mouse_label = "Mouse 3", gut_label = "Co")
}
if (!("Cage_2_Mouse_1" %in% strain_trajectory_df$mouse_number)) {
  strain_trajectory_df <- strain_trajectory_df %>%
    add_row(species = species, strain = "1", mouse_number = "Mouse 1", region = "Colon", cage = "Cage 2", sample = "Col_2_1", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_2_Mouse_1", mouse_label = "Mouse 4", gut_label = "Co")
}
if (!("Cage_2_Mouse_2" %in% strain_trajectory_df$mouse_number)) {
  strain_trajectory_df <- strain_trajectory_df %>%
    add_row(species = species, strain = "1", mouse_number = "Mouse 2", region = "Colon", cage = "Cage 2", sample = "Col_2_2", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_2_Mouse_2", mouse_label = "Mouse 5", gut_label = "Co")
}
if (!("Cage_2_Mouse_3" %in% strain_trajectory_df$mouse_number)) {
  strain_trajectory_df <- strain_trajectory_df %>%
    add_row(species = species, strain = "1", mouse_number = "Mouse 3", region = "Colon", cage = "Cage 2", sample = "Col_2_3", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_2_Mouse_3", mouse_label = "Mouse 6", gut_label = "Co")
}
if (!("Cage_3_Mouse_1" %in% strain_trajectory_df$mouse_number)) {
  strain_trajectory_df <- strain_trajectory_df %>%
    add_row(species = species, strain = "1", mouse_number = "Mouse 1", region = "Colon", cage = "Cage 3", sample = "Col_3_1", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA,  mouse_id = "Cage_3_Mouse_1", mouse_label = "Mouse 7", gut_label = "Co")
}
if (!("Cage_3_Mouse_2" %in% strain_trajectory_df$mouse_number)) {
  strain_trajectory_df <- strain_trajectory_df %>%
    add_row(species = species, strain = "1", mouse_number = "Mouse 2", region = "Colon", cage = "Cage 3", sample = "Col_3_2", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_3_Mouse_2", mouse_label = "Mouse 8", gut_label = "Co")
}
if (!("Cage_3_Mouse_3" %in% strain_trajectory_df$mouse_number)) {
  strain_trajectory_df <- strain_trajectory_df %>%
    add_row(species = species, strain = "1",  mouse_number = "Mouse 3", region = "Colon", cage = "Cage 3", sample = "Col_3_3", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_3_Mouse_3", mouse_label = "Mouse 9", gut_label = "Co")
}

#Factoring
gut_order = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")
gut_label_order = c("D", "J", "I", "Ce", "Co")
mouse_order = c("Cage_1_Mouse_1", "Cage_1_Mouse_2", "Cage_1_Mouse_3", "Cage_2_Mouse_1", "Cage_2_Mouse_2", "Cage_2_Mouse_3", "Cage_3_Mouse_1", "Cage_3_Mouse_2", "Cage_3_Mouse_3")
strain_trajectory_df$region = factor(strain_trajectory_df$region, levels = gut_order)
strain_trajectory_df$gut_label = factor(strain_trajectory_df$gut_label, levels = gut_label_order)
strain_trajectory_df$mouse_id = factor(strain_trajectory_df$mouse_id, levels = mouse_order)
strain_trajectory_df$strain = factor(strain_trajectory_df$strain, levels = strain_order)


# Palette
if (length(strain_order) == 3) {
  pal <- c("#404040","#808080","#C0C0C0" )
} else if (length(strain_order) == 2) {
  pal <- strain_colors
} else if (length(strain_order) == 1) {
  pal <- c("#654321") #brown
}


error_bars <- strain_trajectory_df %>% 
  dplyr::filter(strain == strain_order[1])


strain_order <- strain_trajectory_df %>% 
  group_by(mouse_id, strain) %>% 
  summarise(avg_freq = mean(freq, na.rm = TRUE)) %>% 
  group_by(strain) %>% 
  summarise(avg_strain_freq = mean(avg_freq)) %>% 
  arrange(avg_strain_freq) %>% 
  select(strain) %>% 
  pull()

strain_trajectory_df$strain = factor(strain_trajectory_df$strain, levels = strain_order)

cage_3 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 3"), aes(x = gut_label, y = freq, fill = strain)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_x_discrete(limits = gut_label_order) +
  facet_wrap("mouse_number", scales = "free_x") +
  theme(text = element_text(family = "Helvetica", size = x_axis_font_size),
        axis.title = element_blank(), 
        axis.text.x = element_text(size = x_axis_font_size),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.background = element_rect(fill = scales::alpha("#08519c", 0.85), color = "black"),
        strip.text = element_text(color = "white", size = 10), 
        legend.title = element_text(size = x_axis_font_size),
        legend.text = element_text(size = x_axis_font_size),
        plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef")) +
  scale_fill_manual(values = pal) +
  guides(fill = guide_legend(override.aes = list(size = 3))) 


legend <- get_legend(cage_3)
as_ggplot(legend)

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
code_dir <- "~/Wasney-Briscoe/scripts/figures/"
metadata_dir <- "~/Wasney-Briscoe/metadata/"
data_dir <-  "~/"

source(paste0(code_dir,"helper_scripts/humanized_mouse_utilities.R"))
source(paste0(code_dir, "helper_scripts/mouse_annotation_functions.R"))

# Loading the data
## metaata
metadata = read.csv(paste0(metadata_dir,"metadata.csv"),sep="\t",stringsAsFactors = FALSE)

## Data

data_single <- read.csv(paste0(data_dir, "popgen_stats/SinglePi_SchloissnigPi_cov4_SiteDownsampled.csv")) #%>%


# Filtering the data
data_single <- data_single %>%
  mutate(total_loci = as.numeric(total_loci)) %>%
  # filter(good_species == TRUE) %>%
  filter(total_loci > 5e5) %>%
  rowwise() %>%
  mutate(species = if_else(
    str_detect(species, "Ruminococcus_sp") | str_detect(species, "Faecalibacterium_prausnitzii"),
    paste(str_split(gsub("_", " ", species), " ")[[1]][1:3], collapse = " "),  # Keep 3 elements for "Ruminococcus"
    paste(str_split(gsub("_", " ", species), " ")[[1]][1:2], collapse = " ")   # Keep 2 elements otherwise
  ))

#Choosing good species

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


# Calculate species order
species_order <- data_single %>%
  filter(sample == "TL1gDNAshort") %>%
  arrange(desc(genomewide_pi))  %>%
  select(species) %>%
  pull()

# Merging the data
data_single <- data_single %>% 
  mutate(pi_type = ifelse(sample == "TL1gDNAshort","Inoculum pi","Intra-sample pi"))

all_pi_data <- data_single
# Summary stats
average_labels <- all_pi_data %>% group_by(species,pi_type) %>% 
  summarise(mean_pi = mean(genomewide_pi)) %>%
  rowwise() 
all_pi_data <- all_pi_data %>%
  rowwise() %>%
  mutate(pi_type = ifelse(pi_type == "Intra-sample pi", "Intra-sample pi (average)", pi_type))


# annotating/filtering data

all_pi_data <- all_pi_data %>%
  filter(pi_type == "Intra-sample pi (average)", sample != "TL1gDNAshort") %>% 
  rowwise() %>%
  mutate(gut_segment = gut_site)

#Plotting

all_pi_data$species <- factor(x = all_pi_data$species, levels = species_order)
average_labels$species <- factor(x = average_labels$species, levels = species_order)
levels = c("Inoculum pi",
           "Intra-sample pi (average)")
average_labels$pi_type <- factor(x = average_labels$pi_type, levels = levels)

average_labels <- average_labels %>%
  filter(pi_type %in% c("Intra-sample pi (average)", 
                        "Inoculum pi"))

all_pi_data$gut_segment <- factor(all_pi_data$gut_segment, levels = c("Duodenum",  "Jejunum", "Ileum", "Cecum", "Colon", "Stool"))

ten_highest_and_lowest_species <- c(species_order %>% head(10), species_order %>% tail(10))



species_with_multiple_strains <- average_labels %>%
  filter(mean_pi > 1e-3) %>%
  select(species) %>%
  pull()


font_face <- ifelse(levels(all_pi_data$species) %in% species_with_multiple_strains, "#CD5C5C", "black")

A <- ggplot(all_pi_data, aes(x = species, y =genomewide_pi)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(size = 2,alpha = 0.5, aes(color = mouse), position = position_dodge(width = 0.75)) +
  scale_color_manual(name = "Mouse sample \u03C0",
                     values = c(
                       '1'='#4B0082',
                       '2'='#000000',
                       '3'='#87CEEB',
                       '4'='#FF4500',
                       '5'='#FFD700',
                       '6'='#006400',
                       '7'='#FF69B4',
                       '8'='#1E90FF'
                     )) +
  new_scale_color() +
  geom_text(data = average_labels %>% arrange(pi_type) %>% filter(pi_type == "Inoculum pi"),
            label = "*",
            size = text_size,
            vjust = 0.75,
            hjust = "middle",
            aes(x = species,y = mean_pi, color = pi_type)) +
  scale_color_manual(name = "Inoculum \u03C0",
                     labels = c(""),
                     values = c("brown") #"#33A02C", "#1F78B4",
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = text_size),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(family = "Helvetica", size = subtext_size,
                               angle = 45,
                               colour = font_face,
                               hjust = 1),
    axis.title = element_text(family = "Helvetica", size = subtext_size),  # Adjust size to match the first plot
    legend.text = element_text(family = "Helvetica", size = subtext_size),  # Adjust size to match the first plot
    legend.title = element_text(family = "Helvetica", size = subtext_size),
    # plot.margin = margin(5, 5, 5, 30),
    plot.margin = margin(5, 5, 5, 30),
    axis.title.x = element_blank(),
    legend.margin = margin(0, 0, 0, -2)
  ) +
  scale_y_log10(labels = function(x) parse(text = sprintf('10^{%g}', log10(x)))) +
  labs(y = "Genomewide \u03C0") +
  guides(color = guide_legend(label.wrap = 2))

out_file <- "~/figure_3.png"
ggsave(out_file, A, height = 5, width = 7.5, units = "in", bg = "white", dpi = 300)





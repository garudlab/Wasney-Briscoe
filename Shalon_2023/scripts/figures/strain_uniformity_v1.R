library(dplyr)
library(ggplot2)
library(cowplot)
library(ggtext)



# functions

source("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/figures/helper_scripts/strain_phasing_helper_scripts.R")


############################################################################
# LOOPING THROUGH SPECIES-SUBJECTS                                         #
############################################################################

strain_species_subject_df_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/metadata/strain_species_subject_df.txt"
strain_species_subject_df <- read.csv2(strain_species_subject_df_path, sep = "\t")

for (i in 1:nrow(strain_species_subject_df)) {
  species_ <- strain_species_subject_df[i,"species"]
  subject_ <- strain_species_subject_df[i,"subject"]
  strain_plot <- plot_strains(species = species_, subject = subject_)
  out_path <- paste0("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/figures/strain_phasing/strain_phasing_barplots/", species_, "_subject_", as.character(subject_), "_barplot.png")
  ggsave(out_path, strain_plot, bg = "white", dpi = 300, width = 6, height = 3, units = "in")
}


############################################################################
# FIGURE                                                                   #
############################################################################
## Strains
### Temporal
b_vulgatus_s2 <- plot_strains(species = "Bacteroides_vulgatus_57955", subject = 2, timepoints_to_use = c(3,5,6,7))
R_obeum_s11 <- plot_strains(species = "Ruminococcus_obeum_61472", subject = 11, timepoints_to_use = c(3,5,6,8))
### spatial
P_distasonis_s9 <- plot_strains(species = "Parabacteroides_distasonis_56985", subject = 9, timepoints_to_use = c(1,2,3,5)) 
b_vulgatus_s11 <- plot_strains(species = "Bacteroides_vulgatus_57955", subject = 11, timepoints_to_use = c(3,6,9,10))


## ANOVA
### Load dataframe
strain_anova_df_path <- "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/strain_phasing/strain_anova.csv"
strain_anova_df <- read.csv2(strain_anova_df_path, sep = ",") %>%
  mutate(value = as.numeric(value))

### Create species-subject group

strain_anova_df <- strain_anova_df %>%
  rowwise() %>%
  mutate(species_subject_group = paste0(species, "-", as.character(subject)))


### factor

strain_anova_df$variance_type <- factor(strain_anova_df$variance_type, levels = c("temporal effect", "spatial effect", "residual variance"))

### Plot

font_size = 12
x_axis_font_size = 8

species_A = "Bacteroides_vulgatus_57955"
subject_A = 2
color_A = "red"
species_B = "Ruminococcus_obeum_61472"
subject_B = 11
color_B = "purple"
species_C = "Parabacteroides_distasonis_56985"
subject_C = 9
color_C = "blue"
species_D = "Bacteroides_vulgatus_57955"
subject_D = 11
color_D = "darkgreen"

B <- ggplot(strain_anova_df, aes(x = variance_type, y = value)) +
  geom_boxplot(alpha = 0.5) + #aes(fill = variance_type), 
  geom_point(size = 2) +
  geom_point(
    data = subset(strain_anova_df, species == species_A & subject == subject_A),
    size = 2,
    color = color_A
  ) +
  geom_point(
    data = subset(strain_anova_df, species == species_B & subject == subject_B),
    size = 2,
    color = color_B
  ) +
  geom_point(
    data = subset(strain_anova_df, species == species_C & subject == subject_C),
    size = 2,
    color = color_C
  ) +
  geom_point(
    data = subset(strain_anova_df, species == species_D & subject == subject_D),
    size = 2,
    color = color_D
  ) +
  geom_line(aes(group = species_subject_group)) +
  geom_line(
    data = subset(strain_anova_df, species == species_A & subject == subject_A),
    aes(group = species_subject_group),
    color = color_A,
    size = 1
  ) +
  geom_line(
    data = subset(strain_anova_df, species == species_B & subject == subject_B),
    aes(group = species_subject_group),
    color = color_B,
    size = 1
  ) +
  geom_line(
    data = subset(strain_anova_df, species == species_C & subject == subject_C),
    aes(group = species_subject_group),
    color = color_C,
    size = 1
  ) +
  geom_line(
    data = subset(strain_anova_df, species == species_D & subject == subject_D),
    aes(group = species_subject_group),
    color = color_D,
    size = 1
  ) +
  geom_text(
    data = subset(strain_anova_df, 
                  species == species_A & subject == subject_A & variance_type == "temporal effect"),
    aes(label = "A",fontface = "bold",),
    vjust = 1, 
    hjust = 1.5,
    family = "Helvetica",
    size = x_axis_font_size,
    color = color_A
  ) + 
  geom_text(
    data = subset(strain_anova_df, 
                  species == species_B & subject == subject_B & variance_type == "temporal effect"),
    aes(label = "B",fontface = "bold",),
    vjust = 0.5, 
    hjust = 1.5,
    family = "Helvetica",
    size = x_axis_font_size,
    color = color_B
  ) + 
  geom_text(
    data = subset(strain_anova_df, 
                  species == species_C & subject == subject_C & variance_type == "temporal effect"),
    aes(label = "C",fontface = "bold",),
    vjust = 0, 
    hjust = 1.5,
    family = "Helvetica",
    size = x_axis_font_size,
    color = color_C
  ) + 
  geom_text(
    data = subset(strain_anova_df, 
                  species == species_D & subject == subject_D & variance_type == "temporal effect"),
    aes(label = "D",fontface = "bold",),
    vjust = 1, 
    hjust = 1.5,
    family = "Helvetica",
    size = x_axis_font_size,
    color = color_D
  ) + 
  theme_bw() + 
  theme(text = element_text(family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = font_size),
        axis.text.x = element_text(size = x_axis_font_size),
        axis.text.y = element_text(size = x_axis_font_size),
        legend.position = "none"
  ) +
  scale_fill_manual(values = c("red", "black", "black")) +
  labs(y = "Variance explained (%)")
# B

## Grid

strain_fig_offset <- 0.015

grid <- ggdraw() +
  draw_plot(b_vulgatus_s2, strain_fig_offset <- 0.015, 4/5,1-strain_fig_offset,1/5) +
  draw_plot(R_obeum_s11, strain_fig_offset, 3/5,1-strain_fig_offset,1/5) +
  draw_plot(P_distasonis_s9, strain_fig_offset, 2/5,1-strain_fig_offset,1/5) +
  draw_plot(b_vulgatus_s11, strain_fig_offset, 1/5,1-strain_fig_offset,1/5) +
  draw_plot(B, strain_fig_offset, 0/5,1-strain_fig_offset,1/5) +
  draw_plot_label(c("A", "B", "C", "D", "E"), c(-0.01, -0.01, -0.01, -0.01, -0.01), c(1, 4/5, 3/5, 2/5, 1/5), size = 16, family = "Helvetica")
# grid

out_file <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/figures/revisions/strain_uniformity_V2.png"
ggsave(out_file, grid, bg = "white", dpi = 300, width = 7.5, height = 10, units = "in")


############################################################################
# FIGURE - OLD                                                             #
############################################################################

# Plot 1: with ANOVA

## species plots

b_vulgatus_s2 <- plot_strains(species = "Bacteroides_vulgatus_57955", subject = 2)
b_vulgatus_s11 <- plot_strains(species = "Bacteroides_vulgatus_57955", subject = 11)
b_longum_s5 <- plot_strains(species = "Bifidobacterium_longum_57796", subject = 5)
b_wadsworthia_s9 <- plot_strains(species = "Bilophila_wadsworthia_57364", subject = 9)
B_bacterium_s6 <- plot_strains(species = "Burkholderiales_bacterium_56577", subject = 6) 
C_mitsuokai_s10 <- plot_strains(species = "Catenibacterium_mitsuokai_61547", subject = 10)
P_distasonis_s9 <- plot_strains(species = "Parabacteroides_distasonis_56985", subject = 9) 
R_obeum_s11 <- plot_strains(species = "Ruminococcus_obeum_61472", subject = 11)
A_intestini_s9 <- plot_strains(species = "Acidaminococcus_intestini_54097", subject = 9)
E_rectale_9 <- plot_strains(species = "Eubacterium_rectale_56927", subject = 9, use_simple_timepoint = TRUE)

# b_vulgatus_s2
# b_vulgatus_s11
# b_longum_s5
# b_wadsworthia_s9
# B_bacterium_s6
# C_mitsuokai_s10
# P_distasonis_s9
# R_obeum_s11
# A_intestini_s9
# E_rectale_9

## ANOVA

### Load dataframe
strain_anova_df_path <- "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/strain_phasing/strain_anova.csv"
strain_anova_df <- read.csv2(strain_anova_df_path, sep = ",") %>%
  mutate(value = as.numeric(value))

### Create species-subject group

strain_anova_df <- strain_anova_df %>%
  rowwise() %>%
  mutate(species_subject_group = paste0(species, "-", as.character(subject)))


### factor

strain_anova_df$variance_type <- factor(strain_anova_df$variance_type, levels = c("temporal effect", "spatial effect", "residual variance"))

### Plot

font_size = 12
x_axis_font_size = 8

B <- ggplot(strain_anova_df, aes(x = variance_type, y = value)) +
  geom_boxplot(aes(fill = variance_type), alpha = 0.5) + 
  geom_point(size = 2) +
  geom_line(aes(group = species_subject_group)) +
  theme_bw() + 
  theme(text = element_text(family = "Helvetica"),
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size),
        axis.text.x = element_text(size = x_axis_font_size),
        axis.text.y = element_text(size = x_axis_font_size),
        legend.position = "none"
        ) +
  scale_fill_manual(values = c("red", "black", "black")) +
  labs(x = "Variable", y = "Variance explained (%)")

  
## Grid

strain_fig_offset <- 0.015

grid <- ggdraw() +
  draw_plot(b_vulgatus_s11, strain_fig_offset <- 0.015, 4/5,1-strain_fig_offset,1/5) +
  draw_plot(B_bacterium_s6, strain_fig_offset, 3/5,1-strain_fig_offset,1/5) +
  draw_plot(R_obeum_s11, strain_fig_offset, 2/5,1-strain_fig_offset,1/5) +
  draw_plot(b_longum_s5, strain_fig_offset, 1/5,1-strain_fig_offset,1/5) +
  draw_plot(B, strain_fig_offset, 0/5,1-strain_fig_offset,1/5) +
  draw_plot_label(c("A", "B"), c(-0.01, -0.01), c(1, 1/5), size = 16, family = "Helvetica")
grid

out_file <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/figures/revisions/strain_uniformity.png"
ggsave(out_file, grid, bg = "white", dpi = 300, width = 7.5, height = 10, units = "in")

############################ NEW SCRATCH ############################
species = "Parabacteroides_distasonis_56985"
subject = 9
pal =  c("#654321", "#D2B48C")
path = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/strain_phasing/strain_freqs.csv"
timepoints_to_use = c()
font_size <- 12
x_axis_font_size <- 8


# reassigning variables
species_oi <- species
subject_oi <- subject
strain_pal <- pal

# Loading strain frequency df
strain_freq_df <- read.csv2(path, sep = ",") %>%
  mutate(freq = as.numeric(freq), 
         lower_ci = as.numeric(lower_ci),
         upper_ci = as.numeric(upper_ci)) 

strain_freq_df <- strain_freq_df %>%
  filter(species == species_oi, subject == subject_oi) 

# Switch strain order such that strain 1 is the "major" strain
## switch strains if need be
strain_freq_df <- strain_freq_df %>%
  group_by(strain) %>%
  mutate(mean_freq = mean(freq)) %>%  
  ungroup() %>%
  distinct(strain, mean_freq) %>%     
  arrange(desc(mean_freq)) %>%
  mutate(new_strain_label = row_number()) %>%  
  select(strain, new_strain_label) %>%
  right_join(strain_freq_df, by = "strain") 

# Label with timepoint
strain_freq_df <- strain_freq_df %>%
  distinct(date, time) %>%                    
  arrange(date, time) %>%
  mutate(timepoint = paste0("timepoint ", row_number())) %>%
  right_join(strain_freq_df, by = c("date", "time")) %>%
  arrange(date, time)

# Subsetting timepoints if necessary

if (length(timepoints_to_use) != 0) {
  strain_freq_df <- strain_freq_df %>%
    filter(as.numeric(str_remove(timepoint, "timepoint ")) %in% timepoints_to_use)
  
  strain_freq_df <- strain_freq_df %>%
    select(-timepoint) %>%
    distinct(date, time) %>%                    
    arrange(date, time) %>%
    mutate(timepoint = paste0("timepoint ", row_number())) %>%
    right_join(strain_freq_df %>% select(-timepoint), by = c("date", "time")) %>%
    arrange(date, time)
}



## splitting duplicate names

strain_freq_df <- strain_freq_df %>%
  group_by(timepoint, sample_type_number, strain) %>%
  mutate(dup_index = row_number(),
         no_of_duplicates = n()) %>%
  ungroup() %>%
  mutate(
    sample_type_number = if_else(
      no_of_duplicates == 1,
      as.character(sample_type_number),
      paste0(sample_type_number, " (", dup_index, ")")
    )
  ) %>%
  select(-dup_index, -no_of_duplicates)  # optional: drop helper column


## add capsules that are missing

timepoint_groups <- split(strain_freq_df, list(strain_freq_df$date, strain_freq_df$time, strain_freq_df$timepoint), drop = TRUE)
full_capsule_list <- c("Capsule 1", 
                       "Capsule 2", 
                       "Capsule 3", 
                       "Capsule 4")


for (timepoint_name in names(timepoint_groups)) {
  timepoint_group_subset <- timepoint_groups[[timepoint_name]]
  
  # get sample type list
  sample_type_list <- timepoint_group_subset$sample_type_number %>% unique()
  sample_type_list <- gsub("^Stool \\(\\d+\\)$", "Stool", sample_type_list)
  
  if (("Stool" %in% sample_type_list)) {
    
    next
    
  } else {
    
    # getting date, time and timepoint
    parts <- strsplit(timepoint_name, "\\.")[[1]]
    date_value <- parts[1]
    time_value <- parts[2]
    timepoint_value <- parts[3]
    
    # Identifying missing capsules
    capsules_present <- timepoint_group_subset$sample_type_number %>% unique()
    missing_capsules <- setdiff(full_capsule_list, capsules_present)
    
    for (capsule in missing_capsules) {
      new_capsule_row_1 <- c(date_value, time_value, timepoint_value, 1, 1, "Capsule", NA, as.numeric(0), species_oi, subject_oi, NA, NA, NA, NA, NA, capsule)
      new_capsule_row_1 <- as.data.frame(t(new_capsule_row_1), stringsAsFactors = FALSE)
      colnames(new_capsule_row_1) <- colnames(strain_freq_df)
      
      new_capsule_row_2 <- c(date_value, time_value, timepoint_value, 2, 2, "Capsule", NA, as.numeric(0), species_oi, subject_oi, NA, NA, NA, NA, NA, capsule)
      new_capsule_row_2 <- as.data.frame(t(new_capsule_row_2), stringsAsFactors = FALSE)
      colnames(new_capsule_row_2) <- colnames(strain_freq_df)
      
      strain_freq_df <- rbind(strain_freq_df, new_capsule_row_1, new_capsule_row_2)
      
    }
  }
}

# coerce to numeric

strain_freq_df <- strain_freq_df %>%
  mutate(freq = as.numeric(freq), 
         lower_ci = as.numeric(lower_ci),
         upper_ci = as.numeric(upper_ci)) 

# rename timepoints

strain_freq_df <- strain_freq_df %>%
  rowwise() %>%
  mutate(timepoint_label = paste0("T",str_split(timepoint, " ")[[1]][2]))

day_1 <- as.Date(strain_freq_df %>% filter(timepoint == "timepoint 1") %>% head(1) %>% select(date) %>% pull())


strain_freq_df <- strain_freq_df %>%
  rowwise() %>%
  mutate(day = ifelse(timepoint == "timepoint 1", "Day 1", paste0("Day ", as.character((as.Date(date)-day_1)[[1]] + 1))),
         time_description = ifelse(as.numeric(str_split(time, ":")[[1]][1]) >= 12, "evening", "morning"),
         timepoint_description = paste0(day, ", ", time_description))

# splitting timepoint description duplicate names


strain_freq_df <- strain_freq_df %>%
  group_by(timepoint_description) %>%
  mutate(no_of_repeated_time_descriptions = n_distinct(time),
         repeat_index = dense_rank(time)) %>%
  ungroup() %>%
  mutate(
    timepoint_description = if_else(
      no_of_repeated_time_descriptions == 1,
      as.character(timepoint_description),
      paste0(timepoint_description, " (", repeat_index, ")")
    )
  ) %>%
  select(-no_of_repeated_time_descriptions, -repeat_index)  # optional: drop helper column

# Factor

strain_order <- strain_freq_df$new_strain_label %>% unique() %>% sort() %>% as.character()
strain_freq_df$new_strain_label <- factor(strain_freq_df$new_strain_label, levels = strain_order)
timepoint_order <- strain_freq_df$timepoint_label %>% unique() %>% sort()
timepoint_order <- timepoint_order[order(as.numeric(substring(timepoint_order, 2)))]
strain_freq_df$timepoint_label <- factor(strain_freq_df$timepoint_label, levels = timepoint_order)
timepoint_description_order <- strain_freq_df %>% arrange(day, desc(time)) %>% select(timepoint_description) %>% unique() %>% pull()
strain_freq_df$time_description <- factor(strain_freq_df$time_description, levels = timepoint_description_order)

# Error bar df

error_bars <- strain_freq_df %>% 
  filter(strain == strain_order[2]) %>%
  rowwise() %>%
  mutate(upper_ci = ifelse(strain == new_strain_label, 1 - upper_ci, upper_ci),
         lower_ci = ifelse(strain == new_strain_label, 1 - lower_ci, lower_ci))

# title
plot_title <- paste0(
  "<i>", strsplit(species_oi, "_")[[1]][1], " ",
  strsplit(species_oi, "_")[[1]][2], "</i>, Subject ", subject_oi
)
# Plot
p <- ggplot(strain_freq_df, aes(x = sample_type_number, y = freq, fill = new_strain_label)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.9) +
  geom_errorbar(data = error_bars,
                aes(x = sample_type_number, ymax = upper_ci, ymin = lower_ci), 
                width = 0.2, color = "white") +
  scale_fill_manual(name = "Strain",
                    values = strain_pal) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  theme_bw() + 
  theme(text = element_text(family = "Helvetica"),
        plot.title = element_markdown(size = font_size),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = x_axis_font_size),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = font_size),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(color = "black", size = font_size),
        panel.spacing = unit(0.15, "lines"),
        plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt"),
        legend.text = element_text(size = font_size)
  ) +
  facet_grid(~ timepoint_description, scales = "free", space = "free") +
  labs(title = plot_title, y = "Strain frequency")

p






library(dplyr)
library(ggplot2)
library(cowplot)
library(ggtext)



# functions

source("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/figures/helper_scripts/strain_phasing_helper_scripts.R")


############################################################################
# LOOPING THROUGH SPECIES-SUBJECTS                                         #
############################################################################

# strain_species_subject_df_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/metadata/strain_species_subject_df.txt"
# strain_species_subject_df <- read.csv2(strain_species_subject_df_path, sep = "\t")
# 
# for (i in 1:nrow(strain_species_subject_df)) {
#   species_ <- strain_species_subject_df[i,"species"]
#   subject_ <- strain_species_subject_df[i,"subject"]
#   strain_plot <- plot_strains(species = species_, subject = subject_)
#   out_path <- paste0("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/figures/strain_phasing/strain_phasing_barplots/", species_, "_subject_", as.character(subject_), "_barplot.png")
#   ggsave(out_path, strain_plot, bg = "white", dpi = 300, width = 6, height = 3, units = "in")
# }


############################################################################
# FIGURE                                                                   #
############################################################################
R_obeum_s11 <- plot_strains(species = "Ruminococcus_obeum_61472", subject = 11, timepoints_to_use = c(3,5,6,8), pal = c("#204552", "#7f7e7e"))
out_path = "/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/figures/strain_phasing/Ruminococcus_obeum_61472_QCBio.png"
ggsave(out_path, R_obeum_s11, bg = "white", dpi = 300, width = 7.5, height = 2.5, units = "in")

## Strains
### Temporal
b_vulgatus_s2 <- plot_strains(species = "Bacteroides_vulgatus_57955", subject = 2, timepoints_to_use = c(3,5,6,7))
R_obeum_s11 <- plot_strains(species = "Ruminococcus_obeum_61472", subject = 11, timepoints_to_use = c(3,5,6,8))
### spatial
P_distasonis_s9 <- plot_strains(species = "Parabacteroides_distasonis_56985", subject = 9, timepoints_to_use = c(1,2,3,5)) 
b_vulgatus_s11 <- plot_strains(species = "Bacteroides_vulgatus_57955", subject = 11, timepoints_to_use = c(3,6,9,10))

View(test %>% filter(species == "Ruminococcus_obeum_61472", subject == 11))


## ANOVA
### Load dataframe
strain_freq_df_path <- "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/strain_phasing/strain_freq_GoodSpeciesSubjectPairs.csv"
strain_freq_df <- read.csv2(strain_freq_df_path, sep = ",")

number_of_species = length(strain_freq_df$species %>% unique())
number_of_subjects = length(strain_freq_df$subject %>% unique())
print(paste0("Number of species:", number_of_species))
print(paste0("Number of subjects:", number_of_subjects))

### Initializing table
anova_results_full <- data.frame(species = NA, variable = NA, VarExpl = NA, ss = NA, df = NA, pval = NA)


############################################################################
# SCRATCH                                                                  #
############################################################################
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(lubridate)

# functions

species = "Ruminococcus_obeum_61472"
subject = 11
timepoints_to_use = c(3,5,6,8)
pal =  c("#654321", "#D2B48C")
use_simple_timepoint = FALSE
path = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/strain_phasing/strain_freqs.csv"



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

# Shift timezones

strain_freq_df <- strain_freq_df %>%
  mutate(
    # combine date and time into a single datetime string
    datetime_utc = ymd_hms(paste(date, time), tz = "UTC"),
    
    # convert to Pacific time
    datetime_pst = with_tz(datetime_utc, tzone = "America/Los_Angeles"),
    
    # split
    date_pst = str_split(datetime_pst, " ")[[1]][1],
    time_pst = str_split(datetime_pst, " ")[[1]][2]
  )




# Label with relative days


day_1 <- as.Date(strain_freq_df %>% filter(timepoint == "timepoint 1") %>% head(1) %>% select(date_pst) %>% pull())


strain_freq_df <- strain_freq_df %>%
  rowwise() %>%
  mutate(day = ifelse(timepoint == "timepoint 1", "Day 1", paste0("Day ", as.character((as.Date(date_pst)-day_1)[[1]] + 1))),
         time_description = ifelse(as.numeric(str_split(time_pst, ":")[[1]][1]) >= 20, "evening", "afternoon"),
         time_description = ifelse(as.numeric(str_split(time_pst, ":")[[1]][1]) <= 12, "morning", time_description),
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
timepoint_description_order <- strain_freq_df %>% arrange(date_pst, time_pst) %>% select(timepoint_description) %>% unique() %>% pull()
strain_freq_df$timepoint_description <- factor(strain_freq_df$timepoint_description, levels = timepoint_description_order)


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



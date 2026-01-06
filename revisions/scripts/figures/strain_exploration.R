 
##################################################################################################
# Delta f method
##################################################################################################
 
# Clearing environment
rm(list=ls()) 

# packages

library(dplyr)
library(compositions) 



# Humanized mice ######

# parameters

text_size = 12
subtext_size = 10

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

project_folder <- "/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/"

all_strain_abundance <- lapply(species_list, function(species) {
  
  strain_abundance_path <- paste0(
    project_folder,
    "strain_phasing/strain_clusters/", species, "/", species, "_strain_frequency.csv"
  )
  
  read.csv2(strain_abundance_path, sep = "\t") %>%
    mutate(
      species    = species,
      freq       = as.numeric(freq),
      freq_trans = as.vector(clr(freq))
    )
}) %>% bind_rows()


# Step 1. Filter for strain 1 in mouse samples only (no inoculum)

all_strain_abundance <- all_strain_abundance %>%
  filter(strain == 1, mouse_number != "Inoculum")

all_strain_abundance <- all_strain_abundance %>%
  mutate(
    taxonomy = sapply(species, function(species_oi) {
      parts <- str_split(species_oi, "_", simplify = TRUE)
      if (parts[2] == "sp") {
        paste0(parts[1], " sp.")
      } else {
        paste0(parts[1], " ", parts[2])
      }
    })
  )

# Step 2. Keep species with strain abundances that are (1) non-zero in at least one sample (2) present in at least two gut sites and two mice

species_with_two_mice_two_regions_and_variable_freq <- all_strain_abundance %>% 
  group_by(species, mouse_number) %>% 
  summarise(
    distinct_regions = n_distinct(region),
    .groups = "drop"
  ) %>% 
  filter(distinct_regions >= 2) %>%                 # mice with ≥2 gut sites
  group_by(species) %>% 
  summarise(
    mice_with_2_regions = n_distinct(mouse_number),
    .groups = "drop"
  ) %>% 
  filter(mice_with_2_regions >= 2) %>%               # ≥2 such mice
  inner_join(
    all_strain_abundance %>% 
      group_by(species) %>% 
      summarise(
        has_variable_freq = any(freq > 0 & freq < 1, na.rm = TRUE),
        .groups = "drop"
      ) %>% 
      filter(has_variable_freq),
    by = "species"
  ) %>% 
  pull(species)

filtered_all_strain_abundance <- all_strain_abundance %>% 
  filter(species %in% species_with_two_mice_two_regions_and_variable_freq)

# Step 3. Pick sample pairs
## Within-host: farthest along the gut
## Between host, within same cage: closest to the same site along the gut
## Between host, between cages: closest to the same site along the gut

region_order <- c("Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")

annotated_df <- all_strain_abundance %>% 
  mutate(region_idx = match(region, region_order)) %>% 
  filter(!is.na(region_idx))

### Within-host 

within_host_comparisons <- annotated_df %>% 
  group_by(species,taxonomy,cage, mouse_number) %>% 
  filter(n_distinct(region_idx) >= 2) %>% 
  summarise(
    sample_1 = sample[which.min(region_idx)],
    sample_2 = sample[which.max(region_idx)],
    region_1 = region[which.min(region_idx)],
    region_2 = region[which.max(region_idx)],
    freq_1 = freq[which.min(region_idx)],
    freq_2 = freq[which.max(region_idx)],
    freq_trans_1 = freq_trans[which.min(region_idx)],
    freq_trans_2 = freq_trans[which.max(region_idx)],
    comparison_type = "within_host",
    .groups = "drop"
  ) %>%
  transmute(
    species,
    taxonomy,
    cage_comparison = cage,
    mouse_number_1 = mouse_number,
    mouse_number_2 = mouse_number,
    cage_1 = cage,
    cage_2 = cage,
    sample_1,
    sample_2,
    region_1,
    region_2,
    comparison_type,
    freq_1,
    freq_2,
    freq_trans_1,
    freq_trans_2
  )



## Between mouse, same cage

same_cage_mouse_pairs <- annotated_df %>% 
  distinct(species,taxonomy, cage, mouse_number) %>% 
  inner_join(
    distinct(annotated_df, species, cage, mouse_number),
    by = c("species", "cage"),
    suffix = c("_1", "_2")
  ) %>% 
  filter(mouse_number_1 < mouse_number_2)


same_cage_candidates <- same_cage_mouse_pairs %>% 
  inner_join(
    annotated_df,
    by = c("species","mouse_number_1" = "mouse_number")
  ) %>% 
  rename(
    sample_1 = sample,
    region_1 = region,
    region_idx_1 = region_idx,
    freq_1 = freq,
    freq_trans_1 = freq_trans
  ) %>% 
  inner_join(
    annotated_df,
    by = c("species", "mouse_number_2" = "mouse_number")
  ) %>% 
  rename(
    sample_2 = sample,
    region_2 = region,
    region_idx_2 = region_idx,
    freq_2 = freq,
    freq_trans_2 = freq_trans
  ) %>% 
  mutate(
    region_distance = abs(region_idx_1 - region_idx_2),
    distal_bias = pmax(region_idx_1, region_idx_2)
  )

between_host_same_cage <- same_cage_candidates %>% 
  group_by(species,taxonomy, cage) %>% 
  arrange(region_distance, desc(distal_bias)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  transmute(
    species,
    taxonomy,
    cage_comparison = cage,
    mouse_number_1,
    mouse_number_2,
    cage_1 = cage,
    cage_2 = cage,
    sample_1,
    sample_2,
    region_1,
    region_2,
    comparison_type = "between_host_same_cage",
    freq_1,
    freq_2,
    freq_trans_1,
    freq_trans_2
  )


# Between host, between cage

cage_pairs <- annotated_df %>% 
  distinct(species,taxonomy,cage) %>% 
  inner_join(
    distinct(annotated_df, species, cage),
    by = "species",
    suffix = c("_1", "_2")
  ) %>% 
  filter(cage_1 < cage_2)

between_cage_mouse_pairs <- cage_pairs %>% 
  inner_join(
    distinct(annotated_df, species, cage, mouse_number),
    by = c("species", "cage_1" = "cage")
  ) %>% 
  rename(mouse_number_1 = mouse_number) %>% 
  inner_join(
    distinct(annotated_df, species, cage, mouse_number),
    by = c("species", "cage_2" = "cage")
  ) %>% 
  rename(mouse_number_2 = mouse_number)

between_cage_candidates <- between_cage_mouse_pairs %>% 
  inner_join(
    annotated_df,
    by = c("species","taxonomy", "mouse_number_1" = "mouse_number")
  ) %>% 
  rename(
    sample_1 = sample,
    region_1 = region,
    region_idx_1 = region_idx,
    freq_1 = freq,
    freq_trans_1 = freq_trans
  ) %>% 
  inner_join(
    annotated_df,
    by = c("species", "taxonomy", "mouse_number_2" = "mouse_number")
  ) %>% 
  rename(
    sample_2 = sample,
    region_2 = region,
    region_idx_2 = region_idx,
    freq_2 = freq,
    freq_trans_2 = freq_trans
  ) %>% 
  mutate(
    region_distance = abs(region_idx_1 - region_idx_2),
    distal_bias = pmax(region_idx_1, region_idx_2)
  )

between_host_between_cage <- between_cage_candidates %>%
  group_by(species, taxonomy) %>%
  group_modify(~ {
    used_mice <- character()
    out <- list()
    
    pick_pair <- function(df) {
      df <- df %>%
        arrange(region_distance, desc(distal_bias)) %>%
        filter(
          !(mouse_number_1 %in% used_mice),
          !(mouse_number_2 %in% used_mice)
        )
      
      if (nrow(df) == 0) return(NULL)
      
      row <- df[1, ]
      used_mice <<- c(used_mice, row$mouse_number_1, row$mouse_number_2)
      row
    }
    
    # Cage 1 vs Cage 2
    g12 <- .x %>% filter(cage_1 == "Cage 1", cage_2 == "Cage 2")
    r12 <- pick_pair(g12)
    if (!is.null(r12)) out[[length(out) + 1]] <- r12
    
    # Cage 1 vs Cage 3
    g13 <- .x %>% filter(cage_1 == "Cage 1", cage_2 == "Cage 3")
    r13 <- pick_pair(g13)
    if (!is.null(r13)) out[[length(out) + 1]] <- r13
    
    # Cage 2 vs Cage 3
    g23 <- .x %>% filter(cage_1 == "Cage 2", cage_2 == "Cage 3")
    r23 <- pick_pair(g23)
    if (!is.null(r23)) out[[length(out) + 1]] <- r23
    
    bind_rows(out)
  }) %>%
  ungroup() %>% 
  transmute(
    species,
    taxonomy,
    cage_comparison = paste(cage_1, cage_2, sep = "_vs_"),
    mouse_number_1,
    mouse_number_2,
    cage_1,
    cage_2,
    sample_1,
    sample_2,
    region_1,
    region_2,
    comparison_type = "between_host_between_cage",
    freq_1,
    freq_2,
    freq_trans_1,
    freq_trans_2
  )



## Combining dataframes

final_comparisons_df <- bind_rows(
  within_host_comparisons,
  between_host_same_cage,
  between_host_between_cage
)

# Step 5. Add delta f

final_comparisons_df <- final_comparisons_df %>%
  rowwise() %>%
  mutate(delta_f = abs(freq_2 - freq_1))
  

# Step 6: keep only species with at least one observation in all categories

required_types <- c(
  "within_host",
  "between_host_same_cage",
  "between_host_between_cage"
)

species_with_all_categories <- final_comparisons_df %>% 
  group_by(species) %>% 
  summarise(
    n_types = n_distinct(comparison_type),
    .groups = "drop"
  ) %>% 
  filter(n_types == length(required_types)) %>% 
  pull(species)

final_comparisons_df <- final_comparisons_df %>% 
  filter(species %in% species_with_all_categories)

final_comparisons_df <- final_comparisons_df %>%
  rowwise() %>%
  mutate(comparison_type_label = ifelse(comparison_type == "within_host", "Within host", comparison_type),
         comparison_type_label = ifelse(comparison_type == "between_host_same_cage", "Between host,\nwithin cage", comparison_type_label),
         comparison_type_label = ifelse(comparison_type == "between_host_between_cage", "Between host,\nbetween cage", comparison_type_label))

final_comparisons_df$comparison_type_label <- factor(final_comparisons_df$comparison_type_label, levels = c("Between host,\nbetween cage","Between host,\nwithin cage", "Within host"))


# ASIDE: CALCULATING MEDIAN CHANGE

median_change <- final_comparisons_df %>%
  group_by(comparison_type) %>%
  summarise(median_change = median(delta_f))

median_change

# Step 8. Plotting

p <- ggplot(
  final_comparisons_df,
  aes(
    x = taxonomy,
    y = delta_f,
    fill = comparison_type_label
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.8), 
    alpha = 0.5
  ) +
  geom_jitter(
    aes(color = comparison_type_label),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 2,
    # color = "black"
  ) +
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
  ) +
  labs(
    x = NULL,
    y = expression(Delta~f),
    fill = "Comparison type",
    color = "Comparison type"
  ) +
  scale_fill_manual(values = c("Between host,\nbetween cage" = "#AA98A9", 
                               "Between host,\nwithin cage" = "#BDB5D5", 
                               "Within host" = "#CC0000")) + 
  scale_color_manual(values = c("Between host,\nbetween cage" = "#AA98A9", 
                                "Between host,\nwithin cage" = "#BDB5D5", 
                                "Within host" = "#CC0000")) +
  ylim(0,1)
# p

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/HM_strain_delta_f.png"
ggsave(out_path, p, dpi = 300, height = 4, width = 7)

# Conventional mice ######

# Load species list
strain_phasing_dir <- "/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/strain_phasing/strain_clusters/"
species_list <- list.dirs(strain_phasing_dir, full.names = FALSE)[-1]

# making dataframe
strain_abundance_list <- list()
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

  strain_abundance_list[[species]] <- strain_abundance
}

all_strain_abundance <- bind_rows(strain_abundance_list)

species_metadata_path <- "/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc/metadata.tsv"
species_metadata_df <- read.csv2(species_metadata_path, sep = "\t") 
species_metadata_df$species_id <- as.character(species_metadata_df$species_id) 

all_strain_abundance <- all_strain_abundance %>%
  rowwise() %>% 
  mutate(taxonomy = str_split(species_metadata_df[species_metadata_df['species_id'] == species,"species"], "__")[[1]][2],
         taxonomy = ifelse((length(str_split(taxonomy, " ")[[1]]) == 1) | (substr(str_split(taxonomy, " ")[[1]][2], 1, 2) == "sp"), paste0(str_split(taxonomy, " ")[[1]][1], " sp."), taxonomy))



# Step 1. Filter for strain 1 in mouse samples only (no inoculum) and rename mice

all_strain_abundance <- all_strain_abundance %>%
  filter(strain == 1, mouse_number != "Inoculum")

transform_mouse_number <- function(mouse_str,cage_str) {
  mouse_number <- as.numeric(str_split(mouse_str, " ")[[1]][2])
  cage_number <- as.numeric(str_split(cage_str, " ")[[1]][2])
  real_mouse_str <- paste0("Mouse ", as.character(mouse_number + 3*(cage_number-1)))
  return(real_mouse_str)
}

all_strain_abundance <- all_strain_abundance %>%
  rowwise() %>%
  mutate(mouse_number = transform_mouse_number(mouse_number, cage))

# Step 2. Keep species with strain abundances that are (1) non-zero in at least one sample (2) present in at least two gut sites and two mice

species_with_two_mice_two_regions_and_variable_freq <- all_strain_abundance %>% 
  group_by(species, mouse_number) %>% 
  summarise(
    distinct_regions = n_distinct(region),
    .groups = "drop"
  ) %>% 
  filter(distinct_regions >= 2) %>%                 # mice with ≥2 gut sites
  group_by(species) %>% 
  summarise(
    mice_with_2_regions = n_distinct(mouse_number),
    .groups = "drop"
  ) %>% 
  filter(mice_with_2_regions >= 2) %>%               # ≥2 such mice
  inner_join(
    all_strain_abundance %>% 
      group_by(species) %>% 
      summarise(
        has_variable_freq = any(freq > 0 & freq < 1, na.rm = TRUE),
        .groups = "drop"
      ) %>% 
      filter(has_variable_freq),
    by = "species"
  ) %>% 
  pull(species)


all_strain_abundance <- all_strain_abundance %>% 
  filter(species %in% species_with_two_mice_two_regions_and_variable_freq)


# Step 3. Pick sample pairs
## Within-host: farthest along the gut
## Between host, within same cage: closest to the same site along the gut
## Between host, between cages: closest to the same site along the gut

region_order <- c("Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")

annotated_df <- all_strain_abundance %>% 
  mutate(region_idx = match(region, region_order)) %>% 
  filter(!is.na(region_idx))

### Within-host 

within_host_comparisons <- annotated_df %>% 
  group_by(species,taxonomy,cage, mouse_number) %>% 
  filter(n_distinct(region_idx) >= 2) %>% 
  summarise(
    sample_1 = sample[which.min(region_idx)],
    sample_2 = sample[which.max(region_idx)],
    region_1 = region[which.min(region_idx)],
    region_2 = region[which.max(region_idx)],
    freq_1 = freq[which.min(region_idx)],
    freq_2 = freq[which.max(region_idx)],
    freq_trans_1 = freq_trans[which.min(region_idx)],
    freq_trans_2 = freq_trans[which.max(region_idx)],
    comparison_type = "within_host",
    .groups = "drop"
  ) %>%
  transmute(
    species,
    taxonomy,
    cage_comparison = cage,
    mouse_number_1 = mouse_number,
    mouse_number_2 = mouse_number,
    cage_1 = cage,
    cage_2 = cage,
    sample_1,
    sample_2,
    region_1,
    region_2,
    comparison_type,
    freq_1,
    freq_2,
    freq_trans_1,
    freq_trans_2
  )

  

## Between mouse, same cage

same_cage_mouse_pairs <- annotated_df %>% 
  distinct(species,taxonomy, cage, mouse_number) %>% 
  inner_join(
    distinct(annotated_df, species, cage, mouse_number),
    by = c("species", "cage"),
    suffix = c("_1", "_2")
  ) %>% 
  filter(mouse_number_1 < mouse_number_2)


same_cage_candidates <- same_cage_mouse_pairs %>% 
  inner_join(
    annotated_df,
    by = c("species","mouse_number_1" = "mouse_number")
  ) %>% 
  rename(
    sample_1 = sample,
    region_1 = region,
    region_idx_1 = region_idx,
    freq_1 = freq,
    freq_trans_1 = freq_trans
  ) %>% 
  inner_join(
    annotated_df,
    by = c("species", "mouse_number_2" = "mouse_number")
  ) %>% 
  rename(
    sample_2 = sample,
    region_2 = region,
    region_idx_2 = region_idx,
    freq_2 = freq,
    freq_trans_2 = freq_trans
  ) %>% 
  mutate(
    region_distance = abs(region_idx_1 - region_idx_2),
    distal_bias = pmax(region_idx_1, region_idx_2)
  )

between_host_same_cage <- same_cage_candidates %>% 
  group_by(species,taxonomy, cage) %>% 
  arrange(region_distance, desc(distal_bias)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  transmute(
    species,
    taxonomy,
    cage_comparison = cage,
    mouse_number_1,
    mouse_number_2,
    cage_1 = cage,
    cage_2 = cage,
    sample_1,
    sample_2,
    region_1,
    region_2,
    comparison_type = "between_host_same_cage",
    freq_1,
    freq_2,
    freq_trans_1,
    freq_trans_2
  )


# Between host, between cage

cage_pairs <- annotated_df %>% 
  distinct(species,taxonomy,cage) %>% 
  inner_join(
    distinct(annotated_df, species, cage),
    by = "species",
    suffix = c("_1", "_2")
  ) %>% 
  filter(cage_1 < cage_2)

between_cage_mouse_pairs <- cage_pairs %>% 
  inner_join(
    distinct(annotated_df, species, cage, mouse_number),
    by = c("species", "cage_1" = "cage")
  ) %>% 
  rename(mouse_number_1 = mouse_number) %>% 
  inner_join(
    distinct(annotated_df, species, cage, mouse_number),
    by = c("species", "cage_2" = "cage")
  ) %>% 
  rename(mouse_number_2 = mouse_number)

between_cage_candidates <- between_cage_mouse_pairs %>% 
  inner_join(
    annotated_df,
    by = c("species","taxonomy", "mouse_number_1" = "mouse_number")
  ) %>% 
  rename(
    sample_1 = sample,
    region_1 = region,
    region_idx_1 = region_idx,
    freq_1 = freq,
    freq_trans_1 = freq_trans
  ) %>% 
  inner_join(
    annotated_df,
    by = c("species", "taxonomy", "mouse_number_2" = "mouse_number")
  ) %>% 
  rename(
    sample_2 = sample,
    region_2 = region,
    region_idx_2 = region_idx,
    freq_2 = freq,
    freq_trans_2 = freq_trans
  ) %>% 
  mutate(
    region_distance = abs(region_idx_1 - region_idx_2),
    distal_bias = pmax(region_idx_1, region_idx_2)
  )

between_host_between_cage <- between_cage_candidates %>%
  group_by(species, taxonomy) %>%
  group_modify(~ {
    used_mice <- character()
    out <- list()
    
    pick_pair <- function(df) {
      df <- df %>%
        arrange(region_distance, desc(distal_bias)) %>%
        filter(
          !(mouse_number_1 %in% used_mice),
          !(mouse_number_2 %in% used_mice)
        )
      
      if (nrow(df) == 0) return(NULL)
      
      row <- df[1, ]
      used_mice <<- c(used_mice, row$mouse_number_1, row$mouse_number_2)
      row
    }
    
    # Cage 1 vs Cage 2
    g12 <- .x %>% filter(cage_1 == "Cage 1", cage_2 == "Cage 2")
    r12 <- pick_pair(g12)
    if (!is.null(r12)) out[[length(out) + 1]] <- r12
    
    # Cage 1 vs Cage 3
    g13 <- .x %>% filter(cage_1 == "Cage 1", cage_2 == "Cage 3")
    r13 <- pick_pair(g13)
    if (!is.null(r13)) out[[length(out) + 1]] <- r13
    
    # Cage 2 vs Cage 3
    g23 <- .x %>% filter(cage_1 == "Cage 2", cage_2 == "Cage 3")
    r23 <- pick_pair(g23)
    if (!is.null(r23)) out[[length(out) + 1]] <- r23
    
    bind_rows(out)
  }) %>%
  ungroup() %>% 
  transmute(
    species,
    taxonomy,
    cage_comparison = paste(cage_1, cage_2, sep = "_vs_"),
    mouse_number_1,
    mouse_number_2,
    cage_1,
    cage_2,
    sample_1,
    sample_2,
    region_1,
    region_2,
    comparison_type = "between_host_between_cage",
    freq_1,
    freq_2,
    freq_trans_1,
    freq_trans_2
  )



## Combining dataframes

final_comparisons_df <- bind_rows(
  within_host_comparisons,
  between_host_same_cage,
  between_host_between_cage
)

# Step 5. Add delta f

final_comparisons_df <- final_comparisons_df %>%
  rowwise() %>%
  mutate(delta_f = abs(freq_2 - freq_1))

# Step 6: keep only species with at least one observation in all categories

required_types <- c(
  "within_host",
  "between_host_same_cage",
  "between_host_between_cage"
)

species_with_all_categories <- final_comparisons_df %>% 
  group_by(species) %>% 
  summarise(
    n_types = n_distinct(comparison_type),
    .groups = "drop"
  ) %>% 
  filter(n_types == length(required_types)) %>% 
  pull(species)

final_comparisons_df <- final_comparisons_df %>% 
  filter(species %in% species_with_all_categories)

final_comparisons_df <- final_comparisons_df %>%
  rowwise() %>%
  mutate(comparison_type_label = ifelse(comparison_type == "within_host", "Within host", comparison_type),
         comparison_type_label = ifelse(comparison_type == "between_host_same_cage", "Between host,\nwithin cage", comparison_type_label),
         comparison_type_label = ifelse(comparison_type == "between_host_between_cage", "Between host,\nbetween cage", comparison_type_label))

final_comparisons_df$comparison_type_label <- factor(final_comparisons_df$comparison_type_label, levels = c("Between host,\nbetween cage","Between host,\nwithin cage", "Within host"))


# Step 8. Plotting

p <- ggplot(
  final_comparisons_df,
  aes(
    x = taxonomy,
    y = delta_f,
    fill = comparison_type_label
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.8), 
    alpha = 0.5
  ) +
  geom_jitter(
    aes(color = comparison_type_label),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 2
  ) +
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
  ) +
  labs(
    x = NULL,
    y = expression(Delta~f),
    fill = "Comparison type",
    color = "Comparison type"
  ) +
  scale_fill_manual(values = c("Between host,\nbetween cage" = "#AA98A9", 
                               "Between host,\nwithin cage" = "#BDB5D5", 
                               "Within host" = "#CC0000")) + 
  scale_color_manual(values = c("Between host,\nbetween cage" = "#AA98A9", 
                                "Between host,\nwithin cage" = "#BDB5D5", 
                                "Within host" = "#CC0000")) +
  ylim(0,1)
# p

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/CM_strain_delta_f.png"
ggsave(out_path, p, dpi = 300, height = 4, width = 7)

# Shalon ######

text_size = 12
subtext_size = 10

### Load dataframe
strain_freq_df_path <- "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/strain_phasing/strain_freqs_filtered.csv"
strain_freq_df <- read.csv2(strain_freq_df_path, sep = ",") %>%
  mutate(freq = as.numeric(freq),
         freq_trans = as.vector(clr(freq)))

### Choosing species using the anova prevalence filter
anova_prevalence_filter <- TRUE

if (anova_prevalence_filter) {
  species_subject_list <- data.frame(
    species = c(
      "Acidaminococcus_intestini_54097",
      "Adlercreutzia_equolifaciens_60310",
      "Anaerostipes_hadrus_55206",
      "Anaerostipes_hadrus_55206",
      "Bacteroides_vulgatus_57955",
      "Bacteroides_vulgatus_57955",
      "Bacteroides_vulgatus_57955",
      "Bifidobacterium_adolescentis_56815",
      "Bifidobacterium_longum_57796",
      "Eubacterium_rectale_56927",
      "Guyana_massiliensis_60772",
      "Parabacteroides_distasonis_56985",
      "Ruminococcus_obeum_61472",
      "Ruminococcus_sp_55468"
    ),
    subject = c(
      9, 12, 6, 8, 2, 8, 11, 1, 5, 9, 12, 9, 11, 6
    ),
    stringsAsFactors = FALSE
  )

} else {
  
  ### species list
  
  species_list <- strain_freq_df %>% select(species) %>% unique() %>% pull()
  
  species_subject_list <- strain_freq_df %>% select(species, subject) %>% unique()
  
  
}


# Step 1. Filter for strain 1 and certain prevalence requirements



eligible_species_subjects <- strain_freq_df %>%
  filter(strain == 1) %>%
  
  # find timepoints with >1 sample
  group_by(species, subject, timepoint) %>%
  summarise(
    n_samples = n_distinct(sample_type_number),
    .groups = "drop"
  ) %>%
  filter(n_samples > 1) %>%
  
  # require >1 such timepoint
  group_by(species, subject) %>%
  filter(n_distinct(timepoint) > 1) %>%
  ungroup() %>%
  select(species, subject)

# keep ALL rows for eligible (species, subject)
if (anova_prevalence_filter) {
  all_strain_abundance <- strain_freq_df %>%
    filter(strain == 1) %>%
    semi_join(species_subject_list, 
              by = c("species", "subject"))
} else {
  all_strain_abundance <- strain_freq_df %>%
    filter(strain == 1) %>%
    semi_join(
      eligible_species_subjects,
      by = c("species", "subject")
    )
}

# Step 2. Drop technical replicates in subject 1 and add taxonomy

all_strain_abundance <- all_strain_abundance %>%
  group_by(species, subject, sample_type_number, timepoint) %>%
  slice(1) %>%
  ungroup()

all_strain_abundance <- all_strain_abundance %>%
  mutate(
    taxonomy = sapply(species, function(species_oi) {
      parts <- str_split(species_oi, "_", simplify = TRUE)
      if (parts[2] == "sp") {
        paste0(parts[1], " sp.")
      } else {
        paste0(parts[1], " ", parts[2])
      }
    })
  )


# Step 3. Pick sample pairs
## Within-timepoint: farthest along the gut
## Between timepoint (consecutive): same capsule if possible
## Between timepoint (farthest measurement: same capsule if possible

region_order <- c("Capsule 1", "Capsule 2", "Capsule 3", "Capsule 4")

annotated_df <- all_strain_abundance %>% 
  mutate(region_idx = match(sample_type_number, region_order)) %>% 
  filter(!is.na(region_idx))

# Within timepoint 

within_timepoint_comparisons <- annotated_df %>% 
  group_by(species,taxonomy,subject, date, time, timepoint) %>% 
  filter(n_distinct(region_idx) >= 2) %>% 
  summarise(
    sample_type_number_1 = sample_type_number[which.min(region_idx)],
    sample_type_number_2 = sample_type_number[which.max(region_idx)],
    freq_1 = freq[which.min(region_idx)],
    freq_2 = freq[which.max(region_idx)],
    freq_trans_1 = freq_trans[which.min(region_idx)],
    freq_trans_2 = freq_trans[which.max(region_idx)],
    comparison_type = "within_timepoint",
    .groups = "drop"
  ) %>%
  transmute(
    species,
    taxonomy,
    subject,
    tmpt_comparison = timepoint,
    date_1 = date,
    date_2 = date,
    time_1 = time,
    time_2 = time,
    sample_type_number_1,
    sample_type_number_2,
    comparison_type,
    freq_1,
    freq_2,
    freq_trans_1,
    freq_trans_2
  )

# Between timepoint (consecutive)

candidate_pairs <- annotated_df %>%
  select(
    species, taxonomy, subject, 
    date, time, timepoint,
    sample_type_number,
    region_idx_1 = region_idx,
    freq_1 = freq,
    freq_trans_1 = freq_trans
  ) %>%
  inner_join(
    annotated_df %>%
      select(
        species, taxonomy, subject,
        date, time, timepoint,
        sample_type_number,
        region_idx_2 = region_idx,
        freq_2 = freq,
        freq_trans_2 = freq_trans
      ),
    by = c("species", "subject", "taxonomy"),
    suffix = c("_1", "_2")
  ) %>%
  filter(timepoint_2 == timepoint_1 + 1) %>%
  mutate(
    region_distance = abs(region_idx_1 - region_idx_2),
    comparison_type = "between_timepoint_consecutive"
  )

consecutive_timepoint_pairs <- candidate_pairs %>%
  group_by(species, taxonomy, subject,date_1,date_2, time_1,time_2, timepoint_1, timepoint_2) %>%
  arrange(region_distance) %>%
  slice(1) %>%
  ungroup()

# Longest timepoint

timepoint_bounds <- annotated_df %>%
  group_by(species, subject) %>%
  summarise(
    t_min = min(timepoint),
    t_max = max(timepoint),
    .groups = "drop"
  ) %>%
  filter(t_min < t_max)

farthest_timepoint_candidates <- timepoint_bounds %>%
  inner_join(
    annotated_df,
    by = c("species", "subject")
  ) %>%
  filter(timepoint == t_min | timepoint == t_max) %>%
  select(
    species, subject,
    t_min, t_max,
    date_1=date, time_1=time,timepoint_1 = timepoint,
    sample_type_number_1 = sample_type_number,
    region_idx_1 = region_idx,
    freq_1=freq,
    freq_trans_1=freq_trans
  ) %>%
  inner_join(
    annotated_df %>%
      select(
        species, taxonomy, subject,
        date_2=date, time_2=time, timepoint_2 = timepoint,
        sample_type_number_2 = sample_type_number,
        region_idx_2 = region_idx,
        freq_2 = freq,
        freq_trans_2 = freq_trans
      ),
    by = c("species", "subject")
  ) %>%
  filter(
    timepoint_1 == t_min,
    timepoint_2 == t_max
  ) %>%
  mutate(
    region_distance = abs(region_idx_1 - region_idx_2),
    distal_bias = pmax(region_idx_1, region_idx_2),
    comparison_type = "between_timepoint_long"
  )

farthest_timepoint_comparisons <- farthest_timepoint_candidates %>%
  group_by(species,taxonomy, subject) %>%
  arrange(
    region_distance,
    desc(distal_bias)
  ) %>%
  slice(1) %>%
  ungroup()

# Step 4. Make final plotting df

within_t <- within_timepoint_comparisons %>%
  transmute(species, taxonomy, subject,
            date_1, date_2, 
            time_1, time_2, 
            timepoint_1 = tmpt_comparison,
            timepoint_2 = tmpt_comparison,
            sample_type_number_1,
            sample_type_number_2,
            comparison_type,
            freq_1, freq_2,
            freq_trans_1,freq_trans_2)

between_t_cons <- consecutive_timepoint_pairs %>%
  transmute(species, taxonomy, subject,
            date_1, date_2, 
            time_1, time_2, 
            timepoint_1,
            timepoint_2,
            sample_type_number_1,
            sample_type_number_2,
            comparison_type,
            freq_1, freq_2,
            freq_trans_1,freq_trans_2)

between_t_long <- farthest_timepoint_comparisons %>%
  transmute(species, taxonomy, subject,
            date_1, date_2, 
            time_1, time_2, 
            timepoint_1,
            timepoint_2,
            sample_type_number_1,
            sample_type_number_2,
            comparison_type,
            freq_1, freq_2,
            freq_trans_1,freq_trans_2)

final_comparisons_df <- bind_rows(
  within_t,
  between_t_cons,
  between_t_long
)


# Step 5. Add delta f

final_comparisons_df <- final_comparisons_df %>%
  rowwise() %>%
  mutate(delta_f = abs(freq_2 - freq_1))

# Step 6: keep only species with at least one observation in all categories

required_types <- c(
  "within_timepoint",
  "between_timepoint_consecutive",
  "between_timepoint_long"
)

species_with_all_categories <- final_comparisons_df %>% 
  group_by(species, subject) %>% 
  summarise(
    n_types = n_distinct(comparison_type),
    .groups = "drop"
  ) %>% 
  filter(n_types == length(required_types)) %>% 
  select(species, subject)

final_comparisons_df <- final_comparisons_df %>%
  semi_join(
    species_with_all_categories,
    by = c("species", "subject")
  )

final_comparisons_df <- final_comparisons_df %>%
  rowwise() %>%
  mutate(comparison_type_label = ifelse(comparison_type == "within_timepoint", "Within timepoint,\nbetween gut region", comparison_type),
         comparison_type_label = ifelse(comparison_type == "between_timepoint_consecutive", "Between consecutive\ntimepoints", comparison_type_label),
         comparison_type_label = ifelse(comparison_type == "between_timepoint_long", "Between first and\nlast timepoint", comparison_type_label))

final_comparisons_df$comparison_type_label <- factor(final_comparisons_df$comparison_type_label, levels = c("Between first and\nlast timepoint","Between consecutive\ntimepoints", "Within timepoint,\nbetween gut region"))

# Step 8. Add species, subject label

final_comparisons_df <- final_comparisons_df %>%
  rowwise() %>%
  mutate(x_label = paste0(taxonomy, ",\nSubject ", subject))

x_label_order <- final_comparisons_df %>%
  select(x_label) %>%
  unique() %>%
  arrange() %>% 
  pull()

# Step 9. Plotting

p <- ggplot(
  final_comparisons_df,
  aes(
    x = x_label,
    y = delta_f,
    fill = comparison_type_label
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.8), 
    alpha = 5
  ) +
  geom_jitter(
    aes(color = comparison_type_label),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 2
  ) +
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = text_size),
        axis.text.x = element_text(angle = 45, hjust = 1, size = subtext_size, lineheight = 0.8),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = subtext_size),
        axis.title.y = element_text(size = text_size),
        title = element_text(size = text_size),
        legend.title = element_blank(),
        legend.text = element_text(size=subtext_size),
        plot.margin = unit(c(0.5, 0.5, 0.5, 2.5), "lines"),
  ) +
  labs(
    x = NULL,
    y = expression(Delta~f),
    fill = "Comparison type",
    color = "Comparison type"
  ) +
  scale_fill_manual(values = c("Between first and\nlast timepoint" = "#855168",
                               "Between consecutive\ntimepoints" = "#e0c9d3",
                               "Within timepoint,\nbetween gut region" = "#CC0000")) + 
  scale_color_manual(values = c("Between first and\nlast timepoint" = "#855168",
                               "Between consecutive\ntimepoints" = "#e0c9d3",
                               "Within timepoint,\nbetween gut region" = "#CC0000")) + 
  ylim(0,1)
# p

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/shalon_strain_delta_f_subset.png"
ggsave(out_path, p, dpi = 300, height = 4, width = 8)





##################################################################################################
# Variance method
##################################################################################################

# Clearing environment
rm(list=ls()) 

# packages

library(dplyr)
library(compositions)  

# parameters

text_size = 12
subtext_size = 10

# Humanized mice ######

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



project_folder <- "/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/"

all_strain_abundance <- lapply(species_list, function(species) {
  
  strain_abundance_path <- paste0(
    project_folder,
    "strain_phasing/strain_clusters/", species, "/", species, "_strain_frequency.csv"
  )
  
  read.csv2(strain_abundance_path, sep = "\t") %>%
    mutate(
      species    = species,
      freq       = as.numeric(freq),
      freq_trans = as.vector(clr(freq))
    )
}) %>% bind_rows()


# Step 1. Filter for strain 1 in mouse samples only (no inoculum)

all_strain_abundance <- all_strain_abundance %>%
  filter(strain == 1, mouse_number != "Inoculum")

# Step 2. Keep species with strain abundances that are (1) non-zero in at least one sample (2) present in at least two gut sites and two mice


species_with_two_mice_two_regions_and_variable_freq <- all_strain_abundance %>% 
  group_by(species, mouse_number) %>% 
  summarise(
    distinct_regions = n_distinct(region),
    .groups = "drop"
  ) %>% 
  filter(distinct_regions >= 2) %>%                 # mice with ≥2 gut sites
  group_by(species) %>% 
  summarise(
    mice_with_2_regions = n_distinct(mouse_number),
    .groups = "drop"
  ) %>% 
  filter(mice_with_2_regions >= 2) %>%               # ≥2 such mice
  inner_join(
    all_strain_abundance %>% 
      group_by(species) %>% 
      summarise(
        has_variable_freq = any(freq > 0 & freq < 1, na.rm = TRUE),
        .groups = "drop"
      ) %>% 
      filter(has_variable_freq),
    by = "species"
  ) %>% 
  pull(species)

filtered_all_strain_abundance <- all_strain_abundance %>% 
  filter(species %in% species_with_two_mice_two_regions_and_variable_freq)

# Step 3. Calculate variance

## Within mouse

within_mouse_stats <- filtered_all_strain_abundance %>%
  group_by(species, cage, mouse_number) %>%
  summarise(
    mean_within_mouse = mean(freq),
    var_within_mouse = var(freq),
    .groups = "drop"
  )

## Between mouse

between_mouse_stats <- within_mouse_stats %>%
  group_by(species, cage) %>%
  summarise(
    var_between_mice = var(mean_within_mouse),
    mean_within_cage = mean(mean_within_mouse),
    .groups = "drop"
  )

## Between cage

between_cage_stats <- between_mouse_stats %>%
  group_by(species) %>%
  summarise(
    var_between_cages = var(mean_within_cage),
    .groups = "drop"
  )

## Combine

species_variances <- within_mouse_stats %>%
  group_by(species) %>%
  summarise(
    within_mouse_var = mean(var_within_mouse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    between_mouse_stats %>%
      group_by(species) %>%
      summarise(between_mouse_var = mean(var_between_mice, na.rm = TRUE), .groups = "drop"),
    by = "species"
  ) %>%
  left_join(
    between_cage_stats,
    by = "species"
  )

# Step 4. remove species without all 4 values

species_variances <- species_variances %>%
  drop_na() 

# Step 5. Convert to longform 

species_variances_long <- species_variances %>%
  pivot_longer(
    cols = c(within_mouse_var, between_mouse_var, var_between_cages),
    names_to = "variance_type",
    values_to = "variance"
  ) %>%
  mutate(
    variance_type = recode(
      variance_type,
      "within_mouse_var" = "Within host",
      "between_mouse_var" = "Between host,\nwithin cage",
      "var_between_cages" = "Between host,\nbetween cage"
    ),
    variance_type = factor(
      variance_type,
      levels = c(
        "Between host,\nbetween cage",
        "Between host,\nwithin cage",
        "Within host"
      )
    )
  )

# Step 6. Species titles

species_variances_long <- species_variances_long %>%
  mutate(
    taxonomy = sapply(species, function(species_oi) {
      parts <- str_split(species_oi, "_", simplify = TRUE)
      if (parts[2] == "sp") {
        paste0(parts[1], " sp.")
      } else {
        paste0(parts[1], " ", parts[2])
      }
    })
  )


# Step 7. Plot

p <- ggplot(species_variances_long, aes(x = taxonomy, y = variance, fill = variance_type)) +
  geom_col(position = position_dodge(width = 0.9)) +   # side-by-side bars
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = text_size),
    axis.text.x = element_text(angle = 45, hjust = 1, size = subtext_size),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = subtext_size),
    axis.title.y = element_text(size = text_size),
    title = element_text(size = text_size),
    legend.title = element_blank(),
    legend.text = element_text(size=subtext_size),
    plot.margin = unit(c(0.5, 0.5, 0.5, 2.5), "lines"),
    
  ) +
  labs(
    x = "Species",
    y = "Variance in strain frequency",
    fill = "Variance type"
  ) +
  scale_fill_manual(values = c("Between host,\nbetween cage" = "#AA98A9", 
                                "Between host,\nwithin cage" = "#BDB5D5", 
                                "Within host" = "#CC0000"))

p

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/HM_strain_variance.png"
ggsave(out_path, p, dpi = 300, height = 4, width = 7)









# Libraries
library(dplyr)
library(tidyr)


# Functions

relabel_samples <- function(string) {
  if (string == "TL1gDNAshort") {
    return("Inoculum")
  } else {
    mouse_no <- substr(string, 2, 2)
    
    if (substr(string, 3, 3) == "J") {
      return(paste0("Mouse ", mouse_no, " Jejunum"))
    } else if (substr(string, 3, 3) == "D") {
      return(paste0("Mouse ", mouse_no, " Duodenum"))
    } else if (substr(string, 3, 3) == "I") {
      return(paste0("Mouse ", mouse_no, " Ileum"))
    } else if (substr(string, 3, 4) == "Ce") {
      return(paste0("Mouse ", mouse_no, " Cecum"))
    } else {
      return(paste0("Mouse ", mouse_no, " Colon"))
    }
  }
}

# Data

# Directories

study <- "HumanizedMouse_Batch2" 
code_dir <- "~/Wasney-Briscoe/scripts/figures/"
metadata_dir <- "~/Wasney-Briscoe/metadata/"
data_dir <-  "~/"

pi_df <- read.csv(paste0(data_dir, "popgen_stats/SinglePi_SchloissnigPi_cov4_SiteDownsampled.csv")) %>%
  mutate(genomewide_pi = as.numeric(genomewide_pi), 
         total_loci = as.numeric(total_loci)) %>%
  filter(total_loci > 5e5) %>%
  rowwise() %>%
  mutate(sample_labels = relabel_samples(sample))


# Filters

inoculum_species <- pi_df %>% # Check that this is this case
  filter(sample_labels == "Inoculum") %>%
  select(species) %>%
  pull() %>%
  unique()

prevalent_species <- pi_df %>%
  filter(species %in% inoculum_species, sample_labels != "Inoculum") %>%
  group_by(species) %>%
  filter(n_distinct(mouse) >= 2, n() >= 3) %>%
  ungroup() %>%
  select(species) %>%
  pull() %>%
  unique()

# Filter

pi_df <- pi_df %>%
  filter(species %in% prevalent_species)

pi_wide <- pi_df %>%
  pivot_wider(
    names_from = sample_labels,      # This will create column names from "sample_labels"
    values_from = genomewide_pi,  # This will fill the new columns with "allele_frequency"
    id_cols = c(species)  # These columns will be kept as identifiers
  )

pi_wide['Mouse 2 Duodenum'] <- NA
pi_wide['Mouse 3 Ileum'] <- NA

supplementary_table_2 <- pi_wide %>%
  select("species", 
         "Inoculum",
         "Mouse 1 Duodenum","Mouse 1 Jejunum","Mouse 1 Ileum","Mouse 1 Cecum","Mouse 1 Colon",
         "Mouse 2 Duodenum","Mouse 2 Jejunum","Mouse 2 Ileum","Mouse 2 Cecum","Mouse 2 Colon",
         "Mouse 3 Duodenum","Mouse 3 Jejunum","Mouse 3 Ileum","Mouse 3 Cecum","Mouse 3 Colon",
         "Mouse 4 Duodenum","Mouse 4 Jejunum","Mouse 4 Ileum","Mouse 4 Cecum","Mouse 4 Colon",
         "Mouse 5 Duodenum","Mouse 5 Jejunum","Mouse 5 Ileum","Mouse 5 Cecum","Mouse 5 Colon",
         "Mouse 6 Duodenum","Mouse 6 Jejunum","Mouse 6 Ileum","Mouse 6 Cecum","Mouse 6 Colon",
         "Mouse 7 Duodenum","Mouse 7 Jejunum","Mouse 7 Ileum","Mouse 7 Cecum","Mouse 7 Colon",
         "Mouse 8 Duodenum","Mouse 8 Jejunum","Mouse 8 Ileum","Mouse 8 Cecum","Mouse 8 Colon")

out_path <- "~/supplementary_table_2.tsv"
write.table(supplementary_table_2, out_path, quote = FALSE, row.names = FALSE, sep = "\t")



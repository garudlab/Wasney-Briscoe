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

study <- "HumanizedMouse_Batch2" 
code_dir <- "~/Wasney-Briscoe/scripts/figures/"
metadata_dir <- "~/Wasney-Briscoe/metadata/"
data_dir <-  "~/"

snv_changes_df_path <- paste0(data_dir,"evolutionary_changes/snp_changes.txt.bz2")
snv_freqs_df_path <- paste0(data_dir,"evolutionary_changes/SNV_freqs.txt")

snv_freqs_df <- read.csv2(snv_freqs_df_path, sep = "\t", row.names = 1) %>%
  mutate(depth = as.numeric(depth),
         alt = as.numeric(alt),
         allele_frequency = as.numeric(allele_frequency)) %>%
  filter(depth >= 20)

snv_changes_df <- read.csv2(snv_changes_df_path, sep = ",")

# label with variant type

snv_freqs_df <- snv_freqs_df %>%
  left_join(snv_changes_df %>% select(species, contig, site_pos, variant_type) %>% unique(), by = c("species", "contig", "site_pos")) %>%
  filter(species != "Anaerostipes_hadrus_55206")

# Longform -> wideform

snv_freqs_df <- snv_freqs_df %>%
  rowwise() %>%
  mutate(sample_labels = relabel_samples(sample))

snv_freqs_wide <- snv_freqs_df %>%
  pivot_wider(
    names_from = sample_labels,      # This will create column names from "sample_labels"
    values_from = allele_frequency,  # This will fill the new columns with "allele_frequency"
    id_cols = c(species, contig, site_pos, gene, gene_description, variant_type)  # These columns will be kept as identifiers
  )


snv_freqs_wide["Mouse 2 Duodenum"] <- NA
snv_freqs_wide["Mouse 3 Duodenum"] <- NA
snv_freqs_wide["Mouse 3 Ileum"] <- NA

supplementary_table_3 <- snv_freqs_wide %>%
  select("species","contig","site_pos","gene","gene_description","variant_type",
         "Inoculum",
         "Mouse 1 Duodenum","Mouse 1 Jejunum","Mouse 1 Ileum","Mouse 1 Cecum","Mouse 1 Colon",
         "Mouse 2 Duodenum","Mouse 2 Jejunum","Mouse 2 Ileum","Mouse 2 Cecum","Mouse 2 Colon",
         "Mouse 3 Duodenum","Mouse 3 Jejunum","Mouse 3 Ileum","Mouse 3 Cecum","Mouse 3 Colon",
         "Mouse 4 Duodenum","Mouse 4 Jejunum","Mouse 4 Ileum","Mouse 4 Cecum","Mouse 4 Colon",
         "Mouse 5 Duodenum","Mouse 5 Jejunum","Mouse 5 Ileum","Mouse 5 Cecum","Mouse 5 Colon",
         "Mouse 6 Duodenum","Mouse 6 Jejunum","Mouse 6 Ileum","Mouse 6 Cecum","Mouse 6 Colon",
         "Mouse 7 Duodenum","Mouse 7 Jejunum","Mouse 7 Ileum","Mouse 7 Cecum","Mouse 7 Colon",
         "Mouse 8 Duodenum","Mouse 8 Jejunum","Mouse 8 Ileum","Mouse 8 Cecum","Mouse 8 Colon")




out_path <- "~/supplementary_table_3.tsv"
write.table(supplementary_table_3, out_path, quote = FALSE, row.names = FALSE, sep = "\t")




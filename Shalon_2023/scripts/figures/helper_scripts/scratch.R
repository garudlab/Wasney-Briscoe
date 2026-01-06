# Loading SNVs
snv_data_path <- "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/evolutionary_changes/SNV_freqs_WithinTimepoint.tsv"

snv_data <- read.csv2(snv_data_path, sep = "\t") %>% 
  rename(subject = subject_id, 
         timestamp = timepoint,
         sample_type_number = sample_type) %>%
  mutate(allele_frequency = as.numeric(allele_frequency),
         depth = as.numeric(depth)) %>%
  rowwise() %>%
  mutate(allele_frequency = ifelse(depth >= min_coverage, allele_frequency, NA))

snv_names <- snv_data %>%
  filter(species == species) %>%
  mutate(locus = paste0(contig, ", ", site_pos)) %>%
  select(locus) %>%
  unique() %>%
  pull()

snv_freqs <- snv_data %>%
  filter(species == species_oi,subject == subject_oi) %>%
  mutate(allele_frequency = as.numeric(allele_frequency))

all_snvs <- snv_freqs %>% select(contig, site_pos) %>% unique()



# AFTER STARTING TO MESS WITH WITH STRAIN FREQS
# Label with timepoint
snv_freqs <- snv_freqs %>%
  left_join(
    strain_freq_df %>% distinct(date, time, timepoint),
    by = c("date", "time")
  )

# Subsetting timepoints if necessary
# snvs
snv_freqs <- snv_freqs %>%
  select(-timepoint) %>%  
  left_join(
    strain_freq_df %>% distinct(date, time, timepoint),
    by = c("date", "time")
  ) %>%
  arrange(date, time)

## splitting duplicate names

snv_freqs <- snv_freqs %>%
  left_join(
    strain_freq_df %>% select(sample, sample_type_number) %>% distinct(),
    by = "sample"
  )
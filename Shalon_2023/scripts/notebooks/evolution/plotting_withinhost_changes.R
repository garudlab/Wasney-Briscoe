# Packages
library(dplyr)
library(tidyr)
library(ggplot2)

# Paths
input_path <- paste0("/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/evolutionary_changes/snv_frequencies.txt.bz2")

# Loading the data

allele_frequency_df <- read.csv2(input_path, sep = "\t")
allele_frequency_df$allele_frequency <- as.numeric(allele_frequency_df$allele_frequency)
wierd_sets <- as.character(5:20)
allele_frequency_df <- allele_frequency_df %>%
  rowwise() %>%
  # mutate(sample_set_label = ifelse(sample_set == "Stool", paste0("Stool, collected on ", collection_date), ifelse(sample_set %in% wierd_sets, paste0("Sample set ", sample_set, ", collected on ", swallow_date), paste0("Sample set ", sample_set))))
  mutate(sample_set_label = ifelse(sample_set == "Stool", paste0("Stool"), paste0("Sample set ", sample_set)))
  

##################################################################################################
##################################### INDIVIDUAL SPECIES/HOST ####################################
##################################################################################################

# Species
species_name <- "Bacteroides_vulgatus_57955"

allele_frequency_df_plotting <- allele_frequency_df %>%
  filter(species == species_name)

# Subjects
subjects <-  allele_frequency_df_plotting %>%
  select(subject) %>%
  unique() %>%
  pull()

subject_no <- 8

allele_frequency_df_plotting <- allele_frequency_df_plotting %>%
  filter(subject == subject_no)

# plotting

sample_sets <- allele_frequency_df_plotting %>% 
  select(sample_set_label, sample_set) %>%
  unique() %>%
  rowwise() %>%
  mutate(sample_set_numeric = ifelse(sample_set == "Stool", 1e6, as.numeric(sample_set))) %>%
  arrange(sample_set_numeric) %>%
  select(sample_set_label) %>%
  pull()


allele_frequency_df_plotting$location <- factor(allele_frequency_df_plotting$location, levels = c("Small intestine 1", "Small intestine 2", "Small intestine 3", "Ascending colon", "Stool"))
allele_frequency_df_plotting$sample_set_label <- factor(allele_frequency_df_plotting$sample_set_label, levels = sample_sets)


p <- ggplot(allele_frequency_df_plotting, aes(x = location, y = allele_frequency, color = paste0(contig, ", ",as.character(site_pos)), group = paste0(contig, ", ",as.character(site_pos)))) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none", 
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  facet_wrap(~sample_set_label, scale = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(color = "Genomic locus", x = NULL, y = "Allele frequency") +
  ggtitle(paste0(species_name, ", SNV changes along the gut"))
p

# Saving
output_path <- paste0("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/figures/evolutionary_changes/withinhost_changes/", species_name, "_withinhostchanges.png")
ggsave(output_path, p, dpi = 300)


##################################################################################################
############################################ FOR LOOP ############################################
##################################################################################################

species_with_adaptation <- allele_frequency_df %>%
  select(species) %>%
  unique() %>%
  pull()

for (species_name in species_with_adaptation) {
  print(paste0("Processing ", species_name))
  allele_frequency_df_plotting_species <- allele_frequency_df %>%
    filter(species == species_name)
  
  # Subjects
  subjects <-  allele_frequency_df_plotting_species %>%
    select(subject) %>%
    unique() %>%
    pull()
  
  for (subject_no in subjects) {
    print(paste0("Processing  ", species_name, " in ", subject_no))
    allele_frequency_df_plotting <- allele_frequency_df_plotting_species %>%
      filter(subject == subject_no)
    
    # plotting
    sample_sets <- allele_frequency_df_plotting %>% 
      select(sample_set_label, sample_set) %>%
      unique() %>%
      rowwise() %>%
      mutate(sample_set_numeric = ifelse(sample_set == "Stool", 1e6, as.numeric(sample_set))) %>%
      arrange(sample_set_numeric) %>%
      select(sample_set_label) %>%
      pull()
    
    allele_frequency_df_plotting$location <- factor(allele_frequency_df_plotting$location, levels = c("Small intestine 1", "Small intestine 2", "Small intestine 3", "Ascending colon", "Stool"))
    allele_frequency_df_plotting$sample_set_label <- factor(allele_frequency_df_plotting$sample_set_label, levels = sample_sets)
    
    
    p <- ggplot(allele_frequency_df_plotting, aes(x = location, y = allele_frequency, color = paste0(contig, ", ",as.character(site_pos)), group = paste0(contig, ", ",as.character(site_pos)))) +
      geom_line() +
      geom_point() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.position = "none", 
            strip.text = element_text(size = 8), 
            plot.title = element_text(size = 10),axis.title.y = element_text(size = 10)) +
      facet_wrap(~sample_set_label, scale = "free_x") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(color = "Genomic locus", x = NULL, y = "Allele frequency") +
      ggtitle(paste0(species_name, ", SNV changes along the gut in subject ", subject_no))
    
    print(paste0("Saving plot of ", species_name, " in subject ", subject_no))
    output_path <- paste0("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/figures/evolutionary_changes/withinhost_changes/", species_name, ", subject ", subject_no, "_withinhostchanges.png")
    ggsave(output_path, p, dpi = 300, height = 8, width = 8)
    
    
  }
}


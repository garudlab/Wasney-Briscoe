# Packages
library(dplyr)
library(tidyr)
library(ggplot2)

# Paths
species <- "Anaerostipes_hadrus_55206"
input_path <- paste0("/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/evolutionary_changes/", species, "_withinsetchanges.csv.bz2")

# Loading the data

allele_frequency_df <- read.csv2(input_path, sep = "\t")
allele_frequency_df$allele_frequency <- as.numeric(allele_frequency_df$allele_frequency)
allele_frequency_df$polarized_af <- as.numeric(allele_frequency_df$polarized_af)



# plotting
allele_frequency_df$location <- factor(allele_frequency_df$location, levels = c("Small intestine 1", "Small intestine 2", "Small intestine 3", "Ascending colon"))


p <- ggplot(allele_frequency_df %>% filter(sample_set == "3"), aes(x = location, y = polarized_af, color = paste0(contig, ", ",as.character(site_pos)), group = paste0(contig, ", ",as.character(site_pos)))) +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") +
  facet_wrap(~paste("Sample set",sample_set), scale = "free_x") +
  labs(color = "Genomic locus", x = NULL, y = "Allele frequency") +
  ggtitle(paste0(species, ", SNV changes along the gut"))
p

# Saving
output_path <- paste0("/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/figures/evolutionary_changes/", species, "_withinsetchanges.png")
ggsave(output_path, p, dpi = 300)




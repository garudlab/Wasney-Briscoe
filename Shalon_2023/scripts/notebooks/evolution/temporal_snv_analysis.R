# Packages
library(dplyr)
library(ggplot2)
library(cowplot)

# Loading data
## Paths
snps_df_path <- "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/evolutionary_changes/snp_changes.txt.bz2"
opportunities_df_path <- "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/evolutionary_changes/opportunities.txt.bz2"
## Loading
snps_df <- read.csv2(snps_df_path, sep = ",")
opportunities_df <- read.csv2(opportunities_df_path, sep = ",") 
opportunities_df$opportunities <- as.numeric(opportunities_df$opportunities)




# Filter opportunities dataframe to include comparisons we're interested in
allowed_sample_types <- c("Capsule 1", "Capsule 2", "Capsule 3", "Capsule 4")
allowed_sample_sets <- c("2", "3", "4", "5")

opportunities_df <- opportunities_df %>%
  filter(subject_1 == subject_2, 
         type_1 %in% allowed_sample_types, 
         type_2 %in% allowed_sample_types, 
         sample_set_1 %in% allowed_sample_sets,
         sample_set_2 %in% allowed_sample_sets)

opportunities_df$sample_set_1 <- as.numeric(opportunities_df$sample_set_1)
opportunities_df$sample_set_2 <- as.numeric(opportunities_df$sample_set_2)



# Functions

sum_snv_changes <- function(species_name, subject_id, sample_1, sample_2, snv_dataframe) {
  no_snvs <- nrow(snv_dataframe %>% 
                    filter(species == species_name, 
                           subject == subject_id, 
                           sample1 == sample_1, sample2 == sample_2))
  return(no_snvs)
}

set_to_time <- function(set) {
  if (set == 2) {
    return(0)
  } else if (set == 3) {
    return(6)
  } else if (set == 4) {
    return(24)
  } else if (set == 5) {
    return(30)
  } else {
    stop("Incorrect input to the set_to_time() function: must 2, 3, 4, or 5.")
  }
}

hour_difference <- function(set_1, set_2) {
  time_1 <- set_to_time(set_1)
  time_2 <- set_to_time(set_2)
  return(abs(time_1 - time_2))
}

# Annotation
## Annotate opportunities dataframe with number of SNVs

opportunities_df <- opportunities_df %>%
  rowwise() %>%
  mutate(no_changes = sum_snv_changes(species, subject_1, accession_1, accession_2, snps_df), 
         time_difference = hour_difference(sample_set_1, sample_set_2)) %>%
  mutate(change_rate = no_changes/opportunities)

opportunities_df$time_difference <- factor(opportunities_df$time_difference, levels = opportunities_df %>% select(time_difference) %>% unique() %>% arrange(time_difference) %>% pull())


# Plotting

p1 <- ggplot(opportunities_df, aes(x = time_difference, y = change_rate, color = time_difference)) +
  geom_jitter() +  # Add jittered points
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.8))) +
  ylim(0, (opportunities_df %>% select(change_rate) %>% max())*1.05) +
  labs(x = "Hours between sampling", y = "Rate of within-host SNV changes") +
  ggtitle("With Bacteroides vulgatus")

p2 <- ggplot(opportunities_df %>% filter(species != "Bacteroides_vulgatus_57955"), aes(x = time_difference, y = change_rate, color = time_difference)) +
  geom_jitter() +  # Add jittered points
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.8))) +
  ylim(0, (opportunities_df %>% select(change_rate) %>% max())*1.05) +
  labs(x = "Hours between sampling", y = "") +
  ggtitle("Without Bacteroides vulgatus")

p3 <- ggplot(opportunities_df, aes(x = time_difference, y = change_rate, color = time_difference)) +
  geom_jitter() +  # Add jittered points
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(0.8))) +
  ylim(0, 1e-5) +
  labs(x = "Hours between sampling", y = "Within-host SNV change rate") +
  ggtitle("With Bacteroides vulgatus")

p3



# Make a grid
grid <- ggdraw() +
  draw_plot(p1, 0, 0.0, 0.5, 0.94) +
  draw_plot(p2, 0.5, 0, 0.5, 0.94) +
  # draw_plot_label("Rate of within-host SNV changes over time", c(0.0), c(1),hjust = c(-0.1), size = 15, family = "Helvetica")
  draw_label("Rate of within-host SNV changes over time", fontface = "bold", x = 0,y = 1,hjust = 0, vjust=1)
grid


grid <- plot_grid(p1, p2, labels = "AUTO", align = "hv", rel_heights = c(1, 1)) + 
  draw_label("Rate of within-host SNV changes over time", fontface = "bold", x = 0,y = 1,hjust = 0, vjust=1)
  

plot_row <- plot_grid(p1, p2, labels = "AUTO")
title <- ggdraw() +
  draw_label("Rate of within-host SNV changes over time", fontface = "bold", x = 0,y = 1,hjust = 0, vjust=1.5) +
  theme(plot.margin = margin(0, 0, 0, 7))

grid <- plot_grid(
  title, plot_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.075, 0.925)
)

  
grid

out_file <- paste0("~/project-ngarud/Diversity-Along-Gut/Shalon_2023/figures/evolutionary_changes/temporal_snv_changes.png")
ggsave(out_file, grid, dpi = 300, bg = "white", height = 4, width = 10)



# Zoomed in (p3)
out_file <- paste0("~/project-ngarud/Diversity-Along-Gut/Shalon_2023/figures/evolutionary_changes/temporal_snv_changes_zoomed.png")
ggsave(out_file, p3, dpi = 300, bg = "white", height = 4, width = 10)








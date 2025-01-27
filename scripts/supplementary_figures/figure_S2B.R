library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")

text_size = 12
subtext_size = 10

# Paths
data_dir <-  "~/"
metadata_dir <- "~/Wasney-Briscoe-2024/metadata/"

### Haploid info df path
haploid_info_path <- paste0(metadata_dir, "sample_haploid_status.txt")

haploid_df <- read.csv2(haploid_info_path, sep = ",")

### Group by species and sum the number of (1) haploid and (2) polyploid samples in each grouping

haploid_counts_df <- haploid_df %>%
  group_by(species, ploid) %>%
  summarise(total_count = n())  %>%
  ungroup() %>%
  # filter(total_count > 1) %>%
  rowwise() %>%
  mutate(species = gsub("_", " ", species), 
         ploid = ifelse(ploid == "Haploid", "QP", "Not QP"))

species_order = haploid_df %>%
  group_by(species) %>%
  summarise(total_count = n()) %>%
  arrange(desc(total_count)) %>%
  rowwise() %>%
  mutate(species = gsub("_", " ", species)) %>%
  select(species) %>%
  pull()

### Plot
haploid_counts_df$species <- factor(haploid_counts_df$species, levels = rev(species_order))

haploid_counts_df_filter <- haploid_counts_df %>%
  group_by(species) %>%
  filter(sum(total_count) > 1) %>%
  ungroup()

haploid_counts_df_filter$species <- factor(haploid_counts_df_filter$species, levels = species_order)


p <- ggplot(haploid_counts_df_filter, aes(species, total_count, fill = ploid)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = subtext_size),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = subtext_size),
    legend.title = element_blank(),
    legend.text = element_text(size = subtext_size),
    legend.position = c(0.8, 0.8),
    plot.margin = margin(5, 15, 5, 65)
  ) +
  labs(y = "Number of high coverage samples") +
  scale_fill_manual(
    values = c("Not QP" = "light blue", "QP" = "dark blue")  # Set colors for "Not QP" and "QP"
  )

out_path <- "~/figure_S2B.png"
ggsave(outpath, p, dpi = 300, height = 4, width = 7.5, units = "in")
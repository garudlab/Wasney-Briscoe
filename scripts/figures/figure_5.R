
# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(scales)

# parameters

text_size = 12
subtext_size = 10

# Helper scripts

code_dir <- "~/Wasney-Briscoe-2024/scripts/figures/"
data_dir <-  "~/"
source(paste0(code_dir,"helper_scripts/SNV_plot_helper_scripts.R"))

# Paths

snp_data_path <- paste0(data_directory, "evolutionary_changes/snp_changes.txt.bz2")
gene_data_path <- paste0(data_directory, "evolutionary_changes/gene_changes.txt.bz2")
opportunities_path <- paste0(data_directory, "evolutionary_changes/opportunities.txt.bz2")


# Process data

aggregate_species = FALSE

if (aggregate_species) {
  difference_df <- snp_changes_df %>%
    filter(species != "Anaerostipes_hadrus_55206") %>%
    group_by(sample1, sample2) %>%
    summarise(no_of_changes = n(), .groups = "keep") %>%
    arrange(desc(no_of_changes)) %>%
    rename(sample_1 = sample1, sample_2 = sample2)
  
  opportunities_df <- opportunities_df %>%
    filter(species != "Anaerostipes_hadrus_55206") %>%
    group_by(sample_1, sample_2, host_orientation) %>%
    summarise(opportunities = sum(opportunities), .groups = "keep")
  
  difference_df <- merge(opportunities_df, difference_df, by = c("sample_1", "sample_2"), all.x = TRUE) %>%
    mutate(no_of_changes = ifelse(is.na(no_of_changes), 0, no_of_changes))
} else {
  difference_df <- snp_changes_df %>%
    filter(species != "Anaerostipes_hadrus_55206") %>%
    group_by(species, sample1, sample2) %>%
    summarise(no_of_changes = n(), .groups = "keep") %>%
    arrange(desc(no_of_changes)) %>%
    rename(sample_1 = sample1, sample_2 = sample2)
  
  opportunities_df <- opportunities_df %>%
    filter(species != "Anaerostipes_hadrus_55206")
  
  difference_df <- merge(opportunities_df, difference_df, by = c("species", "sample_1", "sample_2"), all.x = TRUE) %>%
    mutate(no_of_changes = ifelse(is.na(no_of_changes), 0, no_of_changes))
  
  
}



# Bootstrap

summary_of_group_differences <- difference_df %>%
  group_by(host_orientation) %>%
  summarise(no_of_pairs = n())
summary_of_group_differences

if (aggregate_species) {
  bootstrap_n = 50
  bootstrap_it = 1000
} else {
  bootstrap_n = 100
  bootstrap_it = 1000
}


within_changes_vec = c()
between_changes_vec = c()
inoc_changes_vec = c()

within_changes_perm_vec <- c()
between_changes_perm_vec <- c()
inoc_changes_perm_vec <- c()

test_statistics <- difference_df %>% 
  group_by(host_orientation) %>% 
  summarize(no_of_changes = sum(no_of_changes), opportunities = sum(opportunities)) %>% 
  mutate(rate_of_change = no_of_changes/opportunities)

set.seed(43)

for (it in 1:bootstrap_it) {
  if ((it %% 100) == 0) {
    print(paste0(it, " iterations completed"))
  }
  
  within_changes <- difference_df %>% 
    filter(host_orientation == "Within host") %>% 
    slice_sample(n = bootstrap_n, replace = TRUE)
  
  between_changes <- difference_df %>% 
    filter(host_orientation == "Between host") %>%
    slice_sample(n = bootstrap_n, replace = TRUE)
  
  inoc_changes <- difference_df %>% 
    filter(host_orientation == "Between inoculum") %>% 
    slice_sample(n = bootstrap_n, replace = TRUE)
  
  within_rate_of_change <- (within_changes$no_of_changes %>% sum())/(within_changes$opportunities %>% sum())
  between_rate_of_change <- (between_changes$no_of_changes %>% sum())/(between_changes$opportunities %>% sum())
  inoc_rate_of_change <- (inoc_changes$no_of_changes %>% sum())/(inoc_changes$opportunities %>% sum())
  
  within_changes_vec <- c(within_changes_vec, within_rate_of_change)
  between_changes_vec <- c(between_changes_vec, between_rate_of_change)
  inoc_changes_vec <- c(inoc_changes_vec, inoc_rate_of_change)
  
}

bootstrap_df <- data.frame(within_changes = within_changes_vec, 
                           between_changes = between_changes_vec, 
                           inoc_changes = inoc_changes_vec) %>%
  pivot_longer(cols = c(within_changes, 
                        between_changes, 
                        inoc_changes), 
               names_to = "change_type", 
               values_to = "value") %>%
  mutate(permuted = ifelse(grepl("perm", change_type), "Permuted data", "Real data")) %>%
  mutate(change_type = ifelse(grepl("between_changes", change_type), "Between host", change_type)) %>%
  mutate(change_type = ifelse(grepl("within_changes", change_type), "Within host", change_type)) %>%
  mutate(change_type = ifelse(grepl("inoc_changes", change_type), "Inoculum vs. host", change_type)) 


# Plotting

calc_CI <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

CIs <- data.frame(CI_lower = c(unname(calc_CI(within_changes_vec)['ymin']), 
                               unname(calc_CI(between_changes_vec)['ymin']),
                               unname(calc_CI(inoc_changes_vec)['ymin'])), 
                  CI_upper = c(unname(calc_CI(within_changes_vec)['ymax']), 
                               unname(calc_CI(between_changes_vec)['ymax']),
                               unname(calc_CI(inoc_changes_vec)['ymax'])),
                  value = c(unname(calc_CI(within_changes_vec)['middle']),
                            unname(calc_CI(between_changes_vec)['middle']),
                            unname(calc_CI(inoc_changes_vec)['middle'])),
                  change_type = c("Within host", "Between host", "Inoculum vs. host"))


bootstrap_df$change_type <- factor(bootstrap_df$change_type, levels = c("Within host", "Between host", "Inoculum vs. host"))

summarized_data <- bootstrap_df %>%
  group_by(change_type) %>%
  summarize(
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    median = median(value)
  )

summarized_data <- merge(summarized_data, test_statistics %>% rowwise() %>% mutate(host_orientation = ifelse(host_orientation == "Between inoculum", "Inoculum vs. host", host_orientation)), by.x = "change_type", by.y = "host_orientation")

summarized_data$change_type <- factor(summarized_data$change_type, levels = c("Within host", "Between host", "Inoculum vs. host"))

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

A <- ggplot(summarized_data, aes(x = change_type, y = rate_of_change)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.8),
                width = 0.2, size = 2) +
  geom_point(position = position_dodge(width = 0.8), size = 4, color = "red") +
  theme_bw() +
  labs(title = "Rate of SNV changes between samples", 
       y = "SNV change rate", 
       x = "Change Type") +
  theme(text = element_text(family = "Helvetica", size = text_size),
        axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  scale_color_aaas() +
  scale_y_continuous(label=scientific_10)

# SNV plots #
source("/u/project/ngarud/michaelw/Diversity-Along-Gut/HumanizedMouse/figures/figures/helper_scripts/SNV_plot_helper_scripts.R")

P_distasonis <- plot_snvs("Parabacteroides_distasonis_56985")
R_intestinalis <- plot_snvs("Roseburia_intestinalis_56239")
Coprococcus_sp <- plot_snvs("Coprococcus_sp_62244")
Ruminococcus_sp <- plot_snvs("Ruminococcus_sp_58571")

P_distasonis_legend <- plot_snvs("Parabacteroides_distasonis_56985", output = "legend")
R_intestinalis_legend <- plot_snvs("Roseburia_intestinalis_56239", output = "legend")
Coprococcus_sp_legend <- plot_snvs("Coprococcus_sp_62244", output = "legend")
Ruminococcus_sp_legend <- plot_snvs("Ruminococcus_sp_58571", output = "legend")


y_grid_interval <- 1/6
# strain_fig_offset <- 0.012
strain_fig_offset <- 0.015
# legend_offset <- 1/12
legend_offset <- 1/15


figure_5 <- ggdraw() +
  draw_plot(Coprococcus_sp, strain_fig_offset, y_grid_interval*5, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(P_distasonis, strain_fig_offset, y_grid_interval*4, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(R_intestinalis, strain_fig_offset, y_grid_interval*3, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(Ruminococcus_sp, strain_fig_offset, y_grid_interval*2, 1-strain_fig_offset - legend_offset, y_grid_interval) +
  draw_plot(Coprococcus_sp_legend, 1-legend_offset, y_grid_interval*5, legend_offset, y_grid_interval) +
  draw_plot(P_distasonis_legend, 1-legend_offset, y_grid_interval*4, legend_offset, y_grid_interval) +
  draw_plot(R_intestinalis_legend, 1-legend_offset, y_grid_interval*3, legend_offset, y_grid_interval) +
  draw_plot(Ruminococcus_sp_legend, 1-legend_offset, y_grid_interval*2, legend_offset, y_grid_interval) +
  draw_plot(A, 0, 0, 1, y_grid_interval*2) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1,  y_grid_interval*2), size = 16, family = "Helvetica") +
  draw_label("Frequency", x = 0.009, y = y_grid_interval*2 + (1-y_grid_interval*2)/2, angle = 90, size = 12, fontfamily = "Helvetica")

# figure_5

out_file <- "~/figure_5.png"
ggsave(out_file, figure_5, bg = "white", dpi = 300, width = 7.5, height = 10, units = "in")





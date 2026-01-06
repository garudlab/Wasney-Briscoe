library(dplyr)
library(ggplot2)

# Loading data
in_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/tables/snp_changes.tsv"
snp_changes_df <- read.csv2(in_path, sep = "\t")
in_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/tables/snp_change_rate.tsv"
change_rate_df <- read.csv2(in_path, sep = "\t")
change_rate_df <- change_rate_df %>%
  mutate(number_of_changes = as.numeric(number_of_changes),
         opportunities = as.numeric(opportunities),
         change_rate = as.numeric(change_rate))

# Bootstrapping

summary_of_group_differences <- change_rate_df %>%
  group_by(orientation) %>%
  summarise(no_of_pairs = n())
summary_of_group_differences

bootstrap_n = 100
bootstrap_it = 1000

within_changes_vec = c()
between_changes_vec = c()

within_changes_perm_vec <- c()
between_changes_perm_vec <- c()

test_statistics <- change_rate_df %>% 
  group_by(orientation) %>% 
  summarize(number_of_changes = sum(number_of_changes), opportunities = sum(opportunities)) %>% 
  mutate(rate_of_change = number_of_changes/opportunities)
test_statistics

set.seed(44)

for (it in 1:bootstrap_it) {
  if ((it %% 100) == 0) {
    print(paste0(it, " iterations completed"))
  }
  
  within_changes <- change_rate_df %>% 
    filter(orientation == "Within host") %>% 
    slice_sample(n = bootstrap_n, replace = TRUE)
  
  between_changes <- change_rate_df %>% 
    filter(orientation == "Between host") %>%
    slice_sample(n = bootstrap_n, replace = TRUE)
  
  within_rate_of_change <- (within_changes$number_of_changes %>% sum())/(within_changes$opportunities %>% sum())
  between_rate_of_change <- (between_changes$number_of_changes %>% sum())/(between_changes$opportunities %>% sum())

  within_changes_vec <- c(within_changes_vec, within_rate_of_change)
  between_changes_vec <- c(between_changes_vec, between_rate_of_change)

}

bootstrap_df <- data.frame(within_changes = within_changes_vec, 
                           between_changes = between_changes_vec) %>%
  pivot_longer(cols = c(within_changes, 
                        between_changes), 
               names_to = "change_type", 
               values_to = "value") %>%
  mutate(change_type = ifelse(grepl("between_", change_type), "Between host", "Within host")) 

calc_CI <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

CIs <- data.frame(CI_lower = c(unname(calc_CI(within_changes_vec)['ymin']), 
                               unname(calc_CI(between_changes_vec)['ymin'])), 
                  CI_upper = c(unname(calc_CI(within_changes_vec)['ymax']), 
                               unname(calc_CI(between_changes_vec)['ymax'])),
                  value = c(unname(calc_CI(within_changes_vec)['middle']),
                            unname(calc_CI(between_changes_vec)['middle'])),
                  change_type = c("Within host", "Between host"))

bootstrap_df$change_type <- factor(bootstrap_df$change_type, levels = c("Within host", "Between host"))

summarized_data <- bootstrap_df %>%
  group_by(change_type) %>%
  summarize(
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    median = median(value)
  )

summarized_data <- merge(summarized_data, test_statistics, by.x = "change_type", by.y = "orientation")

summarized_data$change_type <- factor(summarized_data$change_type, levels = c("Within host", "Between host", "Inoculum vs. host"))

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# parameters

text_size = 12
subtext_size = 10

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

A

out_path = "/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/figures/evolution/snv_change_rate_bootstrapped.png"
ggsave(out_path, A, dpi = 300, width = 6, height = 5)




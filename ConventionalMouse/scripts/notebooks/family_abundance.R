######################################################################
# FAMILY BARPLOT                                                     #
######################################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# Parameters
text_size = 12
subtext_size = 10

# Functions

source("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/helper_scripts/family_plot_helpers.R")

# Loading data
## Relative abundance

relab_df_path <- "/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/merge/species/species_relative_abundance.tsv"
relab_df <-  read.csv2(relab_df_path, sep = "\t", row.names = 1)
sample_names <- names(relab_df)
relab_df <-  read.csv2(relab_df_path, sep = "\t")
relab_df[sample_names] <- lapply(relab_df[sample_names], function(x) as.numeric(as.character(x)))
## Not all columns add to 1
relab_df[sample_names] <- sweep(relab_df[sample_names], 2, colSums(relab_df[sample_names]), FUN = "/")

## Species metadata

species_metadata_path <- "/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc/metadata.tsv"
species_metadata <- read.csv2(species_metadata_path, sep = "\t")


# annotating with taxonomy information

relab_df <- relab_df %>%
  rowwise() %>%
  mutate(family = extract_family(species_metadata[species_metadata["species_id"] == species_id, "taxonomy"]),
         order = extract_order(species_metadata[species_metadata["species_id"] == species_id, "taxonomy"]),
         order_family = paste0(order, "; ", family))

# Pivot to long form

relab_df_plot <- relab_df %>%
  pivot_longer(
    cols = -c(species_id, family, order, order_family),
    names_to = "sample_id",
    values_to = "relative_abundance"
  ) %>% 
  mutate(relative_abundance = as.numeric(relative_abundance))

# Annotate with cage and mouse

relab_df_plot <- relab_df_plot %>%
  rowwise() %>%
  mutate(cage = paste0("Cage ",str_split(sample_id, "_")[[1]][2]),
         mouse = paste0("Mouse ", as.character(as.numeric(str_split(sample_id, "_")[[1]][3]) + (as.numeric(str_split(sample_id, "_")[[1]][2])-1)*3)),
         mouse_number = as.numeric(str_split(sample_id, "_")[[1]][3]) + (as.numeric(str_split(sample_id, "_")[[1]][2])-1)*3,
         location = extract_location(sample_id))

# create legend

View(relab_df_plot %>% select(order_family) %>% arrange(order_family) %>% unique())

get_colors_order <- function(order) {
  Actinomycetales_colors <- c('Bifidobacteriaceae' = '#fcf4dd')
  
  Bacillales_colors <- c('Bacillaceae_G' = '#fce5cd')
  
  Bacteroidales_colors <- c('Muribaculaceae' = '#fb6f92',
                            'Rikenellaceae' = '#fffbff')
  
  Christensenellales_colors <- c('Borkfalkiaceae' = '#c0e0de',
                                 'CAG-552' = '#a0d6d3',
                                 'UBA3700' = '#7fc9c5')
  
  Clostridiales_colors <- c('Clostridiaceae' = '#9cadce')
  
  Coriobacteriales_colors <- c('Eggerthellaceae' = '#7ec4cf')
  
  Erysipelotrichales_colors <- c('Erysipelatoclostridiaceae' = '#d9b29e',
                                 'Erysipelotrichaceae' = '#d3ab9e')
  
  Haloplasmatales_colors <- c('Turicibacteraceae' = '#f0d6e0')
  
  Lachnospirales_colors <- c('Anaerotignaceae' = '#e7b6c6',
                             'CAG-274' = '#d9a1b8',
                             'Lachnospiraceae' = '#d291bc')
  
  Lactobacillales_colors <- c('Enterococcaceae' = '#caffbf',
                              'Lactobacillaceae' = '#e9ff70')
  
  Monoglobales_A_colors <- c('UBA1381' = '#f0e5a0')
  
  Oscillospirales_colors <- c('Acutalibacteraceae' = '#e8dff5',
                              'Butyricicoccaceae' = '#d1cfe2',
                              'Oscillospiraceae' = '#cdb4db',
                              'Ruminococcaceae' = '#ffc6ff')
  
  Peptostreptococcales_colors <- c('Anaerovoracaceae' = '#f0c0a0',
                                   'Peptostreptococcaceae' = '#e29980')
  
  RF39_colors <- c('CAG-1000' = '#c0d6f0')
  
  Staphylococcales_colors <- c('Staphylococcaceae' = '#f5b0b0')
  
  TANB77_colors <- c('CAG-508' = '#d0c0f0')
  
  switch(order,
         "Actinomycetales" = Actinomycetales_colors,
         "Bacillales" = Bacillales_colors,
         "Bacteroidales" = Bacteroidales_colors,
         "Christensenellales" = Christensenellales_colors,
         "Clostridiales" = Clostridiales_colors,
         "Coriobacteriales" = Coriobacteriales_colors,
         "Erysipelotrichales" = Erysipelotrichales_colors,
         "Haloplasmatales" = Haloplasmatales_colors,
         "Lachnospirales" = Lachnospirales_colors,
         "Lactobacillales" = Lactobacillales_colors,
         "Monoglobales_A" = Monoglobales_A_colors,
         "Oscillospirales" = Oscillospirales_colors,
         "Peptostreptococcales" = Peptostreptococcales_colors,
         "RF39" = RF39_colors,
         "Staphylococcales" = Staphylococcales_colors,
         "TANB77" = TANB77_colors,
         stop(paste("Error: Unknown order input:", order))
  )
}



library(ggnewscale)

all_orders <- relab_df_plot %>%
  arrange(order) %>%
  pull(order) %>%
  unique()

p1 <- ggplot(mapping = aes(x = mouse_number, y = 100*relative_abundance)) 
p2 <- ggplot(mapping = aes(x = mouse_number, y = 100*relative_abundance))
p3 <- ggplot(mapping = aes(x = mouse_number, y = 100*relative_abundance))
p4 <- ggplot(mapping = aes(x = mouse_number, y = 100*relative_abundance))

order_no <- 0
cutoff_1 <- 6 #7
cutoff_2 <- 10 #7
cutoff_3 <- 13 #10
for (order_oi in all_orders) {
  order_no <- order_no + 1
  print(paste0("Processing ", order_oi, " (", order_no, "/", length(relab_df_plot$order %>% unique()), ")"))
  colors_to_use <- get_colors_order(order_oi)
  if (order_no == 1) {
    p1 <- p1 + 
      geom_bar(data = filter(relab_df_plot, order == order_oi),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
  } else if (order_no < cutoff_1) {
    p1 <- p1 + 
      new_scale("fill") +
      geom_bar(data = filter(relab_df_plot, order == order_oi),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
    
    
  } else if (order_no == cutoff_1) {
    p2 <- p2 + 
      geom_bar(data = filter(relab_df_plot, order == order_oi),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
    
  } else if (order_no < cutoff_2) {
    p2 <- p2 + 
      new_scale("fill") +
      geom_bar(data = filter(relab_df_plot, order == order_oi),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
  } else if (order_no == cutoff_2) {
    p3 <- p3 + 
      geom_bar(data = filter(relab_df_plot, order == order_oi),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
  } else if (order_no < cutoff_3) {
    p3 <- p3 + 
      new_scale("fill") +
      geom_bar(data = filter(relab_df_plot, order == order_oi),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
  } else if (order_no == cutoff_3) {
    p4 <- p4 + 
      geom_bar(data = filter(relab_df_plot, order == order_oi),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))    
  } else {
    p4 <- p4 + 
      new_scale("fill") +
      geom_bar(data = filter(relab_df_plot, order == order_oi),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))     
  }
} 

p1 <- p1 +
  theme(text = element_text(size = subtext_size, family = "Helvetica"),
        # legend.position = "bottom",
        legend.key.size = unit(0.8, 'lines'),
        legend.direction = "vertical",
        legend.margin=margin(c(0,0,0,0)),
        legend.spacing.x = unit(0.75, "line"))
p2 <- p2 +
  theme(text = element_text(size = subtext_size, family = "Helvetica"),
        # legend.position = "bottom",
        legend.key.size = unit(0.8, 'lines'),
        legend.direction = "vertical",
        legend.margin=margin(c(0,0,0,0)),
        legend.spacing.x = unit(0.75, "line"))

p3 <- p3 +
  theme(text = element_text(size = subtext_size, family = "Helvetica"),
        # legend.position = "bottom",
        legend.key.size = unit(0.8, 'lines'),
        legend.direction = "vertical",
        legend.margin=margin(c(0,0,0,0)),
        legend.spacing.x = unit(0.75, "line"))

p4 <- p4 +
  theme(text = element_text(size = subtext_size, family = "Helvetica"),
        # legend.position = "bottom",
        legend.key.size = unit(0.8, 'lines'),
        legend.direction = "vertical",
        legend.margin=margin(c(0,0,0,0)),
        legend.spacing.x = unit(0.75, "line"))

legend_1 <- as_ggplot(get_legend(p1))

legend_2 <- as_ggplot(get_legend(p2))

legend_3 <- as_ggplot(get_legend(p3))

legend_4 <- as_ggplot(get_legend(p4))



# Plot

## Main panel

location_order = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")
relab_df_plot$location <- factor(relab_df_plot$location, levels = location_order)
relab_df_plot$mouse_number <- factor(relab_df_plot$mouse_number)

p <- ggplot(relab_df_plot, aes(x = mouse_number, y = 100*relative_abundance, fill = order_family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ location, nrow = 1, ncol = 5) +
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = text_size), #Change size
    axis.text = element_text(size = subtext_size),
    axis.title.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "white"),
    legend.position = "none"
  ) +
  labs(
    x = "Mouse number",
    y = "Relative Abundance (%)"
  ) + 
  scale_fill_manual(values = family_order_colors, name = "Taxonomic Family")

# p

## plot grid

proportion_legend <- 1/3
x_proportion_legend <- 1/4
y_proportion_legend <- 2/5
plot_offset_y <- 1/30



B <- ggdraw() + 
  draw_plot(p,0,y_proportion_legend + plot_offset_y, 1, 1-(y_proportion_legend + plot_offset_y)) +
  draw_plot(legend_1, 0.0, 0.0, 1/4, y_proportion_legend) +
  draw_plot(legend_2, 1/4, 0.0, 1/4,y_proportion_legend) +
  draw_plot(legend_3, 2/4, 0.0, 1/4,y_proportion_legend) +
  draw_plot(legend_4, 3/4, 0.0, 1/4, y_proportion_legend) +
  draw_label("Mouse", x = 0.5, y = y_proportion_legend + plot_offset_y/2, vjust = 0, size = text_size, angle = 0, fontfamily = "Helvetica")

# grid

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/figures/ecology/family_barplot.png"
ggsave(out_path,grid,dpi=300,width=6.5,height=7.5, units = "in", bg="white")


######################################################################
# alpha diversity                                                    #
######################################################################

# rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(vegan)

# Parameters
text_size = 12
subtext_size = 10

# Functions

source("/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/helper_scripts/family_plot_helpers.R")

# Loading data
## Relative abundance

relab_df_path <- "/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/merge/species/species_relative_abundance.tsv"
relab_df <-  read.csv2(relab_df_path, sep = "\t", row.names = 1)
sample_names <- names(relab_df)
relab_df[sample_names] <- lapply(relab_df[sample_names], function(x) as.numeric(as.character(x)))
## Not all columns add to 1
relab_df[sample_names] <- sweep(relab_df[sample_names], 2, colSums(relab_df[sample_names]), FUN = "/")


# Calculate alpha diversity 

shannon_div = diversity(relab_df, index = "shannon", MARGIN = 2)
to_plot =  data.frame(Sample.Name = names(shannon_div), ShannonDiversity = shannon_div)

# annotate plotting df

to_plot <- to_plot %>%
  rowwise() %>%
  mutate(location = extract_location(Sample.Name),
         mouse = extract_mouse(Sample.Name),
         cage = extract_cage(Sample.Name))

# Plot
to_plot$location <- factor(to_plot$location,levels = c("Duodenum","Jejunum","Ileum","Cecum","Colon", "Inoculum"))

my_comparisons <- list(
  c("Cecum", "Ileum"),
  c("Cecum", "Jejunum"),
  c("Cecum", "Duodenum"),
  c("Colon", "Ileum"),
  c("Colon", "Jejunum"),
  c("Colon", "Duodenum")
)


position_jittered <- position_jitter(width = 0.25, height = 0)
to_plot$ShannonDiversity_jittered <- jitter(to_plot$ShannonDiversity, amount = 0.25)

A <- ggplot(to_plot, aes(x=location,y=ShannonDiversity, label = mouse)) +
  geom_boxplot() +
  geom_jitter(size = 3,alpha = 0.5, position=position_jittered) +
  theme_bw() +
  theme(text = element_text(size=text_size, family="Helvetica"), 
        axis.text.x = element_text(size = subtext_size),
        axis.text.y = element_text(size = subtext_size),
        axis.title.x = element_blank(),
        legend.text = element_text(size = subtext_size),
        # legend.text = element_text(size = 12),
        legend.position = "right") + 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif",
                     vjust=0.2,hide.ns = TRUE,
                     method.args = list(alternative = "greater"), paired=TRUE) +
  labs(y = "Shannon diversity") 
A 


out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/figures/ecology/shannon_div.png"
ggsave(out_path,A,dpi=300,width=6,height=5, units = "in", bg="white")


######################################################################
# alpha diversity                                                    #
######################################################################

supp_fig_11 <- ggdraw() + 
  draw_plot(A, 0,3/4,1,1/4) + 
  draw_plot(B, 0,0/4,1,3/4) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 3/4), size = 16, family = "Helvetica")

out_path <- "/u/project/ngarud/michaelw/Diversity-Along-Gut/revisions/figures/CM_supp_11.png"
ggsave(out_path,supp_fig_11,width=7.5,height=10, units = "in", bg="white")






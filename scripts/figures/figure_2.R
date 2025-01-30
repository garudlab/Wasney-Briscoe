# Packages
require(vegan)
require("ggrepel")
require(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
require(ggpubr)
require(ggfortify)
library(stringr) 
require(stats)
library(factoextra)
library(ggsci)
library(hrbrthemes)
library(gridExtra)
library(cowplot)

# parameters

text_size = 12
subtext_size = 10

# Directories

study <- "HumanizedMouse_Batch2" 
code_dir <- "~/Wasney-Briscoe/scripts/figures/"
metadata_dir <- "~/Wasney-Briscoe/metadata/"
data_dir <-  "~/"


source(paste0(code_dir,"helper_scripts/humanized_mouse_utilities.R"))

# FIGURE 2A. SHANNON DIVERSITY #################################################

# Read in the data
relab_data = read.csv(paste0(data_dir,"merged_data/species/relative_abundance.txt.bz2"),sep="\t",stringsAsFactors = TRUE,row.names = 1)

metadata = read.csv(paste0(metadata_dir,"metadata.csv"),sep="\t",stringsAsFactors = FALSE)

row.names(metadata) = metadata$accession
colnames(relab_data) = metadata[chartr(".", "-", (colnames(relab_data))), "sample_id"]
shannon_div = diversity(relab_data, index = "shannon", MARGIN = 2)
to_plot =  data.frame(Sample.Name = names(shannon_div), ShannonDiversity = shannon_div)
head(to_plot)
head(shannon_div)
intersect(metadata$sample_id, to_plot$Sample.Name)
to_plot = merge(to_plot,metadata, by.x = "Sample.Name", by.y = "sample_id")
#shannon_div %>% tibble::rownames_to_column(var = "Sample") %>% pivot_longer()
#to_plot = data %>% tibble::rownames_to_column(var = "Source")%>% pivot_longer(-Source )
to_plot$TissueGroup = paste(to_plot$tissue_type,to_plot$Diet)

to_plot <- to_plot %>%
  rowwise() %>%
  mutate(tissue_type = ifelse(accession == "TL1gDNAshort", "Inoculum", tissue_type), 
         Diet = ifelse(accession == "TL1gDNAshort", "Inoculum", Diet),
         subject_id = ifelse(subject_id == 0, "Inoculum", as.character(subject_id)))


to_plot$tissue_type <- factor(to_plot$tissue_type,levels = c("Duodenum","Jejunum","Ileum","Cecum","Colon", "Inoculum"))

my_comparisons = inter_tissue_comparisons()
my_comparisons = my_comparisons[c(2,3,4,5,6,7)]

position_jittered <- position_jitter(width = 0.25, height = 0)
to_plot$ShannonDiversity_jittered <- jitter(to_plot$ShannonDiversity, amount = 0.25)

A <- ggplot(to_plot, aes(x=tissue_type,y=ShannonDiversity, label = subject_id)) +
  geom_boxplot() +
  geom_jitter(size = 3,alpha = 0.5, position=position_jittered) +
  theme_bw() +
  theme(text = element_text(size=text_size, family="Helvetica"), 
        axis.text.x = element_text(angle= 45,hjust=1, size = subtext_size),
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

# FIGURE 2B. FAMILY BARPLOT #################################################

# Helper script
source( paste0(code_dir, "helper_scripts/family_plot_helpers.R"))

# constant variabes

study = "HumanizedMouse_Batch2"

### Load the data
data <- read.csv(paste0(data_dir,"merged_data/species/relative_abundance.txt.bz2"),sep="\t",stringsAsFactors = TRUE)
metadata <- read.csv(paste0(metadata_dir,"metadata.csv"),sep="\t",stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(tissue_type = ifelse(accession == "TL1gDNAshort", "Inoculum", tissue_type))
taxonomy <- read.csv(paste0(code_dir,"genome_taxonomy.txt"),sep="\t",header=TRUE)
taxInfo <- read.csv(paste0(code_dir,"genome_info.txt"),sep="\t",header = TRUE)


### Clean the data
data_sub = data[,2:ncol(data)]
rownames(data_sub) = data$species_id
range(data_sub)
species_sufficient = names(which(rowSums(data_sub > 0.001) > 0))
data_sub = data_sub[species_sufficient,]
data_sub$species_id = row.names(data_sub)
data = data_sub

### processing taxonomy dataframes

taxInfo  = taxInfo %>% filter(rep_genome > 0)
taxInfo = taxInfo %>% distinct(genome_name, species_id)

to_plot = data %>% pivot_longer(!species_id,names_to = "sample", values_to = "abundance")
to_plot$accession = to_plot$sample
to_plot <- to_plot %>% 
  mutate(accession = chartr(".", "-", accession))
intersect(to_plot$accession, metadata$accession)
to_plot = merge(to_plot, metadata, by = "accession")
nrow(to_plot)
to_plot = merge(to_plot, taxInfo, by = "species_id")
nrow(to_plot)

taxonomy_distinct = taxonomy %>% distinct(genome_name, phylum, class, order, family, genus)
to_plot = merge(to_plot,taxonomy_distinct, by = "genome_name")
to_plot$tissue_type = factor(to_plot$tissue_type,levels = c("Inoculum", "Duodenum","Jejunum","Ileum","Cecum","Colon"))

## Make stacked barplots
to_plot_ = to_plot %>% group_by(subject_id,family, order, tissue_type) %>%
  summarize(TotalAbundance = sum(abundance,na.rm=TRUE))
to_plot_stats = to_plot_ %>% group_by(order,family) %>% 
  summarize(MeanAbundance = mean(TotalAbundance,na.rm=TRUE))

#length(unique(to_plotstats$family))
length(unique(to_plot$family))

# to_plot_$subject_id = ifelse(to_plot_$subject_id != 0 ,paste0("Mouse ", to_plot_$subject_id), "Inoculum")
to_plot_$subject_id = ifelse(to_plot_$subject_id != 0 ,paste0(to_plot_$subject_id), "Inoculum")
sortie = to_plot_stats %>% arrange(desc(MeanAbundance))
to_plot_ = to_plot_ %>% arrange(order, family)
to_plot_$order[to_plot_$order == ""] = "Unidentified"
to_plot_$family[to_plot_$family == ""] = "unidentified"
to_plot_$family_order = paste0(to_plot_$order, "; ",to_plot_$family)
# to_plot_$family_order = factor(to_plot_$family_order,levels = unique((to_plot_ %>% arrange(family_order))$family_order %>% unique()))
to_plot_$order = factor(to_plot_$order, levels = (to_plot_$order %>% unique() %>% sort()))

# adjusting inoculum bar

inoculum_factor <- to_plot_ %>% filter(tissue_type == "Inoculum") %>% 
  pull(TotalAbundance) %>% 
  sum()

to_plot_ <- to_plot_ %>%
  rowwise() %>%
  mutate(TotalAbundance = ifelse(tissue_type == "Inoculum", TotalAbundance*(1/inoculum_factor), TotalAbundance),
         subject_id = ifelse(subject_id == "TL1", "Inoculum", subject_id))

## legends

library(ggnewscale)

all_orders <- to_plot_ %>%
  arrange(ifelse(order == "Unidentified", 1, 0), order) %>%
  pull(order) %>%
  unique()

p1 <- ggplot(mapping = aes(x = subject_id, y = 100*TotalAbundance)) 
p2 <- ggplot(mapping = aes(x = subject_id, y = 100*TotalAbundance))
p3 <- ggplot(mapping = aes(x = subject_id, y = 100*TotalAbundance))
p4 <- ggplot(mapping = aes(x = subject_id, y = 100*TotalAbundance))


order_no <- 0
cutoff_1 <- 3 #7
cutoff_2 <- 5 #7
cutoff_3 <- 9 #10
for (order_oi in all_orders) {
  order_no <- order_no + 1
  print(paste0("Processing ", order_oi, " (", order_no, "/", length(to_plot_$order %>% unique()), ")"))
  no_of_colors_needed <- length((to_plot_ %>% filter(order == order_oi))$family %>% unique()) - 1
  colors_to_use <- get_colors_order(order_oi)
  if (order_no == 1) {
    p1 <- p1 + 
      geom_bar(data = filter(to_plot_, order == order_oi, subject_id != "Inoculum"),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
  } else if (order_no < cutoff_1) {
    p1 <- p1 + 
      new_scale("fill") +
      geom_bar(data = filter(to_plot_, order == order_oi, subject_id != "Inoculum"),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
    
    
  } else if (order_no == cutoff_1) {
    p2 <- p2 + 
      geom_bar(data = filter(to_plot_, order == order_oi, subject_id != "Inoculum"),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
    
  } else if (order_no < cutoff_2) {
    p2 <- p2 + 
      new_scale("fill") +
      geom_bar(data = filter(to_plot_, order == order_oi, subject_id != "Inoculum"),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
  } else if (order_no == cutoff_2) {
    p3 <- p3 + 
      geom_bar(data = filter(to_plot_, order == order_oi, subject_id != "Inoculum"),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
  } else if (order_no < cutoff_3) {
    p3 <- p3 + 
      new_scale("fill") +
      geom_bar(data = filter(to_plot_, order == order_oi, subject_id != "Inoculum"),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))
  } else if (order_no == cutoff_3) {
    p4 <- p4 + 
      geom_bar(data = filter(to_plot_, order == order_oi, subject_id != "Inoculum"),
               aes(fill = family),
               stat = "identity", 
               colour = "black") + 
      scale_fill_manual(order_oi, values = colors_to_use, guide = guide_legend(order = order_no))    
  } else {
    p4 <- p4 + 
      new_scale("fill") +
      geom_bar(data = filter(to_plot_, order == order_oi, subject_id != "Inoculum"),
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

###


C.1 = ggplot(as.data.frame(to_plot_) %>% filter(tissue_type != "Inoculum"), aes(x = subject_id, fill = family_order, y = 100*TotalAbundance)) +
  geom_bar(stat = "identity", colour = "black") + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = text_size), #Change size
        legend.position = "none",
        axis.text = element_text(size = subtext_size),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect(color = "black", fill = "white") 
  ) +
  scale_y_continuous(expand = c(0,0),limits = c(0,100.0000001)) + 
  labs(x = "Mouse", y = "Relative Abundance (%)", fill = "Taxonomic Order") + 
  facet_wrap(~tissue_type, nrow = 1,scale = "free") +
  scale_fill_manual(values = family_order_colors, name = "Taxonomic Family") +
  guides(fill = guide_legend(title = "Order; Family", ncol = 1))

# C.1

C.2 <- ggplot(to_plot_ %>% filter(tissue_type == "Inoculum") %>% mutate(tissue_type = "In", subject_id = "In"), aes(x = subject_id, fill = family_order, y = 100*TotalAbundance)) +
  geom_bar(stat = "identity", colour = "black") + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica", size = text_size),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = subtext_size),
        strip.background = element_rect(color = "black", fill = "white")
  ) +
  scale_y_continuous(expand = c(0,0),limits = c(0,100.0001)) + 
  labs(x = "Mouse", y = "Relative Abundance (%)", fill = "Taxonomic Order") + 
  facet_wrap(~tissue_type, nrow = 1,scale = "free") +
  scale_fill_manual(values = family_order_colors, name = "Taxonomic Family") +
  guides(fill = guide_legend(title = NULL))

# C.2
proportion_inoc <- 0.11
# proportion_legend <- 1/3
x_proportion_legend <- 1/4
y_proportion_legend <- 1/3
plot_offset_y <- 1/30


C <- ggdraw() +
  draw_plot(C.2, 0.0, y_proportion_legend + plot_offset_y, proportion_inoc, 1 - (y_proportion_legend + plot_offset_y)) +
  draw_plot(C.1, 0.0 + proportion_inoc, y_proportion_legend + plot_offset_y, 1 - proportion_inoc, 1 - (y_proportion_legend + plot_offset_y)) +
  draw_plot(legend_1, 0.0, 0.0, 1/4, y_proportion_legend) +
  draw_plot(legend_2, 1/4, 0.0, 1/4, y_proportion_legend) +
  draw_plot(legend_3, 2/4, 0.0, 1/4, y_proportion_legend) +
  draw_plot(legend_4, 3/4, 0.0, 1/4, y_proportion_legend) +
  draw_label("Mouse", x = 0.5, y = y_proportion_legend + plot_offset_y/2, vjust = 0, size = text_size, angle = 0, fontfamily = "Helvetica") +
  draw_line(x = c(0, 1), y = c(y_proportion_legend, y_proportion_legend), color = "black", size = 1)

# FINAL GRID #

legend_adjustment <- 0.05
c_size <- 0.6
prop_c <- 0.6

figure_2 <- ggdraw() +
  draw_plot(A, 0,prop_c,1,1-prop_c) +
  draw_plot(C, 0,0,1,prop_c) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, prop_c), size = 16, family = "Helvetica")


out_file <- "~/figure_2.png"
ggsave(out_file,figure_2,dpi=300,width=7.5,height=10, units = "in", bg="white")


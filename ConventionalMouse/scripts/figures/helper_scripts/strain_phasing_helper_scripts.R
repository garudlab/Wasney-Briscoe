# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(scales)
library(ggpubr)
library(grid)
library(ggtext)


plot_strains <- function(species, output = "normal", strain_colors = c("#654321", "#D2B48C"), species_metadata = "/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc/metadata.tsv") {
  
  print(paste0("Processing ", species, "."))
  
  # Setting parameters
  
  font_size <- 12
  subtext_size <- 10
  x_axis_font_size <- 8
  
  # Loading data
  
  strain_trajectory_path <- paste0("/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/strain_phasing/strain_clusters/", 
                                   species, 
                                   "/",
                                   species, 
                                   "_strain_frequency.csv")
  
  
  #Loading strain frequencies
  strain_trajectory_df <- read.csv2(strain_trajectory_path, sep = "\t")
  strain_trajectory_df$freq <- as.numeric(strain_trajectory_df$freq)
  strain_trajectory_df$quantile_25 <- as.numeric(strain_trajectory_df$quantile_25)
  strain_trajectory_df$quantile_75 <- as.numeric(strain_trajectory_df$quantile_75)
  strain_trajectory_df$upper_ci <- as.numeric(strain_trajectory_df$upper_ci)
  strain_trajectory_df$lower_ci <- as.numeric(strain_trajectory_df$lower_ci)
  strain_trajectory_df$species <- as.character(strain_trajectory_df$species)
  strain_trajectory_df$strain <- as.character(strain_trajectory_df$strain)
  strain_trajectory_df <- strain_trajectory_df %>%
    mutate(mouse_id = gsub(" ", "_", paste(cage, mouse_number)))
  
  # Getting species names
  species_metadata_path <- "/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc/metadata.tsv"
  species_metadata_df <- read.csv2(species_metadata, sep = "\t") 
  species_metadata_df$species_id <- as.character(species_metadata_df$species_id) 
  
  taxonomy <- str_split(species_metadata_df[species_metadata_df['species_id'] == species,"species"], "__")[[1]][2]
  
  if (length(str_split(taxonomy, " ")[[1]]) == 1) {
    taxonomy <- paste0(taxonomy, " sp.")
  }
  
  
  ## Changing mouse number
  
  change_mouse_number <- function(mouse_id) {
    cage <- as.numeric(str_split(mouse_id, "_")[[1]][2])
    mouse <- as.numeric(str_split(mouse_id, "_")[[1]][4])
    new_mouse_number = (cage - 1)*3 + mouse
    mouse_label = paste0("Mouse ", as.character(new_mouse_number))
    return(mouse_label)
  }
  
  strain_trajectory_df$mouse_label <- sapply(strain_trajectory_df$mouse_id, change_mouse_number)
  
  ## Identify strain order
  
  mice_present <- as.vector(strain_trajectory_df %>% filter(!is.na(freq)) %>% select(mouse_id) %>% unique() %>% pull())
  
  strain_order <- strain_trajectory_df %>% 
    group_by(mouse_id, strain) %>% 
    summarise(avg_freq = mean(freq, na.rm = TRUE)) %>% 
    group_by(strain) %>% 
    summarise(avg_strain_freq = mean(avg_freq)) %>% 
    arrange(desc(avg_strain_freq)) %>% 
    select(strain) %>% 
    pull() 
  
  #Shortening gut labels
  shorten_gut_label <- function(label) {
    if (is.na(label)) {
      return("Inoculum")
    } else if ((label == "Duodenum") | (label == "Jejunum") | (label == "Ileum")) {
      return(substr(label, 1, 1))
    } else if ((label == "Cecum") | (label == "Colon")) {
      return(substr(label, 1, 2))
    } else if (label == "Inoculum") {
      return(label)
    }
  }
  
  strain_trajectory_df$gut_label <- sapply(strain_trajectory_df$region, shorten_gut_label)
  
  # Add empty rows into mice in which the strain is absent
  
  if (!("Cage_1_Mouse_1" %in% strain_trajectory_df$mouse_id)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = "1", mouse_number = "Mouse 1", region = "Colon", cage = "Cage 1", sample = "Col_1_1", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_1_Mouse_1", mouse_label = "Mouse 1", gut_label = "Co")
  }
  if (!("Cage_1_Mouse_2" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = "1", mouse_number = "Mouse 2", region = "Colon", cage = "Cage 1", sample = "Col_1_2", freq = 0, quantile_25 = NA, quantile_75 = NA,upper_ci = NA, lower_ci = NA, mouse_id = "Cage_1_Mouse_2", mouse_label = "Mouse 2", gut_label = "Co")
  }
  if (!("Cage_1_Mouse_3" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = "1", mouse_number = "Mouse 3", region = "Colon", cage = "Cage 1", sample = "Col_1_3", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_1_Mouse_3", mouse_label = "Mouse 3", gut_label = "Co")
  }
  if (!("Cage_2_Mouse_1" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = "1", mouse_number = "Mouse 1", region = "Colon", cage = "Cage 2", sample = "Col_2_1", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_2_Mouse_1", mouse_label = "Mouse 4", gut_label = "Co")
  }
  if (!("Cage_2_Mouse_2" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = "1", mouse_number = "Mouse 2", region = "Colon", cage = "Cage 2", sample = "Col_2_2", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_2_Mouse_2", mouse_label = "Mouse 5", gut_label = "Co")
  }
  if (!("Cage_2_Mouse_3" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = "1", mouse_number = "Mouse 3", region = "Colon", cage = "Cage 2", sample = "Col_2_3", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_2_Mouse_3", mouse_label = "Mouse 6", gut_label = "Co")
  }
  if (!("Cage_3_Mouse_1" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = "1", mouse_number = "Mouse 1", region = "Colon", cage = "Cage 3", sample = "Col_3_1", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA,  mouse_id = "Cage_3_Mouse_1", mouse_label = "Mouse 7", gut_label = "Co")
  }
  if (!("Cage_3_Mouse_2" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = "1", mouse_number = "Mouse 2", region = "Colon", cage = "Cage 3", sample = "Col_3_2", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_3_Mouse_2", mouse_label = "Mouse 8", gut_label = "Co")
  }
  if (!("Cage_3_Mouse_3" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = "1",  mouse_number = "Mouse 3", region = "Colon", cage = "Cage 3", sample = "Col_3_3", freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, mouse_id = "Cage_3_Mouse_3", mouse_label = "Mouse 9", gut_label = "Co")
  }
  
  #Factoring
  gut_order = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")
  gut_label_order = c("D", "J", "I", "Ce", "Co")
  mouse_order = c("Cage_1_Mouse_1", "Cage_1_Mouse_2", "Cage_1_Mouse_3", "Cage_2_Mouse_1", "Cage_2_Mouse_2", "Cage_2_Mouse_3", "Cage_3_Mouse_1", "Cage_3_Mouse_2", "Cage_3_Mouse_3")
  strain_trajectory_df$region = factor(strain_trajectory_df$region, levels = gut_order)
  strain_trajectory_df$gut_label = factor(strain_trajectory_df$gut_label, levels = gut_label_order)
  strain_trajectory_df$mouse_id = factor(strain_trajectory_df$mouse_id, levels = mouse_order)
  strain_trajectory_df$strain = factor(strain_trajectory_df$strain, levels = strain_order)
  
  
  # Palette
  if (length(strain_order) == 3) {
    pal <- c("#404040","#808080","#C0C0C0" )
  } else if (length(strain_order) == 2) {
    pal <- strain_colors
    # pal <- c("#204552", "#7f7e7e")
  } else if (length(strain_order) == 1) {
    pal <- c("#654321") #brown
  }
  
  
  error_bars <- strain_trajectory_df %>% 
    dplyr::filter(strain == strain_order[1])
  
  # plot_title <- paste0(
  #   "<i>", taxonomy, "</i>, conventional mouse cohort"
  # )
  
  if ((str_split(taxonomy, " ")[[1]][2] == "sp.") | (str_split(taxonomy, " ")[[1]][2] == "sp000403395")) {
    plot_title <- paste0(
      "<i>", str_split(taxonomy, " ")[[1]][1], "</i>", " sp., conventional mice"
    )
  } else {
    plot_title <- paste0(
      "<i>", taxonomy, "</i>, conventional mice"
    )
  }
  
  
  
  if (output == "normal") {
    
    p <-  ggplot(strain_trajectory_df, aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      geom_errorbar(data = error_bars,
                    aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            plot.title = element_markdown(size = font_size),
            axis.title = element_text(size = font_size), 
            axis.title.y = element_text(size = font_size), 
            axis.title.x = element_blank(),
            # legend.position = "none",
            axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_text(size = x_axis_font_size),
            strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(color = "black", size = x_axis_font_size),
            # axis.ticks.y = element_blank(), 
            # strip.text = element_text(color = "white", size = 10),
            plot.background = element_rect(fill = "#FFFFFF80", color = NA),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt")
      ) +
      scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap("mouse_label", scales = "free_x", nrow = 1) +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal, labels = c("1", "2")) +
      labs(title = plot_title, fill = "Strain", y = "Strain frequency")
    
    
  } else if (output == "blank") {
    
    p <-  ggplot(strain_trajectory_df, aes(x = gut_label, y = freq, fill = strain)) +
      geom_blank() +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_text(size = font_size), 
            plot.title = element_markdown(size = font_size),
            # legend.position = "none",
            axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(),
            strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(color = "black", size = x_axis_font_size),
            axis.ticks.y = element_blank(), 
            # strip.text = element_text(color = "white", size = 10),
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef"),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt")
      ) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap("mouse_label", scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal) +
      labs(title = plot_title, fill = "Strain", y = "Strain frequency")
    
  } else {
    stop("You must provide a valid value for `output` (normal, legend, blank).")
  }
  
  if (output == "legend") {
    
    return(as_ggplot(legend))
  } else {
    # 
    # rect1 <- rectGrob(
    #   x = 0.0275, y = 0.0,            # center of panels 1-3
    #   width = 0.31, height = 0.96,
    #   gp = gpar(fill = "#e7e1ef", col = NA),
    #   just = c("left", "bottom")
    # )
    # 
    # rect2 <- rectGrob(
    #   x = 0.335, y = 0.0,            # center of panels 4-6
    #   width = 0.31, height = 0.96,
    #   gp = gpar(fill = "#AA98A9", col = NA),
    #   just = c("left", "bottom")
    # )
    # 
    # rect3 <- rectGrob(
    #   x = 0.6425, y = 0.0,            # center of panels 7-9
    #   width = 0.3075, height = 0.96,
    #   gp = gpar(fill = "#e7e1ef", col = NA),
    #   just = c("left", "bottom")
    # )
    # 
    # # Put them in a list (optional)
    # rects <- list(rect1, rect2, rect3)
    # 
    # grid <- ggdraw() +
    #   draw_grob(rect1) +
    #   draw_grob(rect2) +
    #   draw_grob(rect3) +
    #   draw_plot(p, 0, 0, 1, 1)   # main plot on top
    # 
    # return(grid)
    return(p)
    
  }
}


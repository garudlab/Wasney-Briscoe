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

font_size <- 12

data_dir <-  "~/"


plot_strains <- function(species, output = "normal", polarize_by_inoculum = TRUE, show_diet = FALSE, strain_colors = c("#654321", "#D2B48C")) {
  # output: "normal", "inoculum", "inoculum only", "blank" "legend"
  
  print(paste0("Processing ", species, "."))
  
  # Loading data
  
  strain_trajectory_path <- paste0(data_dir,"strain_phasing/strain_clusters/", 
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
  
  
  ## Identify strain order
  
  mice_present <- as.vector(strain_trajectory_df %>% filter(!is.na(freq)) %>% select(mouse_number) %>% unique() %>% pull())
  
  if ((polarize_by_inoculum) & ("Inoculum" %in% strain_trajectory_df$mouse_number)) {
    strain_order <- strain_trajectory_df %>% 
      filter(mouse_number == "Inoculum") %>% 
      arrange(desc(freq)) %>% 
      select(strain) %>% 
      pull()
  } else {
    strain_order <- strain_trajectory_df %>% 
      group_by(mouse_number, strain) %>% 
      summarise(avg_freq = mean(freq, na.rm = TRUE)) %>% 
      group_by(strain) %>% 
      summarise(avg_strain_freq = mean(avg_freq)) %>% 
      arrange(desc(avg_strain_freq)) %>% 
      select(strain) %>% 
      pull()
  }
  
  
  
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
  
  strain_trajectory_df <- strain_trajectory_df %>%
    rowwise() %>%
    mutate(gut_label = shorten_gut_label(region))
  
  if (!("1" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = 1, mouse_number = "1", region = "Colon", diet = "Control diet", cage = "Cage 1", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, gut_label = "Co")
  }
  if (!("2" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = 1, mouse_number = "2", region = "Colon", diet = "Control diet", cage = "Cage 1", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA,upper_ci = NA, lower_ci = NA, gut_label = "Co")
  }
  if (!("3" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = 1, mouse_number = "3", region = "Colon", diet = "Control diet", cage = "Cage 1", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, gut_label = "Co")
  }
  if (!("4" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = 1, mouse_number = "4", region = "Colon", diet = "Guar gum diet", cage = "Cage 2", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, gut_label = "Co")
  }
  if (!("5" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = 1, mouse_number = "5", region = "Colon", diet = "Guar gum diet", cage = "Cage 2", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, gut_label = "Co")
  }
  if (!("6" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = 1, mouse_number = "6", region = "Colon", diet = "Guar gum diet", cage = "Cage 3", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, gut_label = "Co")
  }
  if (!("7" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = 1, mouse_number = "7", region = "Colon", diet = "Guar gum diet", cage = "Cage 3", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, gut_label = "Co")
  }
  if (!("8" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = 1, mouse_number = "8", region = "Colon", diet = "Guar gum diet", cage = "Cage 3", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, gut_label = "Co")
  }
  if (!("Inoculum" %in% strain_trajectory_df$mouse_number)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species, strain = 1, mouse_number = "Inoculum", region = "Colon", diet = "Human diet", cage = "Inoculum", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, upper_ci = NA, lower_ci = NA, gut_label = "Co")
  }
  
  strain_trajectory_df <- as.data.frame(strain_trajectory_df)
  
  
  
  #Factoring
  gut_order = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")
  gut_label_order = c("D", "J", "I", "Ce", "Co")
  mouse_order = c("Inoculum", "1", "2", "3", "4", "5", "6", "7", "8")
  strain_trajectory_df$region = factor(strain_trajectory_df$region, levels = gut_order)
  strain_trajectory_df$gut_label = factor(strain_trajectory_df$gut_label, levels = gut_label_order)
  strain_trajectory_df$mouse_number = factor(strain_trajectory_df$mouse_number, levels = mouse_order)
  strain_trajectory_df$strain = factor(strain_trajectory_df$strain, levels = strain_order)
  
  # Palette
  if (length(strain_order) == 3) {
    # pal <- c("#C0C0C0","#808080", "#404040" )
    pal <- c("#404040","#808080","#C0C0C0" )
  } else if (length(strain_order) == 2) {
    # pal <- c("#404040", "#C0C0C0")
    # pal <- c("#4B0082", "#D8BFD8") #purple
    # pal <- c("#654321", "#D2B48C") #brown
    # pal <- c("#204552", "#7f7e7e")
    # pal <- c("#014421","#C1E1C1") #green
    # pal <- c("#996515","#FFD700") #brown.yellow
    # pal <- c("#00008B", "#ADD8E6")
    pal <- strain_colors
    
  } else if (length(strain_order) == 1) {
    pal <- c("#404040")
    pal <- c("#654321") #brown
  }
  
  # Plotting
  
  error_bars <- strain_trajectory_df %>% 
    filter(strain == strain_order[1])
  
  if ((output == "normal") | (output == "inoculum only") | (output == "inoculum only (blank)")) {

    x_axis_font_size <- 8
    
    if (output == "inoculum only") {
      inoculum <- ggplot(strain_trajectory_df %>% filter(mouse_number == "Inoculum"), aes(x = gut_label, y = freq, fill = strain)) +
        geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
        geom_errorbar(data = error_bars %>% filter(cage == "Inoculum"),
                      aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                      width = 0.2, color = "white") +
        theme_bw() +
        theme(text = element_text(family = "Helvetica"),
              axis.title = element_text(size = 12, face = "bold.italic"), 
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              # axis.text.x = element_text(size = 5, color = "#644117"),
              axis.text.x = element_text(size = x_axis_font_size, face = "bold"),
              # axis.text.x = element_text(size = 5),
              legend.position = "none",
              strip.background = element_rect(color = "black", fill = "white"), 
              strip.text = element_text(size = 10),
              plot.title = element_text(size = 8, face = "bold")) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        scale_x_discrete(labels = c("Inoc")) +
        facet_wrap(~substr(mouse_number,0,0), scale = "free_x") +
        #scale_fill_brewer(palette = pal)
        scale_fill_manual(values = pal) +
        labs(title = paste(strsplit(species, "_")[[1]][1], strsplit(species, "_")[[1]][2]))
      # inoculum
    } else if (output == "inoculum only (blank)") {
      inoculum <- ggplot(strain_trajectory_df %>% filter(mouse_number == "Inoculum"), aes(x = gut_label, y = freq, fill = strain)) +
        geom_blank() +
        theme_bw() +
        theme(text = element_text(family = "Helvetica"),
              # axis.title = element_blank(), 
              axis.title = element_text(size = 12, face = "bold.italic"), 
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              # axis.text.x = element_text(size = 5, color = "#644117"),
              axis.text.x = element_text(size = x_axis_font_size, face = "bold"),
              # axis.text.x = element_text(size = 5),
              legend.position = "none",
              strip.background = element_rect(color = "black", fill = "white"), 
              strip.text = element_text(size = 10),
              plot.title = element_text(size = 8, face = "bold")) +
        scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
        scale_x_discrete(labels = c("Inoc")) +
        facet_wrap(~substr(mouse_number,0,0), scale = "free_x") +
        #scale_fill_brewer(palette = pal)
        scale_fill_manual(values = pal) +
        labs(title = paste(strsplit(species, "_")[[1]][1], strsplit(species, "_")[[1]][2]))
      
    } else {
      
      color_names <- c("black", "red", "black", "black","red","black")
    
      
      inoculum <- ggplot(strain_trajectory_df %>% filter(mouse_number == "Inoculum"), aes(x = gut_label, y = freq, fill = strain)) +
        geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
        geom_errorbar(data = error_bars %>% filter(cage == "Inoculum"),
                      aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                      width = 0.2, color = "white") +        
        theme_bw() +
        theme(text = element_text(family = "Helvetica"),
              axis.title = element_blank(), 
              axis.text.x = element_text(size = x_axis_font_size, face = "bold"),
              axis.text.y = element_text(size = x_axis_font_size),
              legend.position = "none",
              strip.background = element_rect(color = "black", fill = "white"), 
              strip.text = element_text(size = font_size),
              plot.margin = margin(r = 2, l = 5 , t = 4, b = 4, unit = "pt")
              ) +
        scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8, 1)) + #, limits = c(0,1)
        coord_cartesian(ylim = c(0, 1)) +
        scale_x_discrete(labels = c("In")) +
        facet_wrap(~substr(mouse_number,0,0), scale = "free_x") +
        #scale_fill_brewer(palette = pal)
        scale_fill_manual(values = pal)
      # inoculum
    }
    
    
    
    cage_1 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 1"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      geom_errorbar(data = error_bars %>% filter(cage == "Cage 1"),
                    aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            legend.position = "none",
            axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(),
            strip.background = element_rect(fill = ifelse(show_diet, scales::alpha("gold", 0.75), "white"), color = "black"),
            # strip.background = element_rect(fill = scales::alpha("#08519c", 0.75), color = "black"),
            strip.text = element_text(color = "black", size = font_size),
            axis.ticks.y = element_blank(), 
            # strip.text = element_text(color = "white", size = 10),
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef"),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt")
            ) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal) +
      labs(fill = "Strain")
    # cage_1
    
    cage_2 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 2"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      geom_errorbar(data = error_bars %>% filter(cage == "Cage 2"),
                    aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      theme_bw() +
      theme(text = element_text(family = "Helvetica", size = font_size),
            axis.title = element_blank(), 
            #legend.key.width = unit(0.12, "npc"), 
            legend.position = "none",
            # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
            # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
            axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.background = element_rect(fill = ifelse(show_diet,scales::alpha("#08519c", 0.85), "white"), color = "black"),
            strip.text = element_text(color = ifelse(show_diet,"white","black"), size = font_size),
            # plot.background = element_rect(fill = "white", color = "black"),
            plot.background = element_rect(fill = "#AA98A9", color = "#AA98A9"),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt") 
      ) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal)
    # cage_2
    
    cage_3 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 3"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      geom_errorbar(data = error_bars %>% filter(cage == "Cage 3"),
                    aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      theme_bw() +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            legend.position = "none", 
            # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
            # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
            axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.background = element_rect(fill = ifelse(show_diet,scales::alpha("#08519c", 0.85), "white"), color = "black"),
            strip.text = element_text(color = ifelse(show_diet,"white","black"), size = font_size),
            # plot.background = element_rect(fill = "#e7e1ef", color = "black"),
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef"),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt")
            ) +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal)
    # cage_3
  } else if (output == "inoculum") {
    x_axis_font_size <- 7.5
    
    inoculum <- ggplot(strain_trajectory_df %>% filter(mouse_number == "Inoculum"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      geom_errorbar(data = error_bars %>% filter(cage == "Inoculum"),
                    aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            # axis.text.x = element_text(size = 5, color = "#644117"),
            # axis.text.x = element_text(size = x_axis_font_size, face = "bold"),
            # axis.text.x = element_text(size = 5),
            legend.position = "none",
            strip.background = element_rect(color = "black", fill = "white"), 
            strip.text = element_text(size = 10)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(labels = c("Inoc")) +
      facet_wrap(~substr(mouse_number,0,0), scale = "free_x") +
      #scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal)
    # inoculum
    
    cage_1 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 1"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_blank() +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            legend.position = "none",
            # legend.key.width = unit(0.1, "npc"),
            # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
            # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
            # axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.background = element_rect(fill = scales::alpha("gold", 0.75), color = "black"),
            # strip.background = element_rect(fill = scales::alpha("#08519c", 0.75), color = "black"),
            strip.text = element_text(color = "black", size = 10),
            # strip.text = element_text(color = "white", size = 10),
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef")) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal) +
      labs(fill = "Strain")
    # cage_1
    
    cage_2 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 2"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_blank() +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            #legend.key.width = unit(0.12, "npc"), 
            legend.position = "none",
            # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
            # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
            # axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.background = element_rect(fill = scales::alpha("#08519c", 0.85), color = "black"),
            strip.text = element_text(color = "white", size = 10)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal)
    # cage_2
    
    cage_3 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 3"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_blank() +
      theme_bw() +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            legend.position = "none", 
            # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
            # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
            # axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.background = element_rect(fill = scales::alpha("#08519c", 0.85), color = "black"),
            strip.text = element_text(color = "white", size = 10), 
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef")) +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal)
    # cage_3
  } else if (output == "blank") {
    x_axis_font_size <- 7.5
    
    inoculum <- ggplot(strain_trajectory_df %>% filter(mouse_number == "Inoculum"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_blank() +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            # axis.text.x = element_text(size = 5, color = "#644117"),
            axis.text.x = element_text(size = x_axis_font_size, face = "bold"),
            # axis.text.x = element_text(size = 5),
            legend.position = "none",
            strip.background = element_rect(color = "black", fill = "white"), 
            strip.text = element_text(size = 10)) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      scale_x_discrete(labels = c("Inoc")) +
      facet_wrap(~substr(mouse_number,0,0), scale = "free_x") +
      #scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal)
    # inoculum
    
    cage_1 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 1"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_blank() +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            legend.position = "none",
            # legend.key.width = unit(0.1, "npc"),
            # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
            # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
            axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.background = element_rect(fill = scales::alpha("gold", 0.75), color = "black"),
            # strip.background = element_rect(fill = scales::alpha("#08519c", 0.75), color = "black"),
            strip.text = element_text(color = "black", size = 10),
            # strip.text = element_text(color = "white", size = 10),
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef")) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal) +
      labs(fill = "Strain")
    # cage_1
    
    cage_2 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 2"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_blank() +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            #legend.key.width = unit(0.12, "npc"), 
            legend.position = "none",
            # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
            # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
            # axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.background = element_rect(fill = scales::alpha("#08519c", 0.85), color = "black"),
            strip.text = element_text(color = "white", size = 10)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal)
    # cage_2
    
    cage_3 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 3"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_blank() +
      theme_bw() +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            legend.position = "none", 
            # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
            # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
            # axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.background = element_rect(fill = scales::alpha("#08519c", 0.85), color = "black"),
            strip.text = element_text(color = "white", size = 10), 
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef")) +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal)
    # cage_3
  } else if (output == "legend") {
    x_axis_font_size <- 8
    
    strain_order <- strain_trajectory_df %>% 
      group_by(mouse_number, strain) %>% 
      summarise(avg_freq = mean(freq, na.rm = TRUE)) %>% 
      group_by(strain) %>% 
      summarise(avg_strain_freq = mean(avg_freq)) %>% 
      arrange(avg_strain_freq) %>% 
      select(strain) %>% 
      pull()
    
    strain_trajectory_df$strain = factor(strain_trajectory_df$strain, levels = strain_order)
    
    cage_3 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 3"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      theme_bw() +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse_number), scales = "free_x") +
      theme(text = element_text(family = "Helvetica", size = x_axis_font_size),
            axis.title = element_blank(), 
            # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
            # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
            axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.background = element_rect(fill = scales::alpha("#08519c", 0.85), color = "black"),
            strip.text = element_text(color = "white", size = 10), 
            legend.title = element_text(size = x_axis_font_size),
            legend.text = element_text(size = x_axis_font_size),
            # legend.key.size = unit(1, 'lines'),
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef")) +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal) +
      guides(fill = guide_legend(override.aes = list(size = 3))) 
    
  }
  
  
  
  if (output == "legend") {
    # legend <- get_legend(cage_3) 
    # grid.newpage() 
    legend <- get_legend(cage_3)
    return(as_ggplot(legend))
  } else if ((output == "inoculum only") | (output == "inoculum only (blank)")) {
    return(inoculum)
  } else {
    inoc_size <- 0.075 #0.0675
    height <- 0.88
    
    grid <- ggdraw() +
      draw_plot(inoculum, 0, 0.0, inoc_size, height) +
      draw_plot(cage_1, inoc_size, 0.0, (1-inoc_size)*(3/8), height) +
      draw_plot(cage_2, inoc_size+(1-inoc_size)*(3/8), 0.0, (1-inoc_size)*(2/8), height) +
      draw_plot(cage_3, inoc_size+(1-inoc_size)*(3/8)+(1-inoc_size)*(2/8), 0.0, (1-inoc_size)*(3/8), height) +
      draw_plot_label(c(paste(strsplit(species, "_")[[1]][1], strsplit(species, "_")[[1]][2])), c(0.045), c(1),hjust = 0,vjust = 1.25, size = font_size, family = "Helvetica", fontface = "italic")
    
    return(grid)
  }
}


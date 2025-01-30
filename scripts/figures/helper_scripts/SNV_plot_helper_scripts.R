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

data_dir <- "~/"


plot_snvs <- function(species_name, output = "normal",min_coverage = 20, snv_filter = FALSE, focal_contig = NA, focal_site_pos = NA, polarize_by_inoculum = TRUE) {
  # output: "normal", "legend"
  
  # Loading SNVs
  snv_data_path <- paste0(data_dir,
                          "evolutionary_changes/SNV_freqs.txt")

  snv_data <- read.csv2(snv_data_path, sep = "\t", row.names = "X") %>%
    mutate(allele_frequency = as.numeric(allele_frequency),
           alt = as.numeric(alt),
           depth = as.numeric(depth)) %>%
    rowwise() %>%
    mutate(allele_frequency = ifelse(depth >= min_coverage, allele_frequency, NA))
  
  species_with_adaptation <- snv_data %>% select(species) %>% unique() %>% pull()
  
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
  strain_trajectory_df <- rename(strain_trajectory_df, mouse = mouse_number)
  strain_trajectory_df$mouse <- as.character(strain_trajectory_df$mouse)
  strain_trajectory_df$upper_ci <- as.numeric(strain_trajectory_df$upper_ci)
  strain_trajectory_df$lower_ci <- as.numeric(strain_trajectory_df$lower_ci)
  
  ## Identify strain order
  
  mice_present <- as.vector(strain_trajectory_df %>% filter(!is.na(freq)) %>% select(mouse) %>% unique() %>% pull())
  
  if ((polarize_by_inoculum) & ("Inoculum" %in% strain_trajectory_df$mouse)) {
    strain_order <- strain_trajectory_df %>% 
      filter(mouse == "Inoculum") %>% 
      arrange(desc(freq)) %>% 
      select(strain) %>% 
      pull()
  } else {
    strain_order <- strain_trajectory_df %>% 
      group_by(mouse, strain) %>% 
      summarise(avg_freq = mean(freq, na.rm = TRUE)) %>% 
      group_by(strain) %>% 
      summarise(avg_strain_freq = mean(avg_freq)) %>% 
      arrange(desc(avg_strain_freq)) %>% 
      select(strain) %>% 
      pull()
  }
  
  # strain_order <- strain_trajectory_df %>% 
  #   group_by(mouse, strain) %>% 
  #   summarise(avg_freq = mean(freq, na.rm = TRUE)) %>% 
  #   group_by(strain) %>% 
  #   summarise(avg_strain_freq = mean(avg_freq)) %>% 
  #   arrange(desc(avg_strain_freq)) %>% 
  #   select(strain) %>% 
  #   pull()
  
  # View(strain_trajectory_df %>% group_by(sample) %>% summarise(strain_sum = sum(freq)))
  
  snv_species <- snv_data %>%
    filter(species == species_name)
  
  
  
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
  
  snv_species <- snv_species %>%
    rowwise() %>%
    mutate(gut_label = shorten_gut_label(gut_site))
  
  #Annotate strain frequencies with evolutionary changes
  snv_names <- snv_data %>%
    filter(species == species_name) %>%
    mutate(locus = paste0(contig, ", ", site_pos)) %>%
    select(locus) %>%
    unique() %>%
    pull()
  # 
  snv_freqs <- snv_data %>%
    filter(species == species_name) %>%
    select(sample, mouse, contig, site_pos, allele_frequency) %>%
    mutate(allele_frequency = as.numeric(allele_frequency))
  
  all_snvs <- snv_freqs %>% select(contig, site_pos) %>% unique()
  
  # Repolarizing SNVs
  for (snv_idx in 1:nrow(all_snvs)) {
    contig_oi <- all_snvs[snv_idx, "contig"]
    site_pos_oi <- all_snvs[snv_idx, "site_pos"]
    mice_present_adaptation <- snv_freqs %>% filter(contig == contig_oi, site_pos == site_pos_oi, !is.na(allele_frequency)) %>% select(mouse) %>% unique() %>% pull()
    mice_present_adaptation <- intersect(mice_present_adaptation, mice_present)
    print(paste0("Contig: ", contig_oi, "; site_pos: ", site_pos_oi))
    # sample_list <- snv_freqs %>% filter(!is.na(allele_frequency), contig == contig_oi, site_pos == site_pos_oi) %>% select(sample) %>% unique() %>% pull()
    sample_list <- strain_trajectory_df %>% select(mouse) %>% unique() %>% pull()
    adaptation_sample_list <- snv_freqs %>% filter(contig==contig_oi, site_pos==site_pos_oi, !is.na(allele_frequency)) %>% select(mouse) %>% unique() %>% pull()
    # if ((species_name == "Ruminococcus_sp_58571") & (contig_oi == "CAHL01000016") & (site_pos_oi == 6907)) {
    #   # next
    # } else 
    if (("Inoculum" %in% sample_list) & ("Inoculum" %in% adaptation_sample_list)) { #If the allele is present in the inoculum
      inoculum_af <- snv_freqs %>% filter(contig == contig_oi, site_pos == site_pos_oi, sample == "TL1gDNAshort") %>% select(allele_frequency) %>% pull()
      if (inoculum_af > 0.5) { #If the allele is over 0.5 AF in the inoculum
        print(paste0("Repolarizing SNV ", contig_oi, ", ", as.character(site_pos_oi), " according to inoculum."))
        snv_freqs <- snv_freqs %>%
          rowwise() %>%
          mutate(allele_frequency = ifelse((contig == contig_oi) & (site_pos == site_pos_oi), 1 - allele_frequency, allele_frequency))
      }
    } else if (sum(snv_freqs %>% filter(contig == contig_oi, site_pos == site_pos_oi) %>% filter(!is.na(allele_frequency)) %>% group_by(mouse) %>% summarise(avg_af = mean(allele_frequency, na.rm = TRUE)) %>% select(avg_af) > 0.5)/length((snv_freqs %>% filter(contig == contig_oi, site_pos == site_pos_oi) %>% filter(!is.na(allele_frequency)) %>% select(mouse) %>% unique() %>% pull())) >= 0.5) {
      if ((length(setdiff(mice_present_adaptation, "Inoculum")) <= 2)) {
        first_mouse <- as.character(sort(as.numeric(setdiff(mice_present_adaptation, "Inoculum")))[1])
        mean_af_in_first_mouse <- snv_freqs %>% 
          filter(mouse == first_mouse, contig == contig_oi, site_pos == site_pos_oi) %>% 
          group_by(mouse) %>%
          summarise(mean_af = mean(allele_frequency, na.rm = TRUE)) %>% 
          pull()
        if (mean_af_in_first_mouse > 0.5) {
          print(paste0("Repolarizing SNV ", contig_oi, ", ", as.character(site_pos_oi), " according to first mouse (one or two mice present)."))
          snv_freqs <- snv_freqs %>%
            rowwise() %>%
            mutate(allele_frequency = ifelse((contig == contig_oi) & (site_pos == site_pos_oi), 1 - allele_frequency, allele_frequency))
        }
      } else if (sum(snv_freqs %>% filter(contig == contig_oi, site_pos == site_pos_oi) %>% filter(!is.na(allele_frequency)) %>% group_by(mouse) %>% summarise(avg_af = mean(allele_frequency, na.rm = TRUE)) %>% select(avg_af) > 0.5)/length((snv_freqs %>% filter(contig == contig_oi, site_pos == site_pos_oi) %>% filter(!is.na(allele_frequency)) %>% select(mouse) %>% unique() %>% pull())) == 0.5) {
        first_mouse <- as.character(sort(as.numeric(setdiff(mice_present_adaptation, "Inoculum")))[1])
        mean_af_in_first_mouse <- snv_freqs %>% 
          filter(mouse == first_mouse, contig == contig_oi, site_pos == site_pos_oi) %>% 
          group_by(mouse) %>%
          summarise(mean_af = mean(allele_frequency, na.rm = TRUE)) %>% 
          pull()
        if (mean_af_in_first_mouse > 0.5) {
          print(paste0("Repolarizing SNV ", contig_oi, ", ", as.character(site_pos_oi), " according to first mouse (polarization evenly split between mice)."))
          snv_freqs <- snv_freqs %>%
            rowwise() %>%
            mutate(allele_frequency = ifelse((contig == contig_oi) & (site_pos == site_pos_oi), 1 - allele_frequency, allele_frequency))
        }
      } else if (sum(snv_freqs %>% filter(contig == contig_oi, site_pos == site_pos_oi) %>% group_by(mouse) %>% summarise(mean_af = mean(allele_frequency, na.rm = TRUE)) %>% select(mean_af) %>% na.omit() > 0.5)/length((snv_freqs %>% filter(contig == contig_oi, site_pos == site_pos_oi) %>% filter(!is.na(allele_frequency)) %>% select(mouse) %>% unique() %>% pull())) > 0.5) { #if more than half of mice have over 0.5 frequency on average, repolarize
        print(paste0("Repolarizing SNV ", contig_oi, ", ", as.character(site_pos_oi), " according to mean across mice."))
        snv_freqs <- snv_freqs %>%
          rowwise() %>%
          mutate(allele_frequency = ifelse((contig == contig_oi) & (site_pos == site_pos_oi), 1 - allele_frequency, allele_frequency))
      }
    }
  } 
  
  # ADD MISSING SAMPLES TO SNV_FREQS
  strain_samples <- strain_trajectory_df %>% select(sample) %>% unique() %>% pull()
  snv_samples <- snv_freqs %>% select(sample) %>% unique() %>% pull()
  non_shared_samples <- setdiff(strain_samples, snv_samples)
  samples_to_drop <- setdiff(snv_samples,strain_samples)
  snv_freqs <- snv_freqs %>%
    filter(!(sample %in% samples_to_drop))
  rows_to_add <- all_snvs %>% 
    slice(rep(1:n(), each = length(non_shared_samples))) 
  rows_to_add$sample <- rep(non_shared_samples, nrow(all_snvs))
  rows_to_add$allele_frequency <- NA
  snv_freqs <- rbind(snv_freqs %>% select(-mouse), rows_to_add)
  ### SOMETHING WRONG HERE
  
  strain_trajectory_df <- merge(snv_freqs, strain_trajectory_df, by = c("sample"), all.y = TRUE, all.x = TRUE)
  
  # Focus on one SNV?

  
  if (snv_filter) {
    strain_trajectory_df <- strain_trajectory_df %>%
      filter(contig %in% focal_contig, site_pos %in% focal_site_pos)
    
    snv_names <- strain_trajectory_df %>%
      filter(contig %in% focal_contig, site_pos %in% focal_site_pos) %>%
      mutate(locus = paste0(contig, ", ", site_pos)) %>%
      select(locus) %>%
      unique() %>%
      pull()
  }
  
  #Add empty rows for mice that aren't present
  if (!("1" %in% strain_trajectory_df$mouse)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species_name, strain = 1, mouse = "1", region = "Colon", diet = "Control diet", cage = "Cage 1", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, gut_label = "Co", contig = NA, site_pos = NA, allele_frequency = NA)
  }
  if (!("2" %in% strain_trajectory_df$mouse)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species_name, strain = 1, mouse = "2", region = "Colon", diet = "Control diet", cage = "Cage 1", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, gut_label = "Co", contig = NA, site_pos = NA, allele_frequency = NA)
  }
  if (!("3" %in% strain_trajectory_df$mouse)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species_name, strain = 1, mouse = "3", region = "Colon", diet = "Control diet", cage = "Cage 1", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, gut_label = "Co", contig = NA, site_pos = NA, allele_frequency = NA)
  }
  if (!("4" %in% strain_trajectory_df$mouse)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species_name, strain = 1, mouse = "4", region = "Colon", diet = "Guar gum diet", cage = "Cage 2", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, gut_label = "Co", contig = NA, site_pos = NA, allele_frequency = NA)
  }
  if (!("5" %in% strain_trajectory_df$mouse)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species_name, strain = 1, mouse = "5", region = "Colon", diet = "Guar gum diet", cage = "Cage 2", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, gut_label = "Co", contig = NA, site_pos = NA, allele_frequency = NA)
  }
  if (!("6" %in% strain_trajectory_df$mouse)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species_name, strain = 1, mouse = "6", region = "Colon", diet = "Guar gum diet", cage = "Cage 3", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA, gut_label = "Co", contig = NA, site_pos = NA, allele_frequency = NA)
  }
  if (!("7" %in% strain_trajectory_df$mouse)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species_name, strain = 1, mouse = "7", region = "Colon", diet = "Guar gum diet", cage = "Cage 3", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA,gut_label = "Co", contig = NA, site_pos = NA, allele_frequency = NA)
  }
  if (!("8" %in% strain_trajectory_df$mouse)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species_name, strain = 1, mouse = "8", region = "Colon", diet = "Guar gum diet", cage = "Cage 3", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA,gut_label = "Co",contig = NA, site_pos = NA, allele_frequency = NA)
  }
  if (!("Inoculum" %in% strain_trajectory_df$mouse)) {
    strain_trajectory_df <- strain_trajectory_df %>%
      add_row(species = species_name, strain = 1, mouse = "Inoculum", region = "Stool", diet = "Human diet", cage = "Inoculum", sample = NA, freq = 0, quantile_25 = NA, quantile_75 = NA,gut_label = "Co", contig = NA, site_pos = NA, allele_frequency = NA)
  }
  
  strain_trajectory_df <- as.data.frame(strain_trajectory_df) 
  strain_trajectory_df$locus <- ifelse(is.na(strain_trajectory_df$contig), NA, paste0(strain_trajectory_df$contig, ", ", strain_trajectory_df$site_pos))

  
  # if (any(is.na(strain_trajectory_df$contig))) {
  #   snv_pal <- c(snv_pal, "black")
  # }
  
  #Factoring
  gut_order = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Colon")
  gut_label_order = c("D", "J", "I", "Ce", "Co")
  mouse_order = c("Inoculum", "1", "2", "3", "4", "5", "6", "7", "8")
  adaptation_order <- strain_trajectory_df %>% select(locus) %>% unique() %>% pull()
  adaptation_order <- adaptation_order[!is.na(adaptation_order)]
  strain_trajectory_df$locus <- factor(strain_trajectory_df$locus, levels = adaptation_order)
  strain_trajectory_df$region = factor(strain_trajectory_df$region, levels = gut_order)
  strain_trajectory_df$gut_label = factor(strain_trajectory_df$gut_label, levels = gut_label_order)
  strain_trajectory_df$mouse = factor(strain_trajectory_df$mouse, levels = mouse_order)
  strain_trajectory_df$strain = factor(strain_trajectory_df$strain, levels = strain_order)
  
  # if (output != "legend") {
  #   strain_trajectory_df$strain = factor(strain_trajectory_df$strain, levels = strain_order)
  # }
  
  
  #Adjust for the number of SNVs
  strain_trajectory_df <- strain_trajectory_df %>% mutate(freq = freq/length(snv_names)) # 
  
  #Choosing palette
  
  if (length(strain_order) == 3) {
    # pal <- c("#C0C0C0","#808080", "#404040" )
    pal <- c("#404040","#808080","#C0C0C0" )
  } else if (length(strain_order) == 2) {
    pal <- c("#404040", "#C0C0C0")
    pal <- c("#654321", "#D2B48C") #brown
  } else if (length(strain_order) == 1) {
    pal <- c("#404040")
    pal <- c("#654321")
  }
  
  # snv_pal <- heat.colors(length(snv_names))
  # snv_pal <- rainbow(n = length(snv_names))
  if (length(snv_names) > 12) {
    snv_pal <- rainbow(n = length(snv_names))
  } else {
    snv_pal <- c('#ffff99','#e31a1c','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#a6cee3','#1f78b4','#33a02c','#cab2d6','#6a3d9a','#b15928')
  }
  # snv_pal <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2')
  # show_col(snv_pal)
  
  if ((snv_filter) & length(contig_oi == 1)) {
    # snv_pal <- c('#F93822')
    # snv_pal <- c("#ffdb58")
    # snv_pal <- c("#00ab41")
    snv_pal <- c("#FFAC1C")
    # snv_pal <- c('white')
  }
  
  
  
  ###############################################################################
  ############################## PLOT 1A. BARPLOT ###############################
  ###############################################################################
  
  error_bars <- strain_trajectory_df %>% 
    filter(strain == strain_order[1])
  
  x_axis_font_size <- 8

  snv_dot_size = 3
  snv_species$fill <- 1
  
  if (output == "normal") {
    color_names <- c("black", "red", "black", "black","red","black")
    inoculum <- ggplot(strain_trajectory_df %>% filter(mouse == "Inoculum"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      geom_errorbar(data = error_bars %>% filter(cage == "Inoculum"),
                    aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      geom_point(data = strain_trajectory_df %>% filter(mouse == "Inoculum", strain  == 1), aes(x = gut_label, y = allele_frequency, color = locus), size = snv_dot_size, show.legend = FALSE, shape = 8, position = position_jitter(width = 0.1, height = 0)) +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            axis.text.x = element_text(size = x_axis_font_size, face = "bold"),
            axis.text.y = element_text(size = x_axis_font_size, colour = color_names),
            legend.position = "none",
            strip.background = element_rect(color = "black", fill = "white"), 
            strip.text = element_text(size = font_size),
            plot.margin = margin(r = 2, l = 5 , t = 4, b = 4, unit = "pt")
            ) +
      scale_y_continuous(breaks = c(0, 0.2,0.4,0.6,0.8, 1), limits = c(0,1)) +
      scale_x_discrete(labels = c("In")) +
      facet_wrap(~substr(mouse,0,0), scale = "free_x") +
      #scale_fill_brewer(palette = pal)
      scale_color_manual(values = snv_pal) +
      scale_fill_manual(values = pal)
    # inoculum
    
    set.seed(100)
    
    cage_1 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 1"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      geom_errorbar(data = error_bars %>% filter(cage == "Cage 1"),
                    aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      geom_point(data = strain_trajectory_df %>% filter(cage == "Cage 1", strain  == 1), aes(y = allele_frequency, color = locus), size = snv_dot_size, show.legend = FALSE, shape = 8, position = position_jitter(width = 0.1, height = 0)) +
      theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.title = element_blank(), 
            legend.position = "none",
            axis.text.x = element_text(size = x_axis_font_size),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            # strip.background = element_rect(fill = "gold", color = "black"),
            # strip.background = element_rect(fill = scales::alpha("gold", 0.75), color = "black"),
            # strip.background = element_rect(fill = scales::alpha("#08519c", 0.85), color = "black"),
            strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(color = "black", size = font_size),
            # strip.text = element_text(color = "white", size = font_size), 
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef"),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt") 
            ) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse), scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal) +
      scale_color_manual(values = snv_pal)
    # labs(fill = "Strain")
    # cage_1
    
    cage_2 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 2"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
      geom_errorbar(data = error_bars %>% filter(cage == "Cage 2"),
                    aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      geom_point(data = strain_trajectory_df %>% filter(cage == "Cage 2", strain  == 1), aes(y = allele_frequency, color = locus), size = snv_dot_size, show.legend = FALSE, shape = 8, position = position_jitter(width = 0.1, height = 0)) +
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
            # strip.background = element_rect(fill = "gold", color = "black"),
            # strip.background = element_rect(fill = scales::alpha("gold", 0.75), color = "black"),
            # strip.background = element_rect(fill = scales::alpha("#08519c", 0.85), color = "black"),
            strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(color = "black", size = font_size),
            # strip.text = element_text(color = "white", size = font_size)
            plot.background = element_rect(fill = "#AA98A9", color = "#AA98A9"),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt") 
            ) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse), scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal) +
      scale_color_manual(values = snv_pal)
    # labs(fill = "Strain")
    # cage_2
    
    cage_3 <-  ggplot(strain_trajectory_df %>% filter(cage == "Cage 3"), aes(x = gut_label, y = freq, fill = strain)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE), show.legend = FALSE) +
      geom_errorbar(data = error_bars %>% filter(cage == "Cage 3"),
                    aes(x = gut_label, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      geom_point(data = strain_trajectory_df %>% filter(cage == "Cage 3", strain == 1), aes(y = allele_frequency, color = locus), size = snv_dot_size, show.legend = FALSE, shape = 8, position = position_jitter(width = 0.1, height = 0)) +
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
            # strip.background = element_rect(fill = scales::alpha("#08519c", 0.85), color = "black"),
            # strip.text = element_text(color = "white", size = font_size), 
            strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(color = "black", size = font_size),
            plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef"),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt") 
            ) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      scale_x_discrete(limits = gut_label_order) +
      facet_wrap(~paste0("Mouse ",mouse), scales = "free_x") +
      # scale_fill_brewer(palette = pal)
      scale_fill_manual(values = pal) +
      scale_color_manual(values = snv_pal)
    # labs(fill = "Strain")
    # cage_3
  } else if (output == "legend") {
    
    # strain_trajectory_df <- strain_trajectory_df %>%
    #   filter(cage == "Cage 3") 
    # 
    unique_snvs <- snv_names
    
    dummy_df <- data.frame(SNV_name = unique_snvs)
    dummy_df['SNV'] <- as.character(1:length(unique_snvs))
    dummy_df['allele_frequency'] <- 1
    dummy_df['dummy_x'] <- "Dummy mouse"
    
    legend_plot <- ggplot(dummy_df) +
      geom_point(data = dummy_df, aes(x = dummy_x, y = allele_frequency, color = SNV), size = snv_dot_size - 1, show.legend = TRUE, shape = 8) +
      theme_bw() +
      theme(text = element_text(family = "Helvetica", size = x_axis_font_size),
            legend.key.width = unit(0.025, "npc"),
            legend.spacing.x = unit(0.025, 'cm'),
            legend.spacing.y = unit(0.025, 'cm'),
            legend.text = element_text(margin = margin(l = -3, unit = "pt")),
            legend.margin = margin(0, 0, 0, 0, "pt"),
            legend.title = element_text(hjust = 0.5)
            ) +
      scale_color_manual(values = snv_pal) +
      guides(color = guide_legend(ncol = if_else(length(unique_snvs) > 5, 2, 1), byrow = FALSE)) 
    
    # strain_trajectory_df <- strain_trajectory_df %>%
    #   rowwise() %>%
    #   mutate(locus_new = paste0("snv ",as.character(which(unique_snvs == locus))))
    # 
    # strain_trajectory_df$locus_new <- factor(strain_trajectory_df$locus_new, levels = strain_trajectory_df %>% pull(locus_new) %>% unique())
    # 
    # 
    # cage_3 <-  ggplot(strain_trajectory_df) +
    #   geom_bar(aes(x = gut_label, y = freq, fill = strain), stat = "identity", position = position_stack(reverse = TRUE), show.legend = FALSE) +
    #   geom_point(data = strain_trajectory_df %>% filter(strain == 1), aes(x = gut_label, y = allele_frequency, color = locus_new), size = snv_dot_size, show.legend = TRUE, shape = 8, position = position_jitter(width = 0.1, height = 0)) +
    #   theme_bw() +
    #   theme(text = element_text(family = "Helvetica", size = x_axis_font_size),
    #         axis.title = element_blank(), 
    #         # legend.position = "none",
    #         legend.key.width = unit(0.025, "npc"),
    #         # axis.text.x = element_text(size = 5, colour = c("#72c8f1","#006c37","#dcc573","#d05c6f","#312c77")),
    #         # axis.text.x = element_text(size = 5, colour = c("darkgreen","darkgreen","darkgreen","goldenrod4","goldenrod4")),
    #         axis.text.x = element_text(size = x_axis_font_size),
    #         legend.spacing.x = unit(0.025, 'cm'),
    #         legend.spacing.y = unit(0.025, 'cm'),
    #         legend.title = element_blank(),
    #         # legend.key.width = unit(0.01, "npc"),
    #         # legend.key.size = unit(0.8, 'lines'),
    #         axis.text.y = element_blank(), 
    #         axis.ticks.y = element_blank(), 
    #         strip.background = element_rect(fill = "white", color = "black"),
    #         strip.text = element_text(color = "black", size = font_size),
    #         plot.background = element_rect(fill = "#e7e1ef", color = "#e7e1ef")) +
    #   scale_y_continuous(breaks = c(0, 0.5, 1)) +
    #   scale_x_discrete(limits = gut_label_order) +
    #   facet_wrap(~paste0("Mouse ",mouse), scales = "free_x") +
    #   # scale_fill_brewer(palette = pal)
    #   scale_fill_manual(values = pal) +
    #   scale_color_manual(values = snv_pal) +
    #   guides(fill = "none")
    # # labs(fill = "Strain")
    # # cage_3
  }
  
  
  if (output == "normal") {
    inoc_size <- 0.075 #0.0675
    height <- 0.88
    
    grid <- ggdraw() +
      draw_plot(inoculum, 0, 0.0, inoc_size, height) +
      draw_plot(cage_1, inoc_size, 0.0, (1-inoc_size)*(3/8), height) +
      draw_plot(cage_2, inoc_size+(1-inoc_size)*(3/8), 0.0, (1-inoc_size)*(2/8), height) +
      draw_plot(cage_3, inoc_size+(1-inoc_size)*(3/8)+(1-inoc_size)*(2/8), 0.0, (1-inoc_size)*(3/8), height) +
      draw_plot_label(c(paste(strsplit(species_name, "_")[[1]][1], strsplit(species_name, "_")[[1]][2])), c(0.045), c(1),hjust = 0,vjust = 1.25, size = font_size, family = "Helvetica", fontface = "italic")
    
    return(grid)
  } else if (output == "legend") {
    legend <- get_legend(legend_plot)
    return(as_ggplot(legend))
  }
  
}

create_label <- function(host_orientation, cage_orientation, diet_orientation, gut_orientation) {
  if (diet_orientation == "Different diet") {
    return("Different cage,\ndifferent diet")
  } else if (cage_orientation == "Different cage") {
    return("Different cage,\nsame diet")
  } else if (host_orientation == "Between host") {
    return("Same cage,\ndifferent host")
  } else {
    return("Same host")
  }
}

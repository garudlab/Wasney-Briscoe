library(dplyr)
library(ggplot2)
library(grid)
library(cowplot)
library(ggtext)
library(stringr)
library(lubridate)

font_size <- 12
x_axis_font_size <- 8

snv_dot_size = 3


plot_snvs <- function(species = "Bacteroides_vulgatus_57955", subject = 8, min_coverage = 20, snv_filter = FALSE, focal_contig = NA, focal_site_pos = NA, sample_to_polarize_by = NA, repolarize = TRUE, output = "normal", pal =  c("#654321", "#D2B48C"), timepoints_to_use = c(), use_simple_timepoint = FALSE, brief_x_labels = TRUE, path = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/strain_phasing/strain_freqs.csv") {
  
  # reassigning variables
  species_oi <- species
  subject_oi <- subject
  strain_pal <- pal
  
  ########################## SNV PROCESSING ##########################
  
  # Loading SNVs
  snv_data_path <- "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/evolutionary_changes/SNV_freqs_WithinTimepoint.tsv"
  
  snv_data <- read.csv2(snv_data_path, sep = "\t") %>% 
    rename(subject = subject_id, 
           timestamp = timepoint,
           sample_type_number = sample_type) %>%
    mutate(allele_frequency = as.numeric(allele_frequency),
           depth = as.numeric(depth)) %>%
    rowwise() %>%
    mutate(allele_frequency = ifelse(depth >= min_coverage, allele_frequency, NA)) %>% 
    select(-timestamp)
  
  if (snv_filter) {
    snv_data <- snv_data %>%
      filter(contig %in% focal_contig, site_pos %in% focal_site_pos)
    
  }
  
  snv_names <- snv_data %>%
    filter(species == species_oi, subject == subject_oi) %>%
    mutate(locus = paste0(contig, ", ", site_pos)) %>%
    select(locus) %>%
    unique() %>%
    pull()
  
  if (length(snv_names) > 12) {
    snv_pal <- rainbow(n = length(snv_names))
  } else { 
    snv_pal <- c('#ffff99','#e31a1c','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#a6cee3','#1f78b4','#33a02c','#cab2d6','#6a3d9a','#b15928')
  }
  
  
  snv_freqs <- snv_data %>%
    filter(species == species_oi,subject == subject_oi) %>%
    mutate(allele_frequency = as.numeric(allele_frequency))
  
  all_snvs <- snv_freqs %>% select(site_id, contig, site_pos) %>% unique()

  ########################## STRAIN PROCESSING ##########################
  
  # Loading strain frequency df
  strain_freq_df <- read.csv2(path, sep = ",") %>%
    mutate(freq = as.numeric(freq), 
           lower_ci = as.numeric(lower_ci),
           upper_ci = as.numeric(upper_ci)) 
  
  # filter strain and SNV dataframes
  strain_freq_df <- strain_freq_df %>%
    filter(species == species_oi, subject == subject_oi) 
  
  
  # Switch strain order such that strain 1 is the "major" strain
  strain_freq_df <- strain_freq_df %>%
    group_by(strain) %>%
    mutate(mean_freq = mean(freq)) %>%  
    ungroup() %>%
    distinct(strain, mean_freq) %>%     
    arrange(desc(mean_freq)) %>%
    mutate(new_strain_label = row_number()) %>%  
    select(strain, new_strain_label) %>%
    right_join(strain_freq_df, by = "strain") %>%
    select(-sample_type, -tissue)
  
  
  ########################## MERGING SNVS WITH STRAINS ##########################
  
  strain_freq_df <- full_join(
    snv_freqs,
    strain_freq_df,
    by = c("sample", "date", "time", "species", "subject", "sample_type_number")) %>%
    mutate(
      strain = if_else(is.na(freq), 1L, strain)
    )
  
  ########################## POST-MERGE ##########################
  
  # Label with timepoint
  strain_freq_df <- strain_freq_df %>%
    ungroup() %>%  
    distinct(date, time) %>%                    
    arrange(date, time) %>%
    mutate(timepoint = paste0("timepoint ", row_number())) %>%
    right_join(strain_freq_df, by = c("date", "time")) %>%
    arrange(date, time)
  
  
  # Subsetting timepoints if necessary
  
  if (length(timepoints_to_use) != 0) {
    
    strain_freq_df <- strain_freq_df %>%
      filter(as.numeric(str_remove(timepoint, "timepoint ")) %in% timepoints_to_use)
    
    strain_freq_df <- strain_freq_df %>%
      select(-timepoint) %>%
      distinct(date, time) %>%                    
      arrange(date, time) %>%
      mutate(timepoint = paste0("timepoint ", row_number())) %>%
      right_join(strain_freq_df %>% select(-timepoint), by = c("date", "time")) %>%
      arrange(date, time)
  }
  
  ## splitting duplicate names
  
  strain_freq_df <- strain_freq_df %>%
    group_by(timepoint, sample_type_number, strain, site_id) %>%
    mutate(dup_index = row_number(),
           no_of_duplicates = n()) %>%
    ungroup() %>%
    mutate(
      sample_type_number = if_else(
        no_of_duplicates == 1,
        as.character(sample_type_number),
        paste0(sample_type_number, " (", dup_index, ")")
      )
    ) %>%
    select(-dup_index, -no_of_duplicates)  # optional: drop helper column
  
  
  ## add capsules that are missing
  
  ### Strains
  timepoint_groups <- split(strain_freq_df, list(strain_freq_df$date, strain_freq_df$time, strain_freq_df$timepoint, strain_freq_df$sample), drop = TRUE)
  full_capsule_list <- c("Capsule 1", 
                         "Capsule 2", 
                         "Capsule 3", 
                         "Capsule 4")
  
  
  for (timepoint_name in names(timepoint_groups)) {
    timepoint_group_subset <- timepoint_groups[[timepoint_name]]
    
    # get sample type list
    sample_type_list <- timepoint_group_subset$sample_type_number %>% unique()
    sample_type_list <- gsub("^Stool \\(\\d+\\)$", "Stool", sample_type_list)
    
    if (("Stool" %in% sample_type_list)) {
      
      next
      
    } else {
      
      # getting date, time and timepoint
      parts <- strsplit(timepoint_name, "\\.")[[1]]
      date_value <- parts[1]
      time_value <- parts[2]
      timepoint_value <- parts[3]
      sample_value <- parts[4]
      
      # Identifying missing capsules
      capsules_present <- timepoint_group_subset$sample_type_number %>% unique()
      missing_capsules <- setdiff(full_capsule_list, capsules_present)
      
      for (capsule in missing_capsules) {
        for (snv_i in 1:nrow(all_snvs)) {
          
          new_capsule_row_1 <- c(date_value, 
                                 time_value, 
                                 timepoint_value, 
                                 species_oi, 
                                 subject_oi, 
                                 capsule, 
                                 all_snvs[snv_i,1][[1]],
                                 all_snvs[snv_i,2][[1]],
                                 all_snvs[snv_i,3][[1]],
                                 substr(all_snvs[snv_i,1][[1]], nchar(all_snvs[snv_i,1][[1]]), nchar(all_snvs[snv_i,1][[1]])), 
                                 sample_value, 
                                 NA,
                                 NA,
                                 1,
                                 1,
                                 NA,
                                 NA,
                                 NA,
                                 NA,
                                 NA)
          new_capsule_row_1 <- as.data.frame(t(new_capsule_row_1), stringsAsFactors = FALSE)
          colnames(new_capsule_row_1) <- colnames(strain_freq_df)
          
          new_capsule_row_2 <- c(date_value, 
                                 time_value, 
                                 timepoint_value, 
                                 species_oi, 
                                 subject_oi, 
                                 capsule, 
                                 all_snvs[snv_i,1][[1]],
                                 all_snvs[snv_i,2][[1]],
                                 all_snvs[snv_i,3][[1]],
                                 substr(all_snvs[snv_i,1][[1]], nchar(all_snvs[snv_i,1][[1]]), nchar(all_snvs[snv_i,1][[1]])), 
                                 sample_value, 
                                 NA,
                                 NA,
                                 2,
                                 2,
                                 NA,
                                 NA,
                                 NA,
                                 NA,
                                 NA)
          new_capsule_row_2 <- as.data.frame(t(new_capsule_row_2), stringsAsFactors = FALSE)
          colnames(new_capsule_row_2) <- colnames(strain_freq_df)
          
          strain_freq_df <- rbind(strain_freq_df, new_capsule_row_1, new_capsule_row_2)
          
        }
      }
    }
  }
  
  # coerce to numeric
  
  strain_freq_df <- strain_freq_df %>%
    mutate(freq = as.numeric(freq), 
           allele_frequency = as.numeric(allele_frequency),
           lower_ci = as.numeric(lower_ci),
           upper_ci = as.numeric(upper_ci)) 
  
  # Shift timezones
  
  strain_freq_df <- strain_freq_df %>%
    rowwise() %>%
    mutate(
      # combine date and time into a single datetime string
      datetime_utc = ymd_hms(paste(date, time), tz = "UTC"),
      
      # convert to Pacific time
      datetime_pst = format(with_tz(datetime_utc, tzone = "America/Los_Angeles"), "%Y-%m-%d %H:%M:%S %Z"),
      
      # split
      date_pst = str_split(datetime_pst, " ")[[1]][1],
      time_pst = str_split(datetime_pst, " ")[[1]][2]
    )
  
  # rename timepoints
  
  strain_freq_df <- strain_freq_df %>%
    rowwise() %>%
    mutate(timepoint_label = paste0("T",str_split(timepoint, " ")[[1]][2]))
  
  day_1 <- as.Date(strain_freq_df %>% filter(timepoint == "timepoint 1") %>% head(1) %>% select(date_pst) %>% pull())
  
  strain_freq_df <- strain_freq_df %>%
    rowwise() %>%
    mutate(day = ifelse(timepoint == "timepoint 1", "Day 1", paste0("Day ", as.character((as.Date(date_pst)-day_1)[[1]] + 1))),
           time_description = ifelse(as.numeric(str_split(time_pst, ":")[[1]][1]) >= 20, "evening", "afternoon"),
           time_description = ifelse(as.numeric(str_split(time_pst, ":")[[1]][1]) <= 12, "morning", time_description),
           timepoint_description = paste0(day, ", ", time_description))
  
  
  # splitting timepoint description duplicate names
  
  
  strain_freq_df <- strain_freq_df %>%
    group_by(timepoint_description) %>%
    mutate(no_of_repeated_time_descriptions = n_distinct(time),
           repeat_index = dense_rank(time)) %>%
    ungroup() %>%
    mutate(
      timepoint_description = if_else(
        no_of_repeated_time_descriptions == 1,
        as.character(timepoint_description),
        paste0(timepoint_description, " (", repeat_index, ")")
      )
    ) %>%
    select(-no_of_repeated_time_descriptions, -repeat_index)  # optional: drop helper column
  
  
  
  # Factor
  strain_order <- strain_freq_df$new_strain_label %>% unique() %>% sort() %>% as.character()
  strain_freq_df$new_strain_label <- factor(strain_freq_df$new_strain_label, levels = strain_order)
  timepoint_order <- strain_freq_df$timepoint_label %>% unique() %>% sort()
  timepoint_order <- timepoint_order[order(as.numeric(substring(timepoint_order, 2)))]
  strain_freq_df$timepoint_label <- factor(strain_freq_df$timepoint_label, levels = timepoint_order)
  timepoint_description_order <- strain_freq_df %>% arrange(date_pst, time_pst) %>% select(timepoint_description) %>% unique() %>% pull()
  strain_freq_df$timepoint_description <- factor(strain_freq_df$timepoint_description, levels = timepoint_description_order)

  
  # title
  
  if (str_split(species_oi, "_")[[1]][2] == "sp") {
    plot_title <- paste0(
      "<i>", str_split(species_oi, "_")[[1]][1], "</i> sp., Subject ", subject_oi, ", healthy human cohort"
    )
  } else {
    plot_title <- paste0(
      "<i>", strsplit(species_oi, "_")[[1]][1], " ",
      strsplit(species_oi, "_")[[1]][2], "</i>, Subject ", subject_oi, ", healthy human cohort"
    )
  }
  
  strain_freq_df <- strain_freq_df %>%
    mutate(xlabel = ifelse(sample_type_number %in% c("Stool", "Saliva"),sample_type_number, str_split(sample_type_number, " ", simplify = TRUE)[, 2]))
  
  # Error bar df
  
  error_bars <- strain_freq_df %>% 
    filter(strain == strain_order[2]) %>%
    rowwise() %>%
    mutate(upper_ci = ifelse(strain == new_strain_label, 1 - upper_ci, upper_ci),
           lower_ci = ifelse(strain == new_strain_label, 1 - lower_ci, lower_ci))
  
  error_bars <- error_bars %>%
    arrange(timepoint, sample_type_number) %>%
    mutate(xlabel = ifelse(sample_type_number %in% c("Stool", "Saliva"),sample_type_number, str_split(sample_type_number, " ", simplify = TRUE)[, 2])) %>%
    select(xlabel, sample_type_number, timepoint_label,timepoint_description, freq, upper_ci, lower_ci, new_strain_label) %>%
    unique()
  
  sample_type_numbers_present <- error_bars %>% select(sample_type_number) %>% unique() %>% pull()
  
  if (("Stool" %in% sample_type_numbers_present) | ("Saliva" %in% sample_type_numbers_present)) {
    brief_x_labels = FALSE
  }

  if (brief_x_labels) {
    strain_freq_df$xlabel <- factor(strain_freq_df$xlabel, levels = c("1", "2", "3", "4", "Stool", "Saliva"))
    error_bars$xlabel <- factor(error_bars$xlabel, levels = c("1", "2", "3", "4", "Stool", "Saliva"))
  }
  
  strain_freq_df_plotting <- strain_freq_df %>% 
    filter(!is.na(freq) | !is.na(allele_frequency)) %>% 
    filter(timepoint_label %in% unique(timepoint_label[!is.na(freq)])) %>%
    filter(site_id == all_snvs[1,1][[1]])
  
  strain_freq_df <- strain_freq_df %>%
    mutate(locus = paste0(contig,"|", site_pos))
  
  # Repolarize
  
  if (repolarize) {
    if (is.na(sample_to_polarize_by)) {
      sample_to_polarize_by <-strain_freq_df %>%
        filter(str_starts(sample_type_number, "Capsule"), !is.na(allele_frequency)) %>%
        arrange(timepoint, sample_type_number) %>%
        head(1) %>%
        select(`sample`) %>%
        pull()
    }
    snvs_to_repolarize <- strain_freq_df %>% 
      filter(sample == sample_to_polarize_by) %>% 
      filter(allele_frequency > 0.5) %>% 
      select(site_id) %>% unique() %>% pull()
    
    strain_freq_df <- strain_freq_df %>%
      rowwise() %>%
      mutate(allele_frequency = ifelse(site_id %in% snvs_to_repolarize, 1 - allele_frequency, allele_frequency))
  }
  
  
  # Plot
  color_names <- c("black", "red", "black", "black","red","black")
  
  if (output == "legend"){ 
    
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
    
    
    } else if (use_simple_timepoint) {
    # Plot
    p <- ggplot(strain_freq_df_plotting, aes(x = sample_type_number, y = freq, fill = new_strain_label)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.9) +
      geom_errorbar(data = error_bars,
                    aes(x = sample_type_number, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      geom_point(data = strain_freq_df %>% filter(strain  == 1), aes(x = sample_type_number, y = allele_frequency, color = locus), size = snv_dot_size, show.legend = FALSE, shape = 8, position = position_jitter(width = 0.1, height = 0))+
      scale_fill_manual(name = "Strain",
                        values = strain_pal) +
      scale_color_manual(values = snv_pal) +
      scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
      theme_bw() + 
      theme(text = element_text(family = "Helvetica"),
            plot.title = element_markdown(size = font_size),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = x_axis_font_size),
            axis.title.x = element_blank(),
            # axis.title.y = element_text(size = font_size),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = x_axis_font_size, colour = color_names),
            strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(color = "black", size = x_axis_font_size),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt"),
            # legend.text = element_text(size = font_size)
            legend.position = "none"
      ) +
      facet_grid(~ timepoint_label, scales = "free", space = "free") +
      labs(title = plot_title, y = "Strain frequency")
  } else if (brief_x_labels) {
    p <- ggplot(strain_freq_df_plotting, aes(x = xlabel, y = freq, fill = new_strain_label)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.9) +
      geom_errorbar(data = error_bars,
                    aes(x = xlabel, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      geom_point(data = strain_freq_df %>% filter(strain  == 1), aes(x = xlabel, y = allele_frequency, color = locus), size = snv_dot_size, show.legend = FALSE, shape = 8, position = position_jitter(width = 0.1, height = 0))+
      scale_fill_manual(name = "Strain",
                        values = strain_pal) +
      scale_color_manual(values = snv_pal) +
      scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
      theme_bw() + 
      theme(text = element_text(family = "Helvetica"),
            plot.title = element_markdown(size = font_size),
            axis.text.x = element_text(hjust = 0.5, vjust = 1, size = x_axis_font_size),
            axis.title.x = element_blank(),
            # axis.title.y = element_text(size = font_size),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = x_axis_font_size, colour = color_names),
            strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(color = "black", size = x_axis_font_size),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt"),
            # legend.text = element_text(size = font_size)
            legend.position = "none"
      ) +
      facet_grid(~ timepoint_description, scales = "free", space = "free", drop = TRUE) +
      labs(title = plot_title, y = "Strain frequency") 
  } else {
    # Plot
    p <- ggplot(strain_freq_df_plotting, aes(x = sample_type_number, y = freq, fill = new_strain_label)) +
      geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.9) +
      geom_errorbar(data = error_bars,
                    aes(x = sample_type_number, ymax = upper_ci, ymin = lower_ci), 
                    width = 0.2, color = "white") +
      geom_point(data = strain_freq_df %>% filter(strain  == 1), aes(x = sample_type_number, y = allele_frequency, color = locus), size = snv_dot_size, show.legend = FALSE, shape = 8, position = position_jitter(width = 0.1, height = 0))+
      scale_fill_manual(name = "Strain",
                        values = strain_pal) +
      scale_color_manual(values = snv_pal) +
      scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
      theme_bw() + 
      theme(text = element_text(family = "Helvetica"),
            plot.title = element_markdown(size = font_size),
            axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1, size = x_axis_font_size),
            axis.title.x = element_blank(),
            # axis.title.y = element_text(size = font_size),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = x_axis_font_size, colour = color_names),
            strip.background = element_rect(fill = "white", color = "black"),
            strip.text = element_text(color = "black", size = x_axis_font_size),
            panel.spacing = unit(0.15, "lines"),
            plot.margin = margin(r = 2, l = 0.05, t = 4, b = 4, unit = "pt"),
            # legend.text = element_text(size = font_size),
            legend.position = "none"
      ) +
      facet_grid(~ timepoint_description, scales = "free", space = "free", drop = TRUE) +
      labs(title = plot_title, y = "Strain frequency") 
    
  }
  
  if (output == "legend") {
    legend <- get_legend(legend_plot)
    return(as_ggplot(legend))
  } else {
    return(p)
  }
  
}
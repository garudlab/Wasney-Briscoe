#library("viridis") 
require(ggplot2)
library(RColorBrewer)



# Functions

extract_species <- function(species_str) {
  species_name <- sub(".*__", "", species_str)
  if (length(str_split(species_name, " ")[[1]]) <= 1) {
    species_name <- paste0(species_name, " sp.")
  }
  return(species_name)
}

extract_family <- function(tax) {
  return(sub(".*f__([^;]+);.*", "\\1", tax))
}

extract_order <- function(tax) {
  return(sub(".*o__([^;]+);.*", "\\1", tax))
}

extract_location <- function(sample_name) {
  sample_type <- str_split(sample_name, "_")[[1]][1]
  if (sample_type == "Col") {
    location <- "Colon"
  } else if (sample_type == "Cec") {
    location <- "Cecum"
  } else if (sample_type == "Ile") {
    location <- "Ileum"
  } else if (sample_type == "Jej") {
    location <- "Jejunum"
  } else if (sample_type == "Duo") {
    location <- "Duodenum"
  } else {
    stop(paste("ValueError: Unknown sample_name prefix:", sample_name))
  }
  return(location)
}

extract_cage <- function(sample_name) {
  cage <- as.numeric(str_split(sample_name, "_")[[1]][2])
  return(cage)
}

extract_mouse <- function(sample_name) {
  cage <- as.numeric(str_split(sample_name, "_")[[1]][2])
  mouse_by_cage <- as.numeric(str_split(sample_name, "_")[[1]][3])
  mouse <- (cage - 1)*3 + mouse_by_cage
  return(mouse)
}


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


# Color palette

family_order_colors <- c(
  'Actinomycetales; Bifidobacteriaceae' = '#fcf4dd',
  'Bacillales; Bacillaceae_G' = '#fce5cd',
  'Bacteroidales; Muribaculaceae' = '#fb6f92',
  'Bacteroidales; Rikenellaceae' = '#ff99aa',
  'Christensenellales; Borkfalkiaceae' = '#c0e0de',
  'Christensenellales; CAG-552' = '#a0d6d3',
  'Christensenellales; UBA3700' = '#7fc9c5',
  'Clostridiales; Clostridiaceae' = '#9cadce',
  'Coriobacteriales; Eggerthellaceae' = '#7ec4cf',
  'Erysipelotrichales; Erysipelatoclostridiaceae' = '#d9b29e',
  'Erysipelotrichales; Erysipelotrichaceae' = "#d3ab9e",
  'Haloplasmatales; Turicibacteraceae' = '#f0d6e0',
  'Lachnospirales; Anaerotignaceae' = '#e7b6c6',
  'Lachnospirales; CAG-274' = '#d9a1b8',
  'Lachnospirales; Lachnospiraceae' = '#d291bc',
  'Lactobacillales; Enterococcaceae' = '#caffbf',
  'Lactobacillales; Lactobacillaceae' = '#e9ff70',
  'Monoglobales_A; UBA1381' = '#f0e5a0',
  'Oscillospirales; Acutalibacteraceae' = '#e8dff5',
  'Oscillospirales; Butyricicoccaceae' = '#d1cfe2',
  'Oscillospirales; Oscillospiraceae' = '#cdb4db',
  'Oscillospirales; Ruminococcaceae' = '#ffc6ff',
  'Peptostreptococcales; Anaerovoracaceae' = '#f0c0a0',
  'Peptostreptococcales; Peptostreptococcaceae' = '#e29980',
  'RF39; CAG-1000' = '#c0d6f0',
  'Staphylococcales; Staphylococcaceae' = '#f5b0b0',
  'TANB77; CAG-508' = '#d0c0f0'
)

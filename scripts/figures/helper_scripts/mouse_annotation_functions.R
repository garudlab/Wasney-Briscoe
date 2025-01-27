require(dplyr)
require(tidyr)

### Fill in the bray-curtis dataframe
extract_mouse <- function(accession, metadata_df) {
  mouse <- metadata_df[metadata_df$accession == accession, "subject_id"]
  return(mouse)
}

extract_mouse_number <- function(accession) {
  mouse <- substr(accession, 2,2)
  if (mouse == "L") {
    mouse = "Inoculum"
  }
  return(mouse)
}

extract_bc <- function(sample_id_1, sample_id_2, bc_dataframe){
  bc_value <- bc_dataframe[sample_id_1, sample_id_2]
  return(bc_value)
}

extract_diet <- function(accession, metadata_df){
  if (accession == "TL1gDNAshort") {
    return("Inoculum")
  } else {
    diet <- metadata_df[metadata_df$accession == accession, "Diet"]
    return(diet)
  }
}

same_or_different_diet <- function(diet_1, diet_2) {
  if(diet_1 == diet_2) {
    if (diet_1 == "Control") {
      return("Same Diet: Control")
    } else {
      return("Same Diet: Guar Gum")
    }
  } else {
    return("Different Diet")
  }
}

extract_diet_comparison <- function(diet_1, diet_2, within_between) {
  if(within_between == "within host") {
    if (diet_1 == "Guargum" | diet_1 == "Guar Gum") {
      return("Guar Gum")
    } else {
      return("Control")
    }
  } else {
    if(diet_1 == "Control" & diet_2 == "Control") {
      return("Control v Control")
    } else if((diet_1 == "Guargum" & diet_2 == "Guargum") | (diet_1 == "Guar Gum" & diet_2 == "Guar Gum")) {
      return("Guar Gum v Guar Gum")
    } else {
      return("Control v Guar Gum")
    }
  }
}

same_or_different_tissue <- function(tissue_1, tissue_2) {
  if (tissue_1 == tissue_2) {
    return("Same tissue")
  } else {
    return("Different tissue")
  }
}

stat_compare_col <- function(same_or_different_diet, tissue_comparison) {
  output <- paste0(same_or_different_diet, ", ", tissue_comparison)
  return(output)
}

extract_mouse_or_between <- function(mouse, orientation) {
  if (orientation == "Between hosts") {
    return("Multiple hosts")
  } else {
    return(mouse)
  }
}

extract_tissue_type <- function(accession, metadata_df = metadata) {
  tissue_type <- metadata_df[metadata_df$accession == accession, "tissue_type"]
  return(tissue_type)
}


extract_coarse_position <- function(accession, metadata_df = metadata) {
  #determine coarse position
  if (!accession %in% metadata_df$accession) {
    stop("Error: incorrect accession.")
  }
  if (accession == "TL1gDNAshort") {
    return("Stool")
  } else if (substr(accession, 3, 3) %in% c("D", "J", "I")) {
    return("Upper gut")
  } else if (substr(accession, 3, 3) == "C") {
    return("Lower gut")
  } else {
    stop("Error: incorrect accession.")
  }
}

extract_coarse_position_comparison <- function(tissue_1, tissue_2, metadata_df = mouse_metadata) {
  coarse_position_1 <- extract_coarse_position(tissue_1, metadata_df = metadata)
  coarse_position_2 <- extract_coarse_position(tissue_2, metadata_df = metadata)
  coarse_position_vec <- c(coarse_position_1, coarse_position_2)
  if (sum(coarse_position_vec=="Upper gut") == 2) {
    return("Upper gut v Upper gut")
  } else if (sum(coarse_position_vec=="Upper gut") == 1) {
    return("Upper gut v Lower gut")
  } else {
    return("Lower gut v Lower gut")
  }
}

extract_metapopulation_conclusion <- function(focal_species, sfs_df) {
  conclusion <- ifelse(SFS_strains[SFS_strains$species == focal_species, "SFS_multiple_strains"], "Multiple strains", "Single strain")[1]
  return(conclusion)
}

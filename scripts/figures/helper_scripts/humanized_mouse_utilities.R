different_host_tissue_mapping <- function(metadata){

  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  sites_short = c("Ce","Co", "D","I","J")
  
  tissue_map = list()
  for(s in 1:length(sites)){
    tissue_map[[sites_short[s]]] = sites[s]
  }
  
  pairs = c()
  samples = metadata$sample_id
  
  for(s in samples){
    #s = samples[1]
    num_char = nchar(s)
    tissue = tissue_map[[substr(s,3,num_char-1)]]
    mouse_num = substr(s,2,2)
    eligible_samples = unlist(metadata %>% filter(subject_id != mouse_num,
                                                  tissue_type != tissue) %>% 
                                select(sample_id))
    for(es in eligible_samples){
      pairs =c( pairs, paste0(s,"_v_",es))
    }
    
    
  }
  return(list(all = pairs))
}

diet_mapping <- function(){
  diet_map = list()
  for(i in c(1:3)){
    diet_map[[paste0("M", i)]] = "C" 
  }
  for(i in c(4:6)){
    diet_map[[paste0("M", i)]] = "G" 
  }
  return(diet_map)
  
  
  
}
tissue_mapping <- function(){
  
  mice = c("M1","M2","M3","M4","M5","M6")
  
  
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  sites_short = c("Ce","Co", "D","I","J")
  diet_map = diet_mapping()
  mapping = list()
  for(s in 1:length(sites)){
    site = sites[s]
    site_short = sites_short[s] 
    for(m in 1:length(mice)){
      
      new_element = paste0(mice[m],site_short, diet_map[[mice[m]]])
      
      if(!(site %in% names(mapping))){
        #print("new key")
        mapping[[site]] =  new_element
      }else{
        mapping[[site]] = c(mapping[[site]] ,new_element )
      }
    }
    
    
    
    
    
    
  }
  return(mapping)
}
intra_tissue_mapping <- function(){
  diet_map = diet_mapping()
  mice = c("M1","M2","M3","M4","M5","M6")
  
  long_combos = combn(mice, m = 2)
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  sites_short = c("Ce","Co", "D","I","J")
  diet_map = diet_mapping()
  mapping = list()
  for(s in 1:length(sites)){
    site = sites[s]
    site_short = sites_short[s] 
    for(l in 1:ncol(long_combos)){
      new_element = paste0(long_combos[1,l],site_short, diet_map[[long_combos[1,l]]],
                           "_v_",
                           long_combos[2,l],site_short, diet_map[[long_combos[2,l]]])
      
      if(!(site %in% names(mapping))){
        #print("new key")
        mapping[[site]] =  new_element
      }else{
        mapping[[site]] = c(mapping[[site]] ,new_element )
      }
    }
     
      
      
    
      
    
  }
  return(mapping)
}

intra_inter_host_comparisons <- function(){
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  treatments = c("Control","Guargum")
  long_combos = combn(sites, m = 2)
  
  inter_host = list()
  for(site in sites){
    inter_host[[site]] = c(paste(site, "Control"), paste(site,"Guargum"))
  }
  
  intra_host = list()
  index = 0
  for(treatment in treatments){
    for(l in 1:ncol(long_combos)){
      index = index + 1
      intra_host[[index]] = c(paste(long_combos[1,l], treatment), 
                              paste(long_combos[2,l], treatment))
      
    }
  }
  my_comparisons = c(inter_host,intra_host)
  return(my_comparisons)
}


intra_inter_host_comparisons2 <- function(){
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  treatments = c("Control","Guar Gum")
  long_combos = combn(sites, m = 2)
  
  inter_host = list()
  for(site in sites){
    inter_host[[site]] = c(paste(site, "Control"), paste(site,"Guargum"))
  }
  
  intra_host = list()
  index = 0
  for(treatment in treatments){
    for(l in 1:ncol(long_combos)){
      index = index + 1
      intra_host[[index]] = c(paste(long_combos[1,l], treatment), 
                              paste(long_combos[2,l], treatment))
      
    }
  }
  my_comparisons = c(inter_host,intra_host)
  return(my_comparisons)
}


inter_tissue_comparisons <- function(){
  inter_tissue  = list()
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  long_combos = combn(sites, m = 2)
  for(l in 1:ncol(long_combos)){
    inter_tissue[[l]] = c(long_combos[1,l],long_combos[2,l])
  }
  return(inter_tissue)
}

tissue_pair_mapping <- function(){
  physio_order =  c( "Duodenum" , "Jejunum"  ,"Ileum", "Cecum","Colon")
  short_to_long = list()
  whole_name_pairs = inter_tissue_comparisons()
  for(key in 1:length(whole_name_pairs)){
    #print(key)
    tissue1 = whole_name_pairs[[key]][1]
    tissue2 = whole_name_pairs[[key]][2]
    if(substr(tissue1,1,1) == "C"){
      short1 = substr(tissue1,1,2)
    }else{
      short1 = substr(tissue1,1,1)
    }
    if(substr(tissue2,1,1) == "C"){
      short2 = substr(tissue2,1,2)
    }else{
      short2 = substr(tissue2,1,1)
    }
    short_pair = paste0(short1,"_v_",short2)
    long_pair =  paste0(tissue1," v ",tissue2)
    if(which(physio_order == tissue1) > which(physio_order == tissue2)){
      long_pair =  paste0(tissue2," v ",tissue1)
    }
    short_to_long[[short_pair]] = long_pair
  }
  return(short_to_long)
  
}

inter_tissue_comparisons_general <- function(sites){
  inter_tissue  = list()
  # sites = c( "Duodenum" , "Jejunum"  ,"Ileum", "Cecum", "Ascending Colon", "Transverse Colon",
  #            "Descending Colon","Sigmoid Colon")
  long_combos = combn(sites, m = 2)
  for(l in 1:ncol(long_combos)){
    inter_tissue[[l]] = c(long_combos[1,l],long_combos[2,l])
  }
  return(inter_tissue)
}

exp_1_strategy_1 <- function(ref_dir){
  
  Experiment1_groups = c()
  Experiment1_Mapping = list()
  Experiment1_RevMapping = list()
  Experiment1_RevMapping[["UpperGI"]] = c()
  Experiment1_RevMapping[["LowerGI"]] = c()
  
  sites = c( "D","J", "Ce","Co","I")
  upper_gi = c( "D","J","I" )
  lower_gi = c("Co","Ce")
  combos = combn(sites, m = 2)
  control_mice = c("M1","M2","M3")
  guargum_mice = c("M4","M5","M6")
  for(ct in control_mice){
    for(gu in guargum_mice){
      for(site in sites){
        group_name_x = paste0(ct,site,"C_v_",gu,site,"G")
        Experiment1_groups = c(Experiment1_groups,group_name_x)
        if(site %in% upper_gi){
          Experiment1_Mapping[[group_name_x]] = "UpperGI"
          Experiment1_RevMapping[["UpperGI"]]  = group_name_x
          
        }else{
          Experiment1_Mapping[[group_name_x]] = "LowerGI"
          Experiment1_RevMapping[["LowerGI"]] = group_name_x
        }
        
      }
      
    }
  }
  saveRDS(Experiment1_Mapping,paste0(ref_dir,"Metadata/Complete_Experiment1_Mapping.rds"))
  saveRDS(Experiment1_RevMapping,paste0(ref_dir,"Metadata/Complete_Experiment1_ReverseMapping.rds"))
  
}

exp_1_strategy_2 <- function(ref_dir){
  Experiment1_groups = c()
  Experiment1_Mapping = list()
  Experiment1_RevMapping = list()
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  treatments = c("Control","Guargum")
  for(treatment in  treatments){
    for(site in sites ){
      
      
      Experiment1_RevMapping[[paste0(site,"_",treatment)]] = c(NA)
    }
  }
  
  
  #sites = c( "D","J", "Ce","Co","I")
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  upper_gi = c("Duodenum","Jejunum","Ileum")
  lower_gi = c("Colon", "Cecum")
  combos = combn(sites, m = 2)
  control_mice = c("M1","M2","M3")
  guargum_mice = c("M4","M5","M6")
  for(ct in control_mice){
    for(gu in guargum_mice){
      for(site in sites){
        if(site %in% c("Cecum","Colon")){
          site_abbrev = substr(site,1,2)
        }else{
          site_abbrev = substr(site,1,1)
        }
        
        
        group_name_x = paste0(ct,site_abbrev,"C_v_",gu,site_abbrev,"G")
        Experiment1_Mapping[[group_name_x]] = site
        Experiment1_RevMapping[[paste0(site,"_",treatment)]]  =  c(Experiment1_RevMapping[[paste0(site,"_",treatment)]] ,
                                                                   group_name_x)
        
      }
    }
  }
  
  saveRDS(Experiment1_Mapping,paste0(ref_dir,"Metadata/Specific_Experiment1_Mapping.rds"))
  saveRDS(Experiment1_RevMapping,paste0(ref_dir,"Metadata/Specific_Experiment1_ReverseMapping.rds"))
  
}
## This strategy avoids replicates of mice in multiple pairings, in case we have some data issue that is problematic
# for one sample it won't infect the whole experiment


exp_1_strategy_3 <- function(ref_dir){
  Experiment1_groups = c()
  Experiment1_Mapping = list()
  Experiment1_RevMapping = list()
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  treatments = c("Control","Guargum")
  for(treatment in  treatments){
    for(site in sites ){
      
      
      Experiment1_RevMapping[[paste0(site,"_",treatment)]] = c(NA)
    }
  }
  
  
  #sites = c( "D","J", "Ce","Co","I")
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  upper_gi = c("Duodenum","Jejunum","Ileum")
  lower_gi = c("Colon", "Cecum")
  combos = combn(sites, m = 2)
  control_mice = c("M1","M2","M3")
  guargum_mice = c("M4","M5","M6")
  for(pair in 1:length(control_mice)){
    ct = control_mice[pair]
    gu = guargum_mice[pair]
    for(site in sites){
      if(site %in% c("Cecum","Colon")){
        site_abbrev = substr(site,1,2)
      }else{
        site_abbrev = substr(site,1,1)
      }
      
      
      group_name_x = paste0(ct,site_abbrev,"C_v_",gu,site_abbrev,"G")
      Experiment1_Mapping[[group_name_x]] = paste0(site)
      Experiment1_RevMapping[[paste0(site)]]  =  c(Experiment1_RevMapping[[paste0(site)]] ,
                                                   group_name_x)
      
    }
  }
  
  saveRDS(Experiment1_Mapping,paste0(ref_dir,"Metadata/Experiment1_Mapping.rds"))
  saveRDS(Experiment1_RevMapping,paste0(ref_dir,"Metadata/Experiment1_ReverseMapping.rds"))
  
}





exp_2_strategy_1 <- function(ref_dir){
  
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  sites_short = c( "Ce","Co","D","I","J")
  upper_gi = c("Duodenum","Jejunum","Ileum")
  lower_gi = c("Colon", "Cecum")
  long_combos = combn(sites, m = 2)
  #long_combos = paste0(long_combos[1,],"_", long_combos[2,])
  short_combos = combn(sites_short, m = 2)
  #short_combos = paste0(short_combos[1,],"_",short_combos[2,])
  
  
  treatments = c("Control","Guargum")
  
  Experiment_groups = c()
  Experiment_Mapping = list()
  Experiment_RevMapping = list()
  for(treatment in  treatments){
    for(cmb in 1:ncol(short_combos )){
      
      Experiment_RevMapping[[paste0(treatment,"_",short_combos[1,cmb], "_v_",short_combos[2,cmb])]] = NA
    }
  }
  
  
  control_mice = c("M1","M2","M3")
  guargum_mice = c("M4","M5","M6")
  all_mice = c(control_mice, guargum_mice)
  for(mice in all_mice){
    for( cmb in 1:ncol(short_combos )){
      if(mice %in% control_mice){
        rev_key = paste0("Control_",short_combos[1,cmb], "_v_",short_combos[2,cmb])
        group_name_x = paste0(mice,short_combos[1,cmb],"C_v_",mice,short_combos[2,cmb],"C")
        
      }else{
        rev_key = paste0("Guargum_",short_combos[1,cmb], "_v_",short_combos[2,cmb])
        group_name_x = paste0(mice,short_combos[1,cmb],"G_v_",mice,short_combos[2,cmb],"G")
        
      }
      Experiment_Mapping[[group_name_x]] = rev_key
      Experiment_RevMapping[[rev_key]] = c(Experiment_RevMapping[[rev_key]],group_name_x)
    }
    
    
    
  }
  
  saveRDS(Experiment_Mapping,paste0(ref_dir,"Metadata/Experiment2_Mapping.rds"))
  saveRDS(Experiment_RevMapping,paste0(ref_dir,"Metadata/Experiment2_ReverseMapping.rds"))
  
}




exp_3_strategy_1 <- function(ref_dir){
  
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  treatments = c("Control","Guargum")
  control_mice = c("M1","M2","M3")
  guargum_mice = c("M4","M5","M6")
  all_mice = c(control_mice, guargum_mice)
  control_combos = combn(control_mice, m = 2)
  guargum_combos = combn(guargum_mice , m = 2)
  
  Experiment1_groups = c()
  Experiment1_Mapping = list()
  Experiment1_RevMapping = list()
  for(treatment in  treatments){
    for(site in sites ){
      
      
      Experiment1_RevMapping[[paste0(site,"_",treatment)]] = c(NA)
    }
  }
  
  #sites = c( "D","J", "Ce","Co","I")
  sites = c("Cecum","Colon", "Duodenum","Ileum","Jejunum")
  
  all_combos = cbind(control_combos,guargum_combos)
  for(cmb in 1:ncol(all_combos)){
    
    for(site in sites){
      if(site %in% c("Cecum","Colon")){
        site_abbrev = substr(site,1,2)
      }else{
        site_abbrev = substr(site,1,1)
      }
      if(cmb < 4){
        group_name_x = paste0(all_combos[1,cmb],site_abbrev,"C_v_",all_combos[2,cmb],site_abbrev,"C")
      }else{
        group_name_x = paste0(all_combos[1,cmb],site_abbrev,"G_v_",all_combos[2,cmb],site_abbrev,"G")
      }
      
      
      Experiment1_Mapping[[group_name_x]] = paste0(site,"_",treatment)
      Experiment1_RevMapping[[paste0(site,"_",treatment)]]  =  c(Experiment1_RevMapping[[paste0(site,"_",treatment)]] ,
                                                                 group_name_x)
      
    }
  }
  
  saveRDS(Experiment1_Mapping,paste0(ref_dir,"Metadata/Experiment3_Mapping.rds"))
  saveRDS(Experiment1_RevMapping,paste0(ref_dir,"Metadata/Experiment3_ReverseMapping.rds"))
  
}




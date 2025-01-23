#library("viridis") 
require(ggplot2)
library(RColorBrewer)



get_colors_order <- function(order) {
  Bacteroidales_colors <- c('Bacteroidaceae' = '#fb6f92',
                            'Porphyromonadaceae' = '#ffb3c6',
                            'Prevotellaceae' = '#fce1e4',
                            'Rikenellaceae' = '#fffbff')
  
  Bifidobacteriales_colors <- c('Bifidobacteriaceae' = '#fcf4dd')
  
  Burkholderiales_colors <- c('Sutterellaceae' = '#c1d3fe',
                              'unidentified' = '#daeaf6')
  
  Clostridiales_colors <- c('Clostridiaceae' = '#9cadce',
                            'Eubacteriaceae' = '#d1cfe2',
                            'Lachnospiraceae' = '#d291bc',
                            'Oscillospiraceae' = '#cdb4db',
                            'Ruminococcaceae' = '#ffc6ff',
                            'unidentified' = '#e8dff5')
  
  Coriobacterineae_colors <- c('Coriobacteriaceae' = '#7ec4cf')
  
  Desulfovibrionales_colors <- c('Desulfovibrionaceae' = '#ebd8d0')
  
  Enterobacteriales_colors <- c('Enterobacteriaceae' = '#52b2cf')
  
  Erysipelotrichales_colors <- c('Erysipelotrichaceae' = '#d3ab9e')
  
  Lactobacillales_colors <- c('Enterococcaceae' = '#caffbf',
                              'Lactobacillaceae' = '#e9ff70',
                              'Leuconostocaceae' = '#f1ffc4',
                              'Streptococcaceae' = '#eaf2d7')
  
  Verrucomicrobiales_colors = c('Verrucomicrobiaceae' = '#ffcaaf')
  
  Unidentified_colors = c('unidentified' = '#dddddd')
  
  if (order == "Bacteroidales") {
    return(Bacteroidales_colors)
  } else if (order == "Bifidobacteriales") {
    return(Bifidobacteriales_colors)
  } else if (order == "Burkholderiales") {
    return(Burkholderiales_colors)
  } else if (order == "Clostridiales") {
    return(Clostridiales_colors)
  } else if (order == "Coriobacterineae") {
    return(Coriobacterineae_colors)
  } else if (order == "Desulfovibrionales") {
    return(Desulfovibrionales_colors)
  } else if (order == "Enterobacteriales") {
    return(Enterobacteriales_colors)
  } else if (order == 'Enterobacteriales') {
    return(Enterobacteriales_colors)
  } else if (order == "Erysipelotrichales") {
    return(Erysipelotrichales_colors)
  } else if (order == "Lactobacillales") {
    return(Lactobacillales_colors)
  } else if (order == "Verrucomicrobiales") {
    return(Verrucomicrobiales_colors)
  } else if (order == "Unidentified") {
    return(Unidentified_colors)
  } else {
    errorCondition("Error: not the correct order input")
  }
}


family_order_colors <- c('Bacteroidales; Bacteroidaceae' = '#fb6f92',
                         'Bacteroidales; Porphyromonadaceae' = '#ffb3c6',
                         'Bacteroidales; Prevotellaceae' = '#fce1e4',
                         'Bacteroidales; Rikenellaceae' = '#fffbff',
                         'Bifidobacteriales; Bifidobacteriaceae' = '#fcf4dd',
                         'Burkholderiales; Sutterellaceae' = '#c1d3fe',
                         'Burkholderiales; unidentified' = '#daeaf6',
                         'Clostridiales; Clostridiaceae' = '#9cadce',
                         'Clostridiales; Eubacteriaceae' = '#d1cfe2',
                         'Clostridiales; Lachnospiraceae' = '#d291bc',
                         'Clostridiales; Oscillospiraceae' = '#cdb4db',
                         'Clostridiales; Ruminococcaceae' = '#ffc6ff',
                         'Clostridiales; unidentified' = '#e8dff5',
                         'Coriobacterineae; Coriobacteriaceae' = '#7ec4cf',
                         'Desulfovibrionales; Desulfovibrionaceae' = '#ebd8d0',
                         'Enterobacteriales; Enterobacteriaceae' = '#52b2cf',
                         'Erysipelotrichales; Erysipelotrichaceae' = '#d3ab9e',
                         'Lactobacillales; Enterococcaceae' = '#caffbf',
                         'Lactobacillales; Lactobacillaceae' = '#e9ff70',
                         'Lactobacillales; Leuconostocaceae' = '#f1ffc4',
                         'Lactobacillales; Streptococcaceae' = '#eaf2d7',
                         'Verrucomicrobiales; Verrucomicrobiaceae' = '#ffcaaf',
                         'Unidentified_; unidentified' = '#dddddd')


















stacked_bar_o <- function(full_data_, tax_level ){

  #full_data_ = to_plot
  #tax_level = "family"
  
  

  if(tax_level == "family"){
    
    
    to_plot_ = full_data_ %>% group_by(subject_id,family, order, tissue_type) %>% 
      summarize(TotalAbundance = sum(abundance,na.rm=TRUE))
    to_plot_stats = to_plot_ %>% group_by(order,family) %>% 
      summarize(MeanAbundance = mean(TotalAbundance,na.rm=TRUE))
  }else{
    to_plot_ = full_data_ %>% group_by(subject_id,order, tissue_type) %>% 
      summarize(TotalAbundance = sum(abundance,na.rm=TRUE))
    to_plot_stats = to_plot_ %>% group_by(order) %>% 
      summarize(MeanAbundance = mean(TotalAbundance,na.rm=TRUE))
  }
  #length(unique(to_plot_stats$family))
  length(unique(to_plot_$family))
  
  to_plot_$subject_id = ifelse(to_plot_$subject_id != 0 ,paste0("Mouse ", to_plot_$subject_id), "Inoculum")
  sortie = to_plot_stats %>% arrange(desc(MeanAbundance))
  to_plot_ = to_plot_ %>% arrange(order, family)
  to_plot_$order[to_plot_$order == ""] = "Unidentified"
  to_plot_$family[to_plot_$family == ""] = "unidentified"
  to_plot_$family_order = paste0(to_plot_$order, "; ",to_plot_$family)
  to_plot_$family_order = factor(to_plot_$family_order,levels = unique(to_plot_$family_order ))
  
  #length(unique(to_plot_$family_order))
  
  #head(to_plot_)
  #test = to_plot_ %>% filter(subject_id == "Mouse 3", tissue_type == "Duodenum") 
  
  
  if(tax_level == "family"){
    mx = ggplot(to_plot_, aes(x = subject_id, fill = family_order, y = 100*TotalAbundance)) 
  }else{
    mx = ggplot(to_plot_, aes(x = subject_id, fill = order, y = 100*TotalAbundance)) 
  }
   mx <- mx +
    geom_bar(stat = "identity", colour = "black") + theme_bw() +
    # theme(axis.text.x = element_text(angle = 45, size = 20, colour = "black", hjust = 1), 
    #       axis.title.y = element_text(size = 20), 
    #       legend.title = element_text(size = 18), 
    #       legend.text = element_text(size = 18, colour = "black"), 
    #       legend = element_blank(),
    #       axis.text.y = element_text(colour = "black", size = 20),
    #       strip.text.x = element_text(size = 21, colour = "black"),
    #       title = element_text(colour = "black", size = 18)
    # ) + 
    scale_y_continuous(expand = c(0,0),limits = c(0,100)) + 
    #labs(x = "", y = "Relative Abundance (%)", fill = "Taxonomic Order") + 
    facet_wrap(~tissue_type, nrow = 1,scale = "free")
   
  
   if(tax_level == "family"){
     # manualcolors<-c('black','forestgreen', 'red2', 'orange', 'cornflowerblue',
     #                        'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue',
     #                        'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
     #                        "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse',
     #                        'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
     #                        "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")
     #                        #"#abc4ff",
     #"#d4afb9",,"#eac9c1","#e8d1c5",,"#efcfe3",
     
     
     manualcolors<-c('#dddddd','#fb6f92','#ffb3c6','#fce1e4',"#fffbff",
                               "#fcf4dd",
                              "#c1d3fe",'#daeaf6',
                      
    
                              "#9cadce", "#d1cfe2","#d291bc","#cdb4db","#ffc6ff", "#e8dff5",
                               "#7ec4cf",
          
                              "#ebd8d0",
                               "#52b2cf",
                                "#d3ab9e",
                               "#caffbf","#e9ff70","#f1ffc4","#eaf2d7",
                               "#ffcaaf")
                               
                               
                               
     # manualcolors<-c(
     # mx<- mx + scale_color_viridis(discrete = TRUE, option = "D")+
     #   scale_fill_viridis(discrete = TRUE)
     mx <-  mx + scale_fill_manual(values = manualcolors, name = "Taxonomic Family")+ guides(fill = guide_legend(title = NULL))
   }else{
     mx <- mx + scale_color_brewer(palette = "Set3") + 
       scale_fill_brewer(palette = "Set3")
   }

  
  return(list(mx,sortie))
  
}

### Annotation functions

def extract_gut_site(sample_name):
    if sample_name[2:4] == "Co":
        return("Colon")
    elif sample_name[2:4] == "Ce":
        return("Cecum")
    elif sample_name[2:3] == "I":
        return("Ileum")
    elif sample_name[2:3] == "J":
        return("Jejunum")
    elif sample_name[2:3] == "D":
        return("Duodenum")
    elif sample_name == "TL1gDNAshort":
        return("Inoculum")
    else:
        raise ValueError('Incorrect sample_name: cannot extract gut site')

def extract_mouse_number(sample_name):
    mouse_number = sample_name[1:2]
    if mouse_number == "L":
        mouse_number = "Inoculum"
    
    return(mouse_number)

def extract_diet(sample_name):
    mouse_num = extract_mouse_number(sample_name)
    if mouse_num in ['1','2','3']:
        return("Control diet")
    elif mouse_num in ['4','5','6','7','8']:
        return("Guar gum diet")
    elif mouse_num == 'Inoculum':
        return("Human diet")
    else:
        raise ValueError('Incorrect sample_name: cannot extract diet')
        
def extract_cage(sample_name):
    mouse_num = extract_mouse_number(sample_name)
    
    if mouse_num in ['1','2','3']:
        return("Cage 1")
    elif mouse_num in ['4','5']:
        return("Cage 2")
    elif mouse_num in ['6','7', '8']:
        return("Cage 3")
    elif mouse_num == 'Inoculum':
        return("Inoculum")
    else:
        raise ValueError('Incorrect sample_name: cannot extract diet')

def extract_region(sample_name):
    if sample_name[2:4] == "Co":
        return "Lower gut"
    elif sample_name[2:4] == "Ce":
        return "Lower gut"
    elif sample_name[2:3] == "I":
        return "Upper gut"
    elif sample_name[2:3] == "J":
        return "Upper gut"
    elif sample_name[2:3] == "D":
        return "Upper gut"
    elif sample_name == "TL1gDNAshort":
        return "Inoculum"
    else:
        raise ValueError("Sample name" + sample_name + " not formatted correctly.")
        
def extract_gene_description(gene, gene_descriptions_dict, centroid_dict):
    if gene in gene_descriptions_dict:
        return gene_descriptions_dict[gene]
    elif gene in centroid_dict:
        if centroid_dict[gene] in gene_descriptions_dict:
            return gene_descriptions_dict[centroid_dict[gene]]
        else:
            return ""
    else:
        return ""
        
        
#Parsing L.B.'s pi outputs



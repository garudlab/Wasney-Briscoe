#################################################################################
# parse_midas_data ##############################################################
#################################################################################

def parse_sample_metadata_map(): 
    import config
    
    sample_metadata_map = {}
    
    # First load mouse metadata
    #file = open(config.scripts_directory+"HMP_ids_order.txt","r")
    file = open(config.analysis_directory+"metadata/shalon_metadata.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        accession_id = items[0].strip()
        subject_id = items[1].strip()
        sample_id = items[2].strip()
        Type = items[3].strip()
        swallow_date_time = items[4].strip()
        collection_date = items[5].strip()
        recover_date_time = items[6].strip()
        recover_date = items[7].strip()
        recover_time = items[8].strip()
        sample_set = items[9].strip()
        sample_type = items[10].strip()
        location = items[11].strip()
        day = items[12].strip()
        
        sample_metadata_map[sample_id] = (subject_id, 
                                          accession_id,
                                          Type, 
                                          swallow_date_time, 
                                          collection_date, 
                                          recover_date_time, 
                                          recover_date, 
                                          recover_time, 
                                          sample_set, 
                                          sample_type, 
                                          location, 
                                          day)
        
    file.close()

    return sample_metadata_map

def parse_subject_sample_map(sample_metadata_map = {}): 
    
    import config
    
    if len(sample_metadata_map)==0:
        # Load it 
        sample_metadata_map = parse_sample_metadata_map() #loading sample_metadata_map (created by 
                                                          #parse_sample_metadata_map() at 1551)

    
    subject_sample_map = {} #a dictionary
    for sample_id in sample_metadata_map:
        subject_id, \
        accession_id,\
        Type, \
        swallow_date_time, \
        collection_date, \
        recover_date_time, \
        recover_date, \
        recover_time, \
        sample_set, \
        sample_type, \
        location, \
        day = sample_metadata_map[sample_id]
    
        if subject_id not in subject_sample_map: #If the subject is not in there already, initialize it as a nested dictionary!
            subject_sample_map[subject_id] = {} 
            
        if sample_id not in subject_sample_map[subject_id]: #if the sample_id is not in the nested subject dictionary, 
                                                            #use the sample_id as a key with a set being the item 
                                                            #associated with it
            subject_sample_map[subject_id][sample_id] = set()
            
        subject_sample_map[subject_id][sample_id].add(accession_id) #append the accession_id to the subject:sample set
 
    return subject_sample_map 

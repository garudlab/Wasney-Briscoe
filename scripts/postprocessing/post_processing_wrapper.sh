#!/bin/bash
#$ -N post_processing_MIDAS_data
#$ -l h_data=8G,h_rt=2:00:00,highp
#$ -o /u/project/ngarud/michaelw/microbiome_evolution/microbiome_evolution_MOUSE/postprocessing_scripts/outputs/
#$ -e /u/project/ngarud/michaelw/microbiome_evolution/microbiome_evolution_MOUSE/postprocessing_scripts/errors/
#$ -t 1-85:1 
#$ -cwd

###########################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s|--species_list / -h|--help]"
   echo "options:"
   echo "-s|--species_list               Path to the list of species to be processed."
   echo
   echo "-m|--microbiome_scripts         Path to the list of species to be processed."
   echo
   echo "-d|--data_directory             Path to the merged midas data."
   echo
   echo "-h|--help                       Display this help message and exit script."
   echo
}

############################################################
# Setup                                                    #
############################################################

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done

#Default arguments
species_list=/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/metadata/species_snps.txt
microbiome_scripts_path=/u/project/ngarud/michaelw/microbiome_evolution/microbiome_evolution_MOUSE/
data_directory=/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/merged_data/

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--species_list)
        species_list="$2"
        shift 
        shift 
        ;;
        -m|--microbiome_scripts)
        microbiome_scripts_path="$2"
        shift 
        shift 
        ;;
        -d|--data_directory)
        data_directory="$2"
        shift 
        shift 
        ;;
    esac
done


############################################################
# execute                                                  #
############################################################

#set species
i=$((SGE_TASK_ID))
species_name=$(sed -n "$i"p $species_list)

echo "Processing" $species_name
echo ""

#Load python 2.7 conda environment
source /u/local/Modules/default/init/modules.sh
module load anaconda3
source activate python27_env

#STEP 0.A: generate core genes! ONLY DO THIS IF ON THE FIRST ITERATION

### CHANGE DIRECOTRY TO MICROBIOME DIRECTORY
# cd $microbiome_scripts_path

# echo "Generating core genes"
# echo ""

# if [ "$i" -eq 1 ]; then
#     ### GENERATING snps/species_snps.txt
#     cd ${data_directory}snps/
#     ls -d */ | sed 's/\/$//' > species_snps.txt
#     # Execute the Python script if 'i' is equal to 1
#     python core_gene_utils.py
# fi

#STEP 0.B: annotate p-values
echo "Annotating p_values"
echo ""

cd ${microbiome_scripts_path}error_filtering_cpp_code
g++ -std=c++11 -O3 *.cpp -o annotate_pvalue


### CHANGING TO MICROBIOME DIRECTORY FOR THE REST OF THE SCRIPT
cd $microbiome_scripts_path

#STEP 1: Calculate marker gene coverage
echo "Calculating marker gene coverage"
echo ""

python calculate_marker_gene_coverage.py ${species_name}

#STEP 2: calculate coverage distribution
echo "Calculating coverage distribution"
echo ""

python calculate_coverage_distribution.py ${species_name}

#STEP 3: Calculate error p-values 
echo "Calculating coverage distribution"
echo ""

python calculate_error_pvalues.py ${species_name}

#STEP 4: Calculate SNP prevalences 
echo "Calculating SNP prevalences"
echo ""

python calculate_snp_prevalences.py ${species_name}

#STEP 5: Calculate within person SFS 
echo "Calculating withing person SFS"
echo ""

python calculate_within_person_sfs.py ${species_name}

#STEP 6: Calculate substitution rates 
echo "Calculating substitution rates"
echo ""

python calculate_substitution_rates.py ${species_name}

#STEP 7: Calculate singletons 
echo "Calculating singletons"
echo ""

python calculate_singletons.py ${species_name}

#STEP 8: Calculate private SNVs 
echo "Calculating private snvs"
echo ""

python calculate_private_snvs.py ${species_name}


#STEP 9: Calculate intersample changes 
echo "Calculating intersample changes"
echo ""

python calculate_intersample_changes.py ${species_name}




















#!/bin/bash
#$ -N core_genes
#$ -l h_data=8G,h_rt=2:00:00,highp
#$ -o /u/project/ngarud/michaelw/microbiome_evolution/microbiome_evolution_MOUSE/postprocessing_scripts/outputs/
#$ -e /u/project/ngarud/michaelw/microbiome_evolution/microbiome_evolution_MOUSE/postprocessing_scripts/errors/
#$ -cwd

###########################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s|--species_list / -h|--help]"
   echo "options:"
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
microbiome_scripts_path=/u/project/ngarud/michaelw/microbiome_evolution/microbiome_evolution_MOUSE/
data_directory=/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/merged_data/

#Arguments that are passed
for arg in "$@"
do
    case $arg in
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

#Load python 2.7 conda environment
source /u/local/Modules/default/init/modules.sh
module load anaconda3
source activate python27_env

#Creating species list if needed
cd ${data_directory}snps/

if [ ! -e "species_snps.txt" ]; then
    # If it doesn't exist, create it using ls and sed
    ls -d */ | sed 's/\/$//' > species_snps.txt
    echo "species_snps.txt created."
else
    echo "species_snps.txt already exists."
fi

#Generating core genes
cd ${microbiome_scripts_path}

echo "Generating core genes"
python core_gene_utils.py

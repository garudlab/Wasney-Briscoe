#!/bin/bash
#$ -N ReadSupport_plot_strain_frequencies
#$ -l h_data=24G,h_rt=4:00:00
#$ -o /u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/strain_phasing/outputs/
#$ -e /u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/strain_phasing/errors/
#$ -t 1-200:1
#$ -cwd

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s/--species_list | -h/--help]"
   echo "options:"
   echo "-s|--species_list     The file specifying which species to process."
   echo 
   echo "-h|--help             Display this help message and exit script."
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
species_file=/u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/metadata/species_snps.txt

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--species_list)
        species_file="$2"
        shift # Remove --species list from processing
        shift # Remove value from processing
        ;;
    esac
done

############################################################
# Main program                                             #
############################################################

source /u/local/Modules/default/init/modules.sh
module load anaconda3
source activate python_env

i=$((SGE_TASK_ID))

species_name=$(sed -n "$i"p $species_file)

python strain_phasing.py --species $species_name

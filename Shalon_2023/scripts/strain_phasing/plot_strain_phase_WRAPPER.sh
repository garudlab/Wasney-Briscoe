#!/bin/bash
#$ -N plot_strain_phase
#$ -l h_data=24G,h_rt=4:00:00,highp
#$ -o /u/project/ngarud/michaelw/Diversity-Along-Gut/HumanizedMouse/scripts/strain_phasing/outputs/
#$ -e /u/project/ngarud/michaelw/Diversity-Along-Gut/HumanizedMouse/scripts/strain_phasing/errors/
#$ -t 1-85:1
#$ -cwd

############################################################
# Purpose                                                  #
############################################################

#This plots the strain phase plots

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s/--species_list | -h/--help]"
   echo "options:"
   echo "-s|--species_list           A path to the list of the species to plot.."
   echo
   echo "-h|--help                   Display this help message and exit script."
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
species_list=/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/metadata/strain_phasing_species.txt

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--species_list)
        species_list="$2"
        shift # Remove --first_acc from processing
        shift # Remove value from processing
        ;;
    esac
done



############################################################
# execute                                                  #
############################################################

#set species
i=$((SGE_TASK_ID))
species_name=$(sed -n "$i"p $species_list)

#Load environment
source /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate python_env

python plot_strain_phase.py --species ${species_name}

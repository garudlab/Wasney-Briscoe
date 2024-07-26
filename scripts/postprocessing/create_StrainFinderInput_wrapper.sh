#!/bin/bash
#$ -N strainfinder_preprocessing
#$ -l h_data=16G,h_rt=24:00:00,highp
#$ -o /u/project/ngarud/michaelw/microbiome_evolution/microbiome_evolution_MOUSE/postprocessing_scripts/outputs/
#$ -e /u/project/ngarud/michaelw/microbiome_evolution/microbiome_evolution_MOUSE/postprocessing_scripts/errors/
#$ -t 1-85:1

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s/--species_list | -o/--outdir | -h/--help]"
   echo "options:"
   echo "-s|--species_list     The file used to select which species to process."
   echo
   echo "-o|--outdir           The strainfinder output directory."
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
species_file=/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/metadata/species.txt
outdir=/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/strain_phasing/input/

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--species_list)
        species_file="$2"
        shift # Remove --species list from processing
        shift # Remove value from processing
        ;;
        -o|--outdir)
        outdir="$2"
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
source activate python27_env


i=$((SGE_TASK_ID))

species_name=$(sed -n "$i"p $species_file)

mkdir -p $outdir

python /u/project/ngarud/michaelw/microbiome_evolution/microbiome_evolution_MOUSE/postprocessing_scripts/create_StrainFinderInput.py --species $species_name --outdir $outdir


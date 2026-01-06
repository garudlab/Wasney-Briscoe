#!/bin/bash
#$ -N mgbc_CM_build_genome_index_REAL
#$ -e /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/midas3_mgbc/errors/
#$ -o /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/midas3_mgbc/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=8G
#$ -pe shared 8
#$ -l time=12:00:00
#4 -js 1000000000000000000000000000000000
#$ -l highp



############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-n|--db_name / -d|--db_dir / -s|--species_profile / -b|--b2index_dir / -h|--help]"
   echo "options:"
   echo "-n|--db_name            Database names"
   echo
   echo "-d|--db_dir             Database directory."
   echo
   echo "-s|--species_profile    Path to species_prevalence.tsv"
   echo
   echo "-b|--b2index_dir    Path to species_prevalence.tsv"
   echo
   echo "-c|--num_cores          Number of threads"
   echo
   echo "-h|--help                Display this help message and exit script."
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
db_name="localdb"
db_dir=/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc
species_profile="/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/merge/species/species_prevalence.tsv"
b2index_dir=/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc/one_bt2_indexes
num_cores=8

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -n|--db_name)
        db_name="$2"
        shift
        shift
        ;;
        -d|--db_dir)
        db_dir="$2"
        shift
        shift
        ;;
        -s|--species_profile)
        species_profile="$2"
        shift
        shift
        ;;
        -b|--b2index_dir)
        b2index_dir="$2"
        shift
        shift
        ;;
        -c|--num_cores)
        num_cores="$2"
        shift
        shift
        ;;
    esac
done

############################################################
# Main program                                             #
############################################################

source /u/local/Modules/default/init/modules.sh
module load miniforge/23.11.0
conda activate midasv3

# # make the sample manifest
# run merge

# midas build_bowtie2db \
#   --midasdb_name $db_name \
#   --midasdb_dir $db_dir \
#   --species_profile $species_profile \
#   --select_by sample_counts \
#   --select_threshold 2 \
#   --bt2_indexes_name repgenomes \
#   --bt2_indexes_dir $b2index_dir \
#   --num_cores $num_cores

midas build_bowtie2db \
  --midasdb_name $db_name \
  --midasdb_dir $db_dir \
  --species_profile $species_profile \
  --select_by sample_counts \
  --select_threshold 2 \
  --bt2_indexes_name repgenomes \
  --bt2_indexes_dir $b2index_dir \
  --num_cores $num_cores


midas build_bowtie2db \
  --midasdb_name $db_name \
  --midasdb_dir $db_dir \
  --species_profile $species_profile \
  --select_by sample_counts \
  --select_threshold 2 \
  --bt2_indexes_name pangenomes \
  --bt2_indexes_dir $b2index_dir \
  --num_cores $num_cores
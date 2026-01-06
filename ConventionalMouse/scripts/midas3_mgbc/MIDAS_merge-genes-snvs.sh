#!/bin/bash
#$ -N CM_merge_genes_snps_step
#$ -e /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/midas3_mgbc/errors/
#$ -o /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/midas3_mgbc/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=6G
#$ -pe shared 6
#$ -l time=48:00:00
#$ -l highp
#$ -t 1-2:1
#$ -js 1000000000000000000000000000000000000


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-s|--sample_manifest / -m|--midasdb / -o|--midas_output / -h|--help]"
   echo "options:"
   echo "-s|--sample_manifest     Path to a sample manifest file list_of_samples.tsv."
   echo
   echo "-m|--midasdb             Path to the MIDASDB."
   echo
   echo "-o|--midas_output        Path to the midas_output directory."
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
sample_manifest=/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/metadata/list_of_samples.tsv
midasdb=/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc
midas_output=/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc
num_cores=6

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -s|--sample_manifest)
        sample_manifest="$2"
        shift
        shift
        ;;
        -m|--midasdb)
        midasdb="$2"
        shift
        shift
        ;;
        -o|--midas_output)
        midas_output="$2"
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

if [ "$SGE_TASK_ID" -eq 1 ]; then
    echo "Processing genes"
    midas merge_genes \
      --samples_list $sample_manifest \
      --midasdb_name localdb \
      --midasdb_dir $midasdb \
      --genome_depth 1 \
      --sample_counts 1 \
      --num_cores $num_cores \
      ${midas_output}/merge

elif [ "$SGE_TASK_ID" -eq 2 ]; then
    echo "Processing snps"
    midas merge_snps \
      --samples_list $sample_manifest \
      --midasdb_name localdb \
      --midasdb_dir $midasdb \
      --snp_type 'any' \
      --snp_maf 0.0 \
      --genome_depth 5 \
      --site_depth 3 \
      --sample_counts 1 \
      --num_cores $num_cores \
      ${midas_output}/merge
      
else
    echo "SGE_TASK_ID must be 1 or 2"
    exit 1
fi


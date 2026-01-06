#!/bin/bash
#$ -N extract_snv_frequencies
#$ -e /u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/job_scripts/errors/
#$ -o /u/project/ngarud/michaelw/Diversity-Along-Gut/Shalon_2023/scripts/job_scripts/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=45G
#$ -l time=168:00:00
#$ -l highp


############################################################
# Main program                                             #
############################################################

#Load python 2.7 conda environment
source /u/local/Modules/default/init/modules.sh
module load anaconda3
source activate python27_env

python extract_snv_freqs.py --within_host_changes
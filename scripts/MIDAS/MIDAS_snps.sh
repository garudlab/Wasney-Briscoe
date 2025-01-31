#!/bin/bash
#$ -N MIDAS_snps
#$ -e ~/Wasney-Briscoe/scripts/MIDAS/errors/
#$ -o ~/Wasney-Briscoe/scripts/MIDAS/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=28G
#$ -l time=168:00:00
#$ -l highp
#$ -tc 60 # Throttle to max 60 tasks at a time
#$ -t 1-2:1

#Load modules
source /u/local/Modules/default/init/modules.sh
module load anaconda3
module load singularity

#Load python 2.7 conda environment
conda activate python27_env

export PYTHONPATH=$PYTHONPATH:~/MIDAS_mod
export PATH=$PATH:~/MIDAS_mod/scripts
export MIDAS_DB=~/midas_db_v1.2

#Accessions
readarray accs < ~/Wasney-Briscoe/scripts/accessions.txt
accs=(null ${accs[@]}) # zero to one start index
acc=${accs[$SGE_TASK_ID]}
echo $acc

#Accession paths
readarray accs_paths < ~/Wasney-Briscoe/scripts/accession_paths.txt
accs_paths=(null ${accs_paths[@]}) # zero to one start index
acc_path=${accs_paths[$SGE_TASK_ID]}
echo $acc_path

# readarray accs_paths < /u/project/ngarud/Garud_lab/HumanizedMouse_Batch2/metadata/accession_paths.txt
# accs_paths=(null ${accs_paths[@]}) # zero to one start index
# acc_path=${accs_paths[$SGE_TASK_ID]}
# echo $acc_path

OUTDIR=~/midas_output/${acc}

fastq1=~/fastqs/${acc_path}_1.fq.gz
fastq2=~/fastqs/${acc_path}_2.fq.gz
species_union=~/species_union/${acc}/species_union.txt

### MIDAS SNPs
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif run_midas.py snps $OUTDIR -1 $fastq1 -2 $fastq2 --extra_species_file $species_union

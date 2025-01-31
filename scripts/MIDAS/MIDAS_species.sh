#!/bin/bash
#$ -N MIDAS_species_step
#$ -e ~/Wasney-Briscoe/scripts/MIDAS/errors/
#$ -o ~/Wasney-Briscoe/scripts/MIDAS/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=26G
#$ -l time=48:00:00
#$ -l highp
#$ -t 1-41:1

source /u/local/Modules/default/init/modules.sh
module load anaconda3
source activate python27_env

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

OUTDIR=~/midas_output/${acc}

mkdir -p $OUTDIR

fastq1=~/fastqs/${acc_path}_1.fq.gz
fastq2=~/fastqs/${acc_path}_2.fq.gz

#run_midas.py species $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp
module load singularity
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif run_midas.py species $OUTDIR -1 $fastq1 -2 $fastq2

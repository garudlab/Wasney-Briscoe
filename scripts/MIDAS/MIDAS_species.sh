#!/bin/bash
#$ -N MIDAS_species_step
#$ -e ~/Wasney-Briscoe-2024/scripts/MIDAS/errors/
#$ -o ~/Wasney-Briscoe-2024/scripts/MIDAS/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=26G
#$ -l time=48:00:00
#$ -l highp
#$ -t 1-11:1

source /u/local/Modules/default/init/modules.sh
module load anaconda3
source activate python27_env

export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2


#Accessions
readarray accs < /u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/metadata/accessions.txt
accs=(null ${accs[@]}) # zero to one start index
acc=${accs[$SGE_TASK_ID]}
echo $acc

#Accession paths
readarray accs_paths < /u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/metadata/accession_paths.txt
accs_paths=(null ${accs_paths[@]}) # zero to one start index
acc_path=${accs_paths[$SGE_TASK_ID]}
echo $acc_path

OUTDIR=/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/midas_output/${acc}

mkdir $OUTDIR

fastq1=/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/fastqs/${acc_path}_1.fq.gz
fastq2=/u/project/ngarud/Garud_lab/HumanizedMouse/HumanizedMouse_Batch2/fastqs/${acc_path}_2.fq.gz

#run_midas.py species $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp
module load singularity
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif run_midas.py species $OUTDIR -1 $fastq1 -2 $fastq2

#!/bin/bash
#$ -N MIDAS_species_step
#$ -e /u/project/ngarud/michaelw/PaulAllen/Shalon_2023/scripts/errors/
#$ -o /u/project/ngarud/michaelw/PaulAllen/Shalon_2023/scripts/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=20G
#$ -l time=48:00:00
#$ -l highp
#$ -tc 100 # Throttle to max 100 tasks at a time
#$ -t 1-302:1

source /u/local/Modules/default/init/modules.sh
module load anaconda3
source activate python27_env

export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2


#Accessions
readarray accs < /u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/metadata/SRA_list.txt
accs=(null ${accs[@]}) # zero to one start index
acc=${accs[$SGE_TASK_ID]}
echo $acc

OUTDIR=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/midas_output/${acc}


mkdir $OUTDIR

fastq1=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/fastqs/${acc}/${acc}_1.fastq.bz2
fastq2=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shalon_2023/fastqs/${acc}/${acc}_2.fastq.bz2

#run_midas.py species $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp
module load singularity
singularity exec $H2_CONTAINER_LOC/MIDAS-mod.sif run_midas.py species $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp

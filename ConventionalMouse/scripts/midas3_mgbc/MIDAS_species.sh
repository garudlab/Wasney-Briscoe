#!/bin/bash
#$ -N CM-MGBC_species_step
#$ -e /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/errors/
#$ -o /u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/scripts/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=4G
#$ -pe shared 4
#$ -l time=23:00:00
#$ -l highp
#4 -js 10000000
#$ -t 1-45:1



############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-a|--accession_list / -f|--fastqs / -m|--midasdb / -o|--midas_output / -c|--num_cores / -h|--help]"
   echo "options:"
   echo "-a|--accession_list      Path to a list of sample/accession names."
   echo
   echo "-f|--fastqs              Path to the fastq directory."
   echo
   echo "-m|--midasdb             Path to the MIDASDB."
   echo
   echo "-n--db_name              Name of the MIDASDB to use."
   echo
   echo "-o|--midas_output        Path to the midas_output directory."
   echo
   echo "-c|--num_cores           Number of cores to use."
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
accession_list=/u/project/ngarud/michaelw/Diversity-Along-Gut/ConventionalMouse/metadata/sample_list.txt
fastqs=/u/scratch/m/michaelw/ConventionalMouse_fastqs/fastqs_merged
midasdb=/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc
db_name=localdb
midas_output=/u/project/ngarud/Garud_lab/ConventionalMouse/midas3_mgbc
num_cores=4

#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -a|--accession_list)
        accession_list="$2"
        shift 
        shift 
        ;;
        -f|--fastqs)
        fastqs="$2"
        shift 
        shift 
        ;;
        -m|--midasdb)
        midasdb="$2"
        shift 
        shift 
        ;;
        -n|--db_name)
        db_name="$2"
        shift 
        shift 
        ;;
        -o|--midas_output)
        midas_output="$2"
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

#Accessions
readarray accs < $accession_list
accs=(null ${accs[@]}) # zero to one start index
acc=${accs[$SGE_TASK_ID]}
echo $acc

fastq1=${fastqs}${acc}_1.fastq.gz
fastq2=${fastqs}${acc}_2.fastq.gz

midas run_species \
    --sample_name $acc \
    -1 $fastq1 \
    -2 $fastq2 \
    --midasdb_name $db_name \
    --midasdb_dir $midasdb \
    --num_cores $num_cores \
    $midas_output
    
    
# rm -r /u/scratch/m/michaelw/midas_output/${acc}/species/temp/


#!/bin/bash
#$ -N MGBC_MIDAS_snvs_step
#$ -e /u/project/ngarud/michaelw/MIDASv3-postprocessing/MIDASv3-mgbc/errors/
#$ -o /u/project/ngarud/michaelw/MIDASv3-postprocessing/MIDASv3-mgbc/outputs/
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=5G
#$ -l time=72:00:00
#$ -l highp
#$ -pe shared 6
#$ -t 1-42:1


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: [-a|--accession_list / -f|--fastqs / -m|--midasdb / -o|--midas_output / -b|--b2index_dir / -h|--help]"
   echo "options:"
   echo "-a|--accession_list      Path to a list of sample/accession names."
   echo
   echo "-f|--fastqs              Path to the fastq directory."
   echo
   echo "-m|--midasdb             Path to the MIDASDB."
   echo
   echo "-o|--midas_output        Path to the midas_output directory."
   echo
   echo "-b|--b2index_dir         Path to genome index directory"
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
accession_list=/u/project/ngarud/michaelw/hsiao_data/metadata/MIDASv3/accessions_noblanks.txt
fastqs=/u/project/ngarud/Garud_lab/hsiao_data/fastqs_renamed/
midasdb=/u/scratch/m/michaelw/assembling_mouse_db/MIDASDB-mgbc
midas_output=/u/project/ngarud/Garud_lab/hsiao_data/MIDASv3-mgbc
b2index_dir=/u/scratch/m/michaelw/assembling_mouse_db/one_bt2_indexes-mgbc


#Arguments that are passed
for arg in "$@"
do
    case $arg in
        -a|--accession_list)
        accession_list="$2"
        shift # Remove --merged_path from processing
        shift # Remove value from processing
        ;;
        -f|--fastqs)
        fastqs="$2"
        shift # Remove --merged_path from processing
        shift # Remove value from processing
        ;;
        -m|--midasdb)
        midasdb="$2"
        shift # Remove --merged_path from processing
        shift # Remove value from processing
        ;;
        -o|--midas_output)
        midas_output="$2"
        shift # Remove --merged_path from processing
        shift # Remove value from processing
        ;;
        -b|--b2index_dir)
        b2index_dir="$2"
        shift # Remove --merged_path from processing
        shift # Remove value from processing
        ;;
    esac
done

############################################################
# Main program                                             #
############################################################

source /u/local/Modules/default/init/modules.sh
module load miniforge/23.11.0
conda activate midasv3

# read species

readarray species_list < $accession_list
species_list=(null ${species_list[@]}) # zero to one start index
species_id=${species_list[$SGE_TASK_ID]}
echo $species

#Accessions
readarray accs < $accession_list
accs=(null ${accs[@]}) # zero to one start index
acc=${accs[$SGE_TASK_ID]}
echo $acc

fastq1=${fastqs}${acc}_1.fastq.gz
fastq2=${fastqs}${acc}_2.fastq.gz
    
# run SNVs
echo "Processing snvs"
midas run_snps \
    --sample_name $acc \
    -1 $fastq1 \
    -2 $fastq2 \
    --midasdb_name localdb \
    --midasdb_dir $midasdb \
    --select_by median_marker_coverage,unique_fraction_covered \
    --select_threshold=2,0.5 \
    --num_cores 6 \
    $midas_output

   #  midas run_snps \
   #  --sample_name $acc \
   #  -1 $fastq1 \
   #  -2 $fastq2 \
   #  --midasdb_name localdb \
   #  --midasdb_dir $midasdb \
   #  --prebuilt_bowtie2_indexes ${b2index_dir}/repgenomes \
   #  --prebuilt_bowtie2_species ${b2index_dir}/repgenomes.species \
   #  --select_threshold=-1 \
   #  --num_cores 6 \
   #  $midas_output
    




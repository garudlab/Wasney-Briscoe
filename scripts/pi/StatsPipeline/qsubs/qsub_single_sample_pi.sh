#!/bin/bash

AccList1=$3
study=$4
first_count_input=$5
calctype=$6
statistic=$7
covmin=$8
summarymode=$9
strainCOUNTER=0
COUNTER="$(($first_count_input-1))"

echo $strainCOUNTER
echo $COUNTER

while IFS= read -r strain;
do
    echo $strain
    strainCOUNTER=$((strainCOUNTER + 1)); 


    while IFS= read -r line1; do
        echo "File 1: $line1"

        accession1=$line1
  
        COUNTER=$((COUNTER + 1)); 

        if [ $2 == "1000" ]
        then 
            echo "1000!"
            lengthfile=1000
        else
            lengthfile=$(head -n $strainCOUNTER $2 | tail -1)
        fi
        #echo $strain
        #echo $lengthfile
        echo $COUNTER

        echo $accession1 $lengthfile $strain $study $calctype $statistic $covmin > inputs/data_$COUNTER.in; 

    done < $AccList1


done < $1

# qsub -cwd -V -N "$statistic"_"$study"_cov"$covmin" -e errors/ -o errors/ -l h_data=10G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_single_sample_pi.sh"

# # With priority
# qsub -cwd -V -N "$statistic"_"$study"_cov"$covmin" -e errors/ -o errors/ -js 10000 -l h_data=10G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_single_sample_pi.sh"

# Higher mem and time
qsub -cwd -V -N "$statistic"_"$study"_cov"$covmin" -e errors/ -o errors/ -l h_data=24G,time=72:00:00 -b y -t $first_count_input:$COUNTER "./run_single_sample_pi.sh"




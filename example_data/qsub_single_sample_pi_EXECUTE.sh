#!/bin/bash

workpath=~/Wasney-Briscoe/
study=SiteDownsampled
calctype=Schloissnig
stat=Pi
covmin=4

${workpath}scripts/pi/StatsPipeline/qsubs/qsub_single_sample_pi.sh ${workpath}metadata/Bacteroides_vulgatus_57955_list.txt ${workpath}metadata/Bacteroides_vulgatus_57955_lengths.txt ${workpath}metadata/accessions.txt $study 1 $calctype $stat $covmin 0 

#!/bin/bash

workpath=~/Wasney-Briscoe/
study=SiteDownsampled
calctype=Schloissnig
stat=Pi
covmin=4

./qsub_single_sample_pi.sh ${workpath}metadata/species_snps.txt ${workpath}metadata/species_snps_lengths.txt ${workpath}metadata/accessions.txt $study 1 $calctype $stat $covmin 0 

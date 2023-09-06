#!/bin/bash

##Lancement du script
# qsub script_analyseDESMAN.sh 

## utilisateur à prévenir

#$ -M ouleye.sidibe@inrae.fr


# Export de toutes les variables d'environnement
#$ -V


# Sortie d'erreurs
#$ -e /home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/DESMAN_iter/analyseDESMAN.err 


conda activate DESMAN_python

directory="/home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/DESMAN_iterAddPost"
gene_name="/home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/genes_list.txt"
freq_dir="/home/osidibe/work/fluGenAvi_fugace/freq_allelic_OSI"
scripts="/home/osidibe/work/fluGenAvi_fugace/scripts"


## nombre de variant analysé avec DESMAN
nb_var=5  
nb_replicate=10
threshold_replicate=8


#Script python pour concatener les haplotypes DES des différents réplicats pour chacun des variant
python3 $scripts/concate_haplo.py -i $directory -g $gene_name -v $nb_var
python3 $scripts/correct_haplo.py -i $directory -f $freq_dir -g $gene_name -r $nb_replicate -s $threshold_replicate


conda deactivate
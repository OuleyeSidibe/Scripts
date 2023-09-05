#!/usr/bin/bash

# Shell à utiliser pour l'execution du job
#$ -S /bin/bash

#Execution du script
#qsub script_mapping_filter.sh <minimum_depth> <coverage_depth> <coverage>

#Utilisateur à avertir
#$ -M ouleye.sidibe@inrae.fr

#Fichier d'erreur
#$ -e /home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/coverage_by_SAMPLE/cover.err

#Fichier de sortie
#$ -o /home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/coverage_by_SAMPLE/cover.out



conda activate test




#Nombre de lectures minimum par base
minimum_depth=$1
depth=$(($minimum_depth-1))

#Pourcentage de couverture de la taille du géne mappée sur laquelle on souhaite calculer la profondeure et couverture total du géne
coverage_Depth=$2
coverage=$3


#Répertoires
script='/home/osidibe/work/fluGenAvi_fugace/scripts'
table_de_comptage='/home/osidibe/work/fluGenAvi_fugace/freq_alleliques_ALA'
directory=`echo '/home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/coverage_by_SAMPLE/out_statbyGENE_'$minimum_depth\_$coverage_Depth\_$coverage`
inputdir='/home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/coverage_by_SAMPLE'


    
# mkdir $inputdir/out_coveragebyGENE

#Générer les tables de couvertures de chaque gène dans les différents échantillons
python3 $script/coverage_byGENE.py -i $table_de_comptage -o $inputdir/out_coveragebyGENE -n $inputdir/gene_size.txt
#Table de prévalence des génes dans les différents échantillons
python3 $script/nbGENES_bySAMPLES.py
#estimation des paramétres de position (moyenne et ecart-type) et couvertures des génes (profondeure, couverture par base et pourcentage de couverture)
mkdir -p $directory/stat_byGENE
python3 $script/stat_byGENE.py -d $depth -c $coverage_Depth -o $directory/stat_byGENE
# generer deux tables statistiques regroupant : nom de géne, nom d'échantillon, moyenne, écart-type et couverture avnt et aprés filtre sur la couverture
python3 $script/table_stat_forGENES.py -c $coverage -i $directory/stat_byGENE -o $directory


conda deactivate
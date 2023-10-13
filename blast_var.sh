#!/bin/bash

## shell à utiliser pour l'execution
#$ -S /bin/bash

#Lancement script
# qsub "script.sh"

## Envoie par mail en debut et fin du job
#$ -M osidibe@inrae.fr

##Nom du job
#$ -N blast_variant

## Nom de la queue
#$ -q short.q

## nombre de CPU 
#$ -pe thread 4 

## export des variables d'env du front au noeud
#$ -V

## Lancer la commande depuis le repertoire où est lancé le script
#$ -cwd


## stdout
#$ -o /home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/Blast_variant/out_blast.txt

##stderr
#$ -e /home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/Blast_variant/err_blast.txt

Directory="/home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/Blast_variant"
databank="/home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/DB"
query="/home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/DB/query"


conda activate blast-2.12.0

blastn -db $databank/resfinder2022Julnohybrids.fsa -query $query/dfrA1.fasta
# blastn -db $databank/resfinder2022Julnohybrids.fsa -query $query/blaTEM/query_2.fasta
# blastn -db $databank/resfinder2022Julnohybrids.fsa -query $query/blaTEM/query_3.fasta

conda deactivate



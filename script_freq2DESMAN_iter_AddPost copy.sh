#! /usr/bin/bash
# Les commentraires qui commencent par '#$' sont
# interprétés par SGE comme des options en ligne

# Execution via qsub
# qsub -q long.q script_freq2DESMAN.sh "<nom du gene>"
# Liste des gènes pour trouver le nom : /work_projet/ala/metachick-fugace/Analyses_ALA/Snakefile_resfinder_55ech_fev2022/gene_names.txt

# Shell à utiliser pour l'execution du job
#$ -S /bin/bash

# Utilisateur à avertir
#$ -M ouleye.sidibe@inrae.fr

# Export de toutes les variables d'environnement
#$ -V

# Avertir au début (b)egin, à la fin (e)nd, à l'éliminaton (a)bort et
# à la suspension (s)uspend d'un job
#$ -m bea

# Vous pouvez utiliser '-j y' pour ajouter stderr avec stdout
#$ -o /home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/DESMAN_iterAddPost/freq_to_iter_desman.out
# Sortie d'erreur (ne pas utiliser cette option avec '-j y')
#$ -e /home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/DESMAN_iterAddPost/freq_to_iter_desman.err

# Lance la commande depuis le répertoire où est lancé le script
#$ -cwd

source /usr/local/genome/Anaconda3/
conda activate DESMAN_python
source ~/work/venv/venvDESMAN_python/bin/activate

export PATH=/home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/desman_fugace/scripts:/home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/desman_fugace/bin:$PATH
export PYTHONPATH="${PYTHONPATH}:/home/osidibe/.local/lib/python3.8/site-packages/"
gene=$1


script='/save_projet/ala/metachick-fugace/flugenavi/snakefile/'
scriptOSI='/home/osidibe/work/fluGenAvi_fugace/scripts'
DESMAN='/home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2'
inputdir='/home/osidibe/work/fluGenAvi_fugace/freq_allelic_OSI'

directory=`echo '/home/osidibe/work/fluGenAvi_fugace/Analyse_postDESMAN_P2/DESMAN_iterAddPost/DESMAN_'$gene`



nb_repet=1000 ## 1000  ## number of Gibbs sampling steps
nb_g=5    ## 10  ## number of variants
nb_repid=10   ## 5  ## number of replicates


mkdir -p $directory
python3 $script/freq_to_desman_AL.py -i $inputdir -o $directory/freq2Desman.txt -g "$gene"
python3 $DESMAN/desman_fugace/desman/Variant_Filter.py  $directory/freq2Desman.txt -o $directory/file.out -p
python3 $script/add_variant_positions.py  -i $directory/freq2Desman.txt -s $directory/file.outsel_var.csv -o $directory/file.outsel_var_corrected.csv -n 200

for g in `seq 1 $nb_g`
do  for repid in `seq 1 $nb_repid`
    do  echo "La variable vaut $g et $repid"
#         desman $directory/file.outsel_var.csv -g $g -e $directory/file.ouir -o $directory/freq2Desman.txt -g "$gene"_df.csv -o $directory/desman.$g.$repid -i $nb_repet -s $repid 
    python3 $DESMAN/desman_fugace/desman.py $directory/file.outsel_var_corrected.csv -g $g -e $directory/file.outtran_df.csv -o $directory/desman.$g.$repid -i $nb_repet -s $repid 
    #Modifier les fichiers Filterd_tau_star
    python3 $scriptOSI/Filtered_Tau.py -i $directory
    python3 $script/convertit_DESMAN_filtredTau_haplotypes.py -i $directory/desman.$g.$repid/Filtered_Tau_star.csv  -o $directory/desman.$g.$repid/Filtered_Tau_star_haplotypes
    python3 $script/find_haplotype_variant.py -p $directory/desman.$g.$repid/Filtered_Tau_star_haplotypes -f $directory/freq2Desman.txt -o $directory/desman.$g.$repid/haplotype.$g.$repid -n $g

         done
done

 cat $directory/*/fit.txt > $directory/all_fit.txt

conda deactivate 
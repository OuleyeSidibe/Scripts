#Ce script est à lancer de façon séquentiel sous bash pour effectuer 
#deux filtres sur l'ensemble des 55 métagénomes mappées ou non sur 771 clusters de GRAs de la base de données resfinder. 

#Trois seuils : "paramétres à définir"
"couverture"
"Profondeur de couverture moyen"
"Redondance dans les 55 échantillons"


#Charger les 55 fichiers.bam dans le repertoire de travail
"path =  /home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/bamfile_alphanum"



#Lancer samtoools coverage sur l'ensemble des 55 fichiers.bam  
for i in ./*.bam; do samtools coverage $i ; done > /work_home/osidibe/fluGenAvi_fugace/analyse_preDESMAN/pre_filter/coverage_by_SAMPLE/coverage55.csv   #Fournit une table avec les résultats samtools des 55 échantillons.



#Premier filre : seuil de couverture "">= 60" et profondeur moyenne de couverture ">=3"
awk -F "\t" '{if ( $6>= 30 && $7 >= 3){print$0}}' coverage55.csv | cut -d "|" -f 2  > coverage55_filtered.csv #table avec la liste des gènes qui présentent une couverture d'au moins 60% ainsi qu'une profondeur 
#de couverture moyenne d'au moins 3 reads par position dans les 55 échantillons.



#Liste non redondante de l'ensemble des gènes qui dépassent le seuil de couverture et de profondeur.
 awk -F "\t" '{if ($0 !~ "#"){print $0}}' coverage55_filtered.csv | awk -F "\t" '!vu[$1]++'| sort -uk1 > geneslist_filtered.csv #table de 65 génes 



#Estimer la redondance des 65 génes dans les 55 echnatillons
awk '{x[$1]++} ; END {for (n in x) print x[n] "\t" n}' coverage55_filtered.csv| sort -uk2  > filtered_list.txt  


#Générer une table finale qui regroupe les statistiques de couverture des gènes dépassant le filtre et leur redondance dans les 55 échantillons
head -n+1 coverage55.csv > genefilter_final.csv | paste geneslist_filtered.csv filtered_list.txt |awk "NF{NF--};1" >> genefilter_final.csv 
#Sur un total de 65 gènes, 38 sont présents dans au moins 10 échantillons
 
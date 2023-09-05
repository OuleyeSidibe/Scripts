                                                                                     
#                                                                                                                                                                                     #Libraries

import pandas as pd
import os                                                      
import sys
import fnmatch
from argparse import ArgumentParser
import glob


######################################################################################################################
##Infos

"""script qui recupére les haplotypes des différents réplicats pour chaque variant du nombre total de variant demandé"""


##lancement script

# python3 concate_haplo.py -i  directory_path -g list_of_gene_name -v variant_number

#######################################################################################################################


#Analyseur d'argument

def config_parameters():
    parser= ArgumentParser()
    parser.add_argument("-i", "--input", dest="directory", default="", help="Desman repository")
    parser.add_argument("-g", "--genes", dest="gene_name", default="", help="list of gene names")
    parser.add_argument("-v", "--variant_nb", dest="nb_var", default="", help="Nombre de variants lancé dans Desman")
    options=parser.parse_args()
    if len(sys.argv)<3:
            sys.exit("warning: wrong number of arguments!!")
    return options.directory, options.gene_name, options.nb_var



def main():

    directory, gene_name, nb_var = config_parameters()
    ##Lecture de la liste des génes
    gene_name =pd.read_csv(gene_name, index_col=[0], header=None)

    ##Pour chaque géne de la liste des génes analysés par DESMAn
    for gene in gene_name.index:
        # print(gene)
        
        ##Extraire les positions variables détectés
        Pos_var = pd.read_csv(f"{directory}/DESMAN_{gene}/file.outsel_var.csv", usecols=[1])
        # stocker les positions variables dans une liste
        list_PosVar = [x for x in Pos_var.values]
        
        ##Pour chaque nombre de variant
        for i in range(1, int(nb_var)+1):

            ##Créer un df avec comme donnée les positions variables du géne
            result = pd.DataFrame(Pos_var)
            # print(result)

            ##Lister les sous-repertoires desman du repertoire parent
            dir_gene = os.path.join(directory + "/DESMAN_" + gene)
            for desman_file in os.listdir(dir_gene):

                ##Si le SR a le pattern "desman.*.*" :
                if fnmatch.fnmatch(desman_file, f"desman.{i}.*"):
                    # print(desman_file)
                    
                    ##Créer un df vide
                    haplo_file = pd.DataFrame()
                    ##Pour chaque position variable d
                    for x in list_PosVar:
                        ##Lire le fichier filtered tau haplotypes du SR s'il exite
                        try:
                            file = pd.read_csv(f"{dir_gene}/{desman_file}/Filtered_Tau_star_haplotypes", sep ="\t", header=None)
                            ##recupére dans file les lignes du filtered tau dont la colonne 1 correspond aux positions variable
                            file = file[file.iloc[:,1].isin(x)]
                            # print(file)
                            ##Concatener toutes les lignes dans haplo_file
                            haplo_file = pd.concat([haplo_file, file], ignore_index=True)
                            ##Ne garder que la derniére colone du fichier haplo_file
                            new_haplofile = haplo_file.iloc[:,2:]

                        ## S'il n'existe pas : ecrit qu'il n'existe pas et passe à l'étape suivante   
                        except FileNotFoundError :
                            print(f"{gene}/{desman_file}/Filtered_Tau_star_haplotypes not found")
                            

                    ##Liste du nombre de colonnes (variants) par fichier
                    list_var = [i for i in range(1, len(new_haplofile.columns)+1)]
                    # print(list_var)

                    ## Si plus de 2 colonnes par fichier :
                    if len(new_haplofile.columns)> 1:
                        ##créer des noms de sous colonnes suivant le nombre de colonnes
                        new_haplofile.columns = pd.MultiIndex.from_product([[desman_file], list_var])
                    else:
                        ##sinon grader comme noms de colonnes le nom du fichier desman
                        new_haplofile.columns = [desman_file]
                        # print(haplo_file)

                    ##Concatener les colonnes de chaque fichier avec le fichier resulte en colonne
                    result = pd.concat([result, new_haplofile], axis=1)
                    # print(result) 
                    result.to_csv(f"{dir_gene}/haplo_{i}Variant.csv")


        ###Concatener en ligne tous les fihier d'un même nombre de variant en un fichier

        ##Chemin vesr les fichiers haplo
        haplo_path = glob.glob(f"{dir_gene}/haplo_*")

        ##Ouvrir un fichier pour écriture
        with open(f"{dir_gene}/all_haplo_{gene}.csv", "w") as outputfile:
            outputfile.write(f"Gene : {dir_gene}" + "\n")
            for file in haplo_path:
                with open (file, "r") as haplo_file:
                    outputfile.write(haplo_file.read() + "\n")


print('process done successfully!!') 
#END MAIN 
     
if __name__=="__main__": 
    main() 
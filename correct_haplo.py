# # Librairies
import os
import pandas as pd
from argparse import ArgumentParser
import sys
import fnmatch
import glob
import numpy as np


######################################################################################################################
# ## Infos 

"""Script pour confirmer les haplotypes correctement reconstruits par DESMAN"""

#Lancement du script
# python3 correct_haplo.py -i directory -f freq_allellic -g gene-list.txt - r "number_replicate" -s "threshold_replicate"

#####################################################################################################################



# Analyseur d'argument

def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="directory", help="Desman directory")
    parser.add_argument("-f", "--freq_allelic", dest="freq_dir", help="freq allelic OSI direcetory")
    parser.add_argument("-g", "--gene_name", dest="gene_name", help="list of gene names")
    parser.add_argument("-r", "--repic", dest="nb_replicate", default="", help="number of replicats use in DESMAN")
    parser.add_argument("-s", "--seuil", dest="threshold_replicate", help="minimal seuil of haplotype redondancy in replicate")

    options = parser.parse_args()
    if len(sys.argv)< 5:
        sys.exit("warning : wrong number of arguments!!")
    return options.directory, options.freq_dir, options.gene_name, options.nb_replicate, options.threshold_replicate





# recupérer le nombre de variant optimal du fichier haplo
def optimal_variant(haplo_file):
    numb_OptVar= int(haplo_file.split("_")[1].replace("Variant.csv",''))
    return numb_OptVar



# recupérer les variants majoritaires sous forme de liste et les positions variables de chaque fichier haplo
def maj_variants(haplo_file_path):

    read_file=pd.read_csv(haplo_file_path, index_col=[0])
    var_pos = read_file["Position"]
    read_file= read_file.drop(columns=["Position"])
    read_file.columns= range(1, len(read_file.columns)+1)
    variants=[]
    for col in read_file.columns:
        variant= read_file[col].tolist()
        variants.append(variant)
    return variants, var_pos





def main():

    directory, freq_dir, gene_name, nb_replicate, threshold_replicate= config_parameters()
    
    ##Lire la liste des noms de gènes
    gene_names = pd.read_csv(gene_name, index_col=[0], header=None)
    
    ##Pour chaque gène de la liste des génes analysés par DESMAN
    nb_gene=0
    for gene in gene_names.index:
        print(gene)
        ##DESMAN gene directory
        
        gene_dir = os.path.join(directory + "/DESMAN_" + gene)
        # Fichier de sortie
        out_path = os.path.join(gene_dir, "correct_haplotype.csv")

        out_file = open(out_path, "w")


        ##Lister les fichier du répertoire
        nb_haplo=0
        for haplo_file in os.listdir(gene_dir):
           
            #Si le pattern du fichier est haplo_*
            if fnmatch.fnmatch(haplo_file, "haplo_*"):

                haplo_file_path = os.path.join(gene_dir, haplo_file)

                variants , var_pos = maj_variants(haplo_file_path)
            

                """ Step_1 : ne pas analyser les fichiers de variants incomplets """

                if (optimal_variant(haplo_file)*int(nb_replicate)) == len(variants):
                    # print(f"{haplo_file} bon nombre de variant")
                    # print(len(variants))


                    """ Step_2 : Garder les haplotypes présents au moins dans 80% des réplicats  """
                
                    variants_not_dup = [list(x) for x in set(tuple(x) for x in variants)]
                
                    variant_file=pd.DataFrame(var_pos)

                    for variant in variants_not_dup:
                        # print(variant)
                        nb_dup=0
                        for var in variants:
                            # print(var)
                            if variant == var:
                                nb_dup+=1
                        
                        if nb_dup >= int(threshold_replicate):
                            # variant = pd.Series(variant, index= var_pos)
                            variant_file =pd.concat([variant_file, pd.Series(variant)], axis=1)
                            variant_file.to_csv(f"{gene_dir}/{haplo_file}")
                            
                    # print(variant_file)
                    nb_haplo+=1

                
                else :
                    os.remove(f"{gene_dir}/{haplo_file}")

        """ Step_3 : créer un fichier des variants uniques et majoritaires  """

        final_variant=pd.DataFrame()
        
        exit_loop = False
        for i in range(1, nb_haplo+1):
            j= i+1

            try:
                variants_1, var_pos_1 = maj_variants(f"{gene_dir}/haplo_{i}Variant.csv")
                variants_2, var_pos_2 = maj_variants(f"{gene_dir}/haplo_{j}Variant.csv")
      

                for var in variants_1:
                    if var in variants_2 and len(variants_2) > len(variants_1):
                        exit_loop = False

                    else:
                        for var1 in variants_1:
                            final_variant = pd.concat([final_variant, pd.Series(var1)], axis=1)
                            final_variant.to_csv(f"{gene_dir}/final_variant.csv")
                        exit_loop = True
                        break
            
            except FileNotFoundError :
                for var1 in variants_1:
                    final_variant = pd.concat([final_variant, pd.Series(var1)], axis=1)
                    final_variant.to_csv(f"{gene_dir}/final_variant.csv")
                exit_loop = True
                
            if exit_loop:
                break


        """ Step_4: vérifier si les variants majoritaires sont correctements reconstruits  """

        # Ouvrir les tables de comptages
        freq_repository = glob.glob(f"{freq_dir}/*")
        final_variant.columns = range(1, len(final_variant.columns)+1)

        nb_var=0
        for col in final_variant.columns:

            nb_base=0
            index_var=0
            for base in final_variant[col].tolist():

                for freq in freq_repository:
                    try:
                        comptage = pd.read_csv(f"{freq}/{gene}.freq", sep="\t", index_col=[0])
                        if var_pos[index_var] in comptage.index:
                            if comptage.loc[var_pos[index_var],base] != 0:
                                nb_base+=1
                                index_var+=1
                                break
                            
                    except FileNotFoundError :
                        print(f"{freq}/{gene}.freq not exist or position not exit in file")

            if nb_base == len(final_variant.index):
                nb_var+=1
               

            else :
                final_variant=final_variant.drop(columns=col, axis=1)

        if nb_var == len(final_variant.columns):
            # print("All variant correct")
            final_variant.to_csv(f"{gene_dir}/final_variant.csv")
            print(final_variant)
            nb_gene+=1


print("Process done succesfully!!")

if __name__ == "__main__":
    main()
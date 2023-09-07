# # Librairies
import os
import pandas as pd
from argparse import ArgumentParser
import sys
import fnmatch
import glob
import numpy as np



###################################################################################################
## Infos :

""" Script pour retrouver le réplicat avec les variants majoritaires corectement reconstruits """

## Lancement du script:

# python3 find_replicate.py -i "directory" -g "gene_name"

###################################################################################################




# analyseur d'argument
def config_parameters():
    parser= ArgumentParser()
    parser.add_argument("-i", "--input", dest="directory", help= "Desman directory")
    parser.add_argument("-g", "--gene", dest="gene_name", help="list of gene_name")
    options=parser.parse_args()

    if len(sys.argv) < 2:
        sys.exit("warning: wrong number of arguments!!")
    
    return options.directory, options.gene_name


# recupérer les variants majoritaires sous forme de liste et les positions variables de chaque fichier haplo
def read_replicate_haplotype(replicate_file_path):
    read_file = pd.read_csv(replicate_file_path, sep="\t", header=None)
    read_file = read_file.loc[:,2:]
    haplo_list = read_file.values.tolist() 
    
    return haplo_list


def main():

    directory, gene_name = config_parameters()

    out_file = open(f"{directory}/replicates.csv", "w")

    gene_names = pd.read_csv(gene_name, index_col=[0])

    nb_gene=0
    for gene in gene_names.index:

        gene_dir = os.path.join(directory + "/DESMAN_" + gene)

        # Ouvrir le fichier des variants majoritaires pour lecture
        read_variant = pd.read_csv(f"{gene_dir}/final_variant.csv", index_col=(0))
        variant_list = read_variant.values.tolist()
        optimal_var = len(read_variant.columns)

        #Ouvrir les fichiers filtered_tau_satar_haplotypes pour rechercher les variants majoritaires correctes
        for replicate_file in os.listdir(gene_dir):
            if fnmatch.fnmatch(replicate_file, f"desman.{optimal_var}.*"):

                replicate_file_path = os.path.join(gene_dir, replicate_file, "Filtered_Tau_star_haplotypes")
                haplo_list = read_replicate_haplotype(replicate_file_path)

                nb_var=0
                for var in variant_list:
                    if var in haplo_list:
                        nb_var+=1

                if nb_var == optimal_var :
                    out_file.write(f"{gene}\t{replicate_file}\n")
                    nb_gene+=1
                    break
                    

                    
print("Process done sucessfully !!")

if __name__ == "__main__":
    main()


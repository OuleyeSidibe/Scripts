# # Librairies
import os
import pandas as pd
from argparse import ArgumentParser
import sys
import fnmatch
import glob



#################################################################################
# ## Infos 

"""Script pour confirmer les haplotypes correctement reconstruits par DESMAN"""

#Lancement du script
# python3 correct_haplo.py -i directory -f freq_allellic -g gene-list.txt

#################################################################################



# Analyseur d'argument

def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="directory", help="Desman directory")
    parser.add_argument("-f", "--freq_allelic", dest="freq_dir", help="freq allelic OSI direcetory")
    parser.add_argument("-g", "--gene_name", dest="gene_name", help="list of gene names")
    parser.add_argument("-r", "--repic", dest="nb_repic", default="", help="number of replicats use in DESMAN")
    options = parser.parse_args()
    if len(sys.argv)< 4:
        sys.exit("warning : wrong number of arguments!!")
    return options.directory, options.freq_dir, options.gene_name, options.nb_repic





# recupérer le nombre de variant optimal du fichier haplo
def optimal_variant(haplo_file):
    numb_OptVar= int(haplo_file.split("_")[1].replace("Variant.csv",''))
    return numb_OptVar


def main():

    directory, freq_dir, gene_name, nb_repic= config_parameters()
    
    ##Lire la liste des noms de gènes
    gene_names = pd.read_csv(gene_name, index_col=[0], header=None)
    
    ##Pour chaque gène de la liste des génes analysés par DESMAN
    for gene in gene_names.index:
        print(gene)
        ##DESMAN gene directory
        
        dir_gene = os.path.join(directory + "/DESMAN_" + gene)
        # Fichier de sortie
        out_path = os.path.join(dir_gene, "correct_haplotype.csv")

        out_file = open(out_path, "w")


        ##Lister les fichier du répertoire
        for haplo_file in os.listdir(dir_gene):
            nb_var=0
            # print(haplo_file)
            
            #Si le pattern du fichier est haplo_*
            if fnmatch.fnmatch(haplo_file, "haplo_*"):

                ##eliminer les fichiers haplo avec un nombre de variant incomplet
                read_haplo =pd.read_csv(f"{dir_gene}/{haplo_file}", index_col=[0])
                haplo = read_haplo.iloc[:, 1:]
                length_haplo = len(haplo.columns)
                

                # numb_OptVar= int(haplo_file.split("_")[1].replace("Variant.csv",''))
                # print(type(numb_OptVar))
                # print(type(nb_repic))
                # nb_varDES = int(numb_OptVar) * int(nb_repic)
                # print(nb_varDES)
                # print(optimal_variant(haplo_file))
                # print(length_haplo)
                print(type(nb_repic))
                print(optimal_variant(haplo_file)*int(nb_repic))
                if (optimal_variant(haplo_file)*int(nb_repic)) == length_haplo:
                    print(f"{haplo_file} bon nombre de variant")
                
                
                    #vecteur vide pour stocker le nombre de bonnes pos var par haplotypes
                    nb_Haplo=0 

                    ##Recupérer les positions variables dans une liste
                    var_pos = read_haplo.Position
                    ## Recupérer le nombre de position variable du géne
                    nb_varpos = len(var_pos)
                    
                    ##Parcourir chaque position variable du gène
                    for i in range (0,nb_varpos):
                        ##Pour chaque position position variable
                        var_pos = read_haplo.iloc[i,0]
                        ##Recupérer les bases uniques pour la position 
                        df_unique = read_haplo.iloc[i,1:].unique()

                        ##Si le nombre de bases unique de la position est inférieur ou égale à 2
                        if len(df_unique) <= 2:

                            ##Lister les chemin des sous repertoires du dossier freq
                            freq_dossier = glob.glob(f"{freq_dir}/*/*.freq")
                            
                            ##Créer un df vide 
                            posVar_file=pd.DataFrame()

                            ## Pour chaque sous repertoire
                            for freq in freq_dossier:
                                ##Recupérer le nom d'échantillon
                                sample_name = os.path.basename(os.path.dirname(freq))
                                ##Définir le pattern du fichier
                                pattern = (f"{freq_dir}/{sample_name}/{gene}.freq")
                                #Si le fichier.freq du géne est présent dans l'échantillon
                                if freq == pattern :
                                    ##Essaye d'ouvrir le fichier pour lecture
                                    try:
                                        read_freqfile = pd.read_csv(pattern, sep="\t", index_col=[0])
                                        ##Recupérer les comptages de la position variable
                                        read_varpos = read_freqfile.loc[var_pos,:]
                                        ##Concatner tous les comptages de la position var suivant les différents échantillons où le géne est retrouvé
                                        posVar_file = pd.concat([posVar_file, read_varpos], axis=1)
                                        ## Si la position variable est absente
                                    except:
                                        #Sinon imprimer :
                                        print(f"{var_pos} not found in {sample_name}")
                                
                            ##Pour plus d'une position variable essayer:
                            try:
                                ##Somme des comptages de la premiére base
                                varpos_1 = sum(posVar_file.loc[df_unique[0],:])
                                ##Somme des comptages de la seconde base
                                varpos_2 = sum(posVar_file.loc[df_unique[1],:])
                                
                                ##Si la somme des comptages des deux bases uniques est différent de 0 pour la mm position
                                if varpos_1 != 0 and varpos_2 != 0:
                                    ##L'haplotype est jugé correcte pour cette position
                                    nb_Haplo += 1
                                ## Si la somme des comptages d'une des bases unique est =0
                                elif varpos_1 == 0 or varpos_2 == 0:
                                    ##L'haplotype est jugé incorrecte pour cette position
                                    nb_Haplo += 0
                                
                                
                            except:
                                varpos_1 = sum(posVar_file.loc[df_unique[0],:])
                                # print(pos_1)
                                if varpos_1 != 0 :
                                    ##L'haplotype est jugé correcte pour cette position
                                    nb_Haplo += 1

                                elif varpos_1 == 0 :
                                    ##L'haplotype est jugé incorrecte pour cette position
                                    nb_Haplo += 0
                        ## Sinon
                        else:
                            ##L'haplotype est jugé incorrecte pour cette position
                            nb_Haplo += 0

                    ##Si le nombre d'haplotye jugé correcte est égale au nombre de position 
                    if nb_Haplo == nb_varpos:
                        #L'haplotype est jugé corecte pour toutes les positions
                        out_file.write(f"{haplo_file} \n")    
                    
                    ##Position variable suivante
                    nb_var +=1
                else:
                    print(f"{haplo_file} mauvais nombre de variant")
                


        
            
print("Process done succesfully!!")

if __name__ == "__main__":
    main()

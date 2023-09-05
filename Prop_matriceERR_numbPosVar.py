import pandas as pd
import numpy as np
import glob
import os
import fnmatch
import sys
from argparse import ArgumentParser		

#Lancement script
#python3 Prop_matriceERR_numbPosVar.py -i {directory} -p {proportion minimum} -g {list des noms de gènes dans directory}


#Analyseurs d'arguments
def config_parameters():
    parser=ArgumentParser()
    parser.add_argument("-i", "--input", dest="directory",default="", help="DESMAN repertory")
    parser.add_argument("-p", "--proportion", dest="proportion", default="", help="minimum proportion")
    parser.add_argument("-g", "--genes", dest="gene_names", default="", help="list of gene_names")
    options =parser.parse_args()
    if len(sys.argv)<3:
        sys.exit("Warning:wrong number of arguments!!")
    return options.directory, options.proportion, options.gene_names




def main():

    directory, proportion, gene_names = config_parameters()

    #Ouvrir gene_names pour lecture
    genes_list = pd.read_table(gene_names, header=None, index_col=[0])

    table= pd.DataFrame(columns=["gene_name", "Num_PosVariable", "Prop_matriceERRexact"])



    for gene in genes_list.index:
        list =[]
        df_sel_var = pd.read_csv(f"{directory}/DESMAN_{gene}/sel_var.csv")
        #Nombre de positions variables du géne
        numb_Pvar = len(df_sel_var)
        
        ##ourvrir chaque sous repertoire du géne
        #Répertoire desman du géne        
        dir_gene = os.path.join(directory, "DESMAN_" + gene)
        #Nom de fichier de la matrice d'erreur
        matriceERR_file = "Eta_star.csv"
        #Pour chaque sous repertoire
        
        for desman_file in os.listdir(dir_gene):
            #Joindre le fichier desman au chemin du repertoire parent
            path = os.path.join(dir_gene, desman_file)
            #renvoie le chemin que si c'est un repertoire
            # if os.path.isdir(path) :
            if fnmatch.fnmatch(desman_file, "desman*"):
                #rajoute le chemin vers le fichier matrice d'erreur
                dir_matriceERR = os.path.join(path, matriceERR_file)
                #print(dir_matriceERR)
                #print(dir_matriceERR)
                
                ##Ouvrir le fichier pour lecture
                df_matriceERR = pd.read_csv(dir_matriceERR)
                #print(df_matriceERR.shape)
                df = np.genfromtxt(dir_matriceERR , delimiter=",", dtype=float, skip_header=True, usecols=(1,2,3,4))
                #print(df)
                prop = float(proportion)

                if np.array((df[0,0] > prop) & (df[1,1] > prop) & (df[2,2] > prop) & (df[3,3] > prop)):
                    #print("yes")
                    list.append(1)
                    
                else :
                    #print("no")
                    list.append(0)
                  
        Prp_matriceERR = sum(list)/ len(list)
     
        
        
        #result = pd.DataFrame({"gene_name":gene, "Num_PosVariable":numb_Pvar,"Prop_matriceERRexact":Prp_matriceERR})
        
        result = pd.DataFrame(data=[[gene, numb_Pvar, Prp_matriceERR]], columns=["gene_name","Num_PosVariable", "Prop_matriceERRexact"])
        table = pd.concat([table, result], ignore_index=True)
        
        
    table.to_csv(f"{directory}/Prop_matriceERRexact.csv", index=False)
        
print('process done successfully!!') 
#END MAIN 
     
if __name__=="__main__": 
    main() 
    
    
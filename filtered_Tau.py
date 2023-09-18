##librairies
import pandas as pd
import os
import sys
import fnmatch
from argparse import ArgumentParser





################################################################################################
# ## Infos 

"""Script pour reconstruire le fichier filtered tau avec uniquement les positions variables"""

#Lancement du script
# python3 filtered_Tau.py -i directory -f freq_alle

################################################################################################



#analyseur d'argument
def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="directory", help="DESMAN directory")
    options = parser.parse_args()
    if len(sys.argv) < 1:
        sys.exit("warning : wrong number of arguments!!")
    return options.directory


def main():

    directory = config_parameters()

    read_selvar = pd.read_csv(f"{directory}/file.outsel_var.csv", index_col=[1])
    pos_var = read_selvar.index

    for desman_file in os.listdir(directory):
        if fnmatch.fnmatch(desman_file, "desman.*"):
            # desman_file_path = os.path.join(directory, desman_file)
            read_filtered_Tau = pd.read_csv(f"{directory}/{desman_file}/Filtered_Tau_star.csv", index_col=[0])
            # print(read_filtered_Tau)
            # Créer un nouveau DataFrame avec les lignes à conserver
            new_filtered_Tau = read_filtered_Tau[read_filtered_Tau["Position"].isin(pos_var)]
            
            print(new_filtered_Tau)
            new_filtered_Tau.to_csv(f"{directory}/{desman_file}/Filtered_Tau_star.csv")

    
                
    print("Process done succesfully!!")

if __name__ == "__main__":
    main()








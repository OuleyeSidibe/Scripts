#!usr/bin/python


#run script
#python3 script_new_freq_allelic.py -i inputdir -o outputdir -g stat_afterFILTER.csv


#Librairies
from argparse import ArgumentParser							
import pandas as pd
import glob
import os
from os import walk
import numpy as np


#Configuration de l'analyseur d'arguments
def config_parameters(): 																								#Fonction de définition des paramétres de config
	parser = ArgumentParser() 																						    #Création de l'analyseur "Parser" qui va exéxuter un prog python avec une liste d'argument
	parser.add_argument("-i", "--input", dest="inputdir",default="", help="input directory - allelic_frequencies") 		#Ajout des arguments "optionnels", nom de l'attribut "dest", bréve description de l'argument "help"
	parser.add_argument("-o", "--output", dest="outputdir",default="", help="output file") 
	parser.add_argument("-g", "--genes", dest="genes", help="retained genes after filter")
	args = parser.parse_args() 	
	return args.inputdir, args.outputdir, args.genes		 					
								

#Main
def main():

    inputdir, outputdir, genes = config_parameters()

    my_data = pd.read_csv(f"{genes}", sep = "\t" , usecols = [0,1], names=["gene_name", "sample_name"], index_col= [0], squeeze=True)
        
    dossier_freq = glob.glob(f"{inputdir}/*")

    for freq in dossier_freq:

        nb_sample = 0
        sample_name = os.path.basename(freq)
        os.mkdir (f"{outputdir}/{sample_name}")
        gene_list = my_data[my_data == sample_name].index
        # print(gene_list)
        for gene in gene_list :
            gene_name = gene
            # print(gene_name)  
            comptage_freq = pd.read_csv(f"{inputdir}/{sample_name}/{gene_name}.freq", sep ='\t', index_col= "coordinate", dtype={'coordinate': np.int32, 'A': np.int32, 'C': np.int32,'G': np.int32, 'T': np.int32} ) 
            comptage_mpileup = pd.read_csv(f"{inputdir}/{sample_name}/{gene_name}.mpileup", index_col=[0]) 
            # print(comptage)
            comptage_freq.to_csv(f"{outputdir}/{sample_name}/{gene_name}.freq", sep='\t')
            comptage_mpileup.to_csv(f"{outputdir}/{sample_name}/{gene_name}.mpileup")
        nb_sample+=1

if __name__=="__main__": 
    main() 


                



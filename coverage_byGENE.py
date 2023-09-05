#!usr/bin/python


#run script
#python3 test_coverage.py -i inputdir -o outputfile -n gene_size.txt





#Librairies
from argparse import ArgumentParser							# lib d'analyseur d'arguments
import random, sys, random , string, csv, os, math			# Lib générateur de nbr aléatoire, arg a fournir en script, opérations sur des chaines, lecture et écriture sur fichiers csv, accés aux fonctions mathématiques
import numpy as np											# Lib de manipulation de fonctions math opérant sur des tableaux ou matrices 
import glob													# Lib de renvoie de ts les chemins de fichiers correspondant à un modéle
import pandas as pd 	    						        # Lib de structurations de données et de maipulations de tableaux numériques




#Configuration de l'analyseur d'arguments
def config_parameters(): 																								#Fonction de définition des paramétres de config
 
	#use = "use : ./test_coverage.py -i freq_dir -o out_dir -w window -g genome" 							            #Description de l'utilsation du programme 
	parser = ArgumentParser() 																						    #Création de l'analyseur "Parser" qui va exéxuter un prog python avec une liste d'argument
	parser.add_argument("-i", "--input", dest="inputdir",default="", help="input directory - allelic_frequencies") 		#Ajout des arguments "optionnels", nom de l'attribut "dest", bréve description de l'argument "help"
	parser.add_argument("-o", "--output", dest="outputfile",default="", help="output file") 							#type :évite d'interpréter les arguments comme des chaines de caractéres "par defaut" mais en entier
#	parser.add_argument("-w", "--windows", dest="windows",type=int, help="windows length") 
#	parser.add_argument("-x", "--interval", dest="interval", type=int, help="interval between begining of windows") 
	parser.add_argument("-n", "--names", dest="names", help="index with names and length of genes") 					
	options = parser.parse_args() 																						#Analyse des arguments
	if len(sys.argv)< 3 :												
		sys.exit("Warning : wrong number of arguments !!")																#Si lors du lancement le nbr d'arg est inf à 5, qiutter le programme et afficher le message d'erreur spécifié	
	return options.inputdir, options.outputfile, options.names															#Analyse des arg ajoutés	

#END config_parameters 



#Fichiers de comptages
def open_inputfile(inputfile):																															#Ouverture du fichié d'entré pour lecture : fichier de comptages

	try:																																				#Instruction d'interception des exceptions
		comptage = pd.read_csv(inputfile, sep='\t', dtype={'coordinate': np.int32, 'A': np.int32, 'C': np.int32,'G': np.int32, 'T': np.int32} )			#Lecture de la table de comptage, contenant des entiers "int32"(4octets(1octet=8bits))
	except IOError:																																		#si rencontre d'une exception spécifiée dans "except"
		print(f"{inputfile} not found")																													#imprimer "nom du fichier" suivie de "not found" puis quitter
		sys.exit()																															
	return comptage																																		#sinon retourner la variable contennant les données lus




#Fichier de coordonnées des gènes
def read_names(gene_size): 																																#Lecture des noms des génes et def de leurs tailles
	try:																																				
		list_of_gene = pd.read_csv(gene_size, sep='\t', dtype={'begin': np.int32, 'end': np.int32, 'name': object})
		#print(list_of_gene.to_string())
	except IOError:																																		#si rencontre d'une exception spécifiée dans "except"
		print(f"{gene_size} not found")																													#imprimer "nom du fichier" suivie de "not found" puis quitter
		sys.exit()																			
	return list_of_gene																																	#retourner la variable contenant la liste des génes et leurs coordonnées 




#Main
def main(): 
	inputdir,  outputfile, genes_list = config_parameters()                    											#Appel de fonction : parametres de config des arguments

	

	dossiers_freq = glob.glob(f"{inputdir}/*/")																			#liste des répertoires "noms d'échantillons" présent dans le directory "freq_ala"

	list_of_gene = read_names(genes_list)																				#Lecture du fichier "genes_list" et le stocké sous forme de liste dans list_of_gene)
	
	for gene in range(list_of_gene.size):																				#Pour chaque gène de la liste des gènes (771)("size":donne un entier représentant le nbr d'elts de l'objet "gene_list" )											
			nb_sample=0																									#Pour un nbr d'échantillon égale à 0 :
			#print(f"{gene}")

			result = pd.DataFrame(columns=['coordinate'], data=range(list_of_gene.iloc[gene]['begin'], list_of_gene.iloc[gene]['end']+1), index=range(list_of_gene.iloc[gene]['end'] -list_of_gene.iloc[gene]['begin']+1))	 # Créer un df des cordonnées pour chaque gene 
			for freq in dossiers_freq: 																																						#Pour chaque elt (échantilon ou ligne) de dossier_freq,
				sample_name=os.path.basename(os.path.dirname(freq)) 																														#Nom de chaque échnatillon (basename : sous repertoire et dirname : repertoire)
				gene_name=list_of_gene.iloc[gene]['name']																																	#nom d chaque gène de la liste des gènes
				#print(gene_name)
				try : 
					df1 = pd.read_csv(f"{inputdir}/{sample_name}/{gene_name}.freq", sep='\t', dtype={'coordinate': np.int32, 'A': np.int32, 'C': np.int32,'G': np.int32, 'T': np.int32} )	#Lecture de chaque fichier de gènes
					df1[sample_name]=df1.iloc[:,1:5].sum(axis=1) 							                                                                                                #Somme des 4 colones (A G T c) dnas une colonne prenant le nom de l'échantillon
					## merge original data frame + new dataframe 
					result = pd.merge(result, right= df1[['coordinate',sample_name]], on="coordinate", how="left")																			#merger les 2 df result sur la colonne 'coordinate de gauche 'result'
					result = result.set_index('coordinate')																																	#retirer l'index				
					nb_sample+=1																																							#passer à l'element "freq" suivant (échantillon)
				except:																																										#si le nom de gene n'est pas présent dans l'echnatillon
					print(f"no {gene_name} for sample {sample_name}")																														#affiche "nom du gène" pas présent  dans "nom de l'echnatillon"s

			if nb_sample >0:	
				# nb_ofSAMPLES = result.columns
				# nb_ofSAMPLES = len(nb_ofSAMPLES)
				# nb_ofSAMPLES																																							#Si le nombre d'échantillon est sup à 0						
#				result = result.fillna(0)																																					#remplacer les Na par 0
#				print(f"{gene_name}\n{result}")																																				#afficher le nom du gene puis retourner à la ligne le resultat
				# result.to_csv(f"{outputfile}/{gene_name}__{nb_ofSAMPLES}")		
				result.to_csv(f"{outputfile}/{gene_name}")																																#fichier csv du résultat avec le nom du géne dans le outputfile
				
	
	print ('process done successfully!!') 
#END MAIN 
     
if __name__=="__main__": 
    main() 
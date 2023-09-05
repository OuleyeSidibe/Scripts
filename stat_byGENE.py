#!usr/bin/python

#run script 
#stat_byGENE.py -d <depth> -c <coverage_Depth> -o outputdir/stat_byGENE


########
#Library
########

import pandas as pd
import numpy as np
import glob
import os
from argparse import ArgumentParser


#############################
#Definition of path directory
#############################

#Directory path of inputfiles
inputdir= "/home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/coverage_by_SAMPLE/out_coveragebyGENE"


def config_parameters():
	parser = ArgumentParser()
	parser.add_argument("-d", "--depth", type= int ,dest="depth", help="Minimum number of readings per base ")
	parser.add_argument("-c", "--coverage_Depth",dest ="coverage_Depth" ,type = int, help="percentage coverage of the size of the mapped area where the depth is applied")
	parser.add_argument("-o", "--outputfile",dest ="outputfile" , help="stats for each gene of inputdir")
	args = parser.parse_args()
	return args.depth, args.coverage_Depth, args.outputfile

#######################################
#Function for reading input files genes 
#######################################

def read_gene_file(inputfile):
    gene_files = pd.read_csv(inputfile, header=0, sep =",", index_col=0,  dtype= {'coordinate' : np.int32})

    return(gene_files)


############################################
#Function for coverage and depth calculation
############################################

def main ():
	depth, coverage_Depth, outputfile  = config_parameters()

    #Put inputfile genes into a list
	gene_files_list = glob.glob(f"{inputdir}/*")
    
    #Print list of inputfile genes
	#print(f"{gene_files_list}")

    #For each gene of genes
	for gene in gene_files_list:
        
		nb_gene=0
		gene_name = os.path.basename(gene)                                                      			#gene name correspond to the basename of gene file
		#print(f"{gene}")                                                                        

		df = pd.read_csv(gene, header=0, sep =",", index_col=0, dtype= {'coordinate' : np.int32})           #read csv gene file
       
		df1 = df
		df1.loc["mean"] = df.mean(axis =0)
		df1.loc["std"] = df.std(axis=0) 
#		print(df1.loc["mean", "CTU_ACR"])

		df2 = df.iloc[:-2]
		df2 = df2.apply(pd.value_counts)                                                          			#Apply to each columns of df a count of each value 
        
		df2.loc["COVbase"]= df2[:].sum()                                                          			#Add a row named "COVbase" whose correspond to a sum of each valu for each columns
        
		df2.loc["Coverage (%)"] = round((df2.loc["COVbase"] / len(df))*100)										#Add a row named "Coverage" which correspond to porcentage of coverage of gene in each sample(columns)

		df3 = df2.iloc[depth:-2,]                 					                            			#for depth calculate, take only values >= "3.0", so don't take the two first index "1.0 and2.0" of df
		df2.loc["Depthbase"] =df3[:].sum()                                                           		#Add depth row df

		base_coverage = ((coverage_Depth*df2.loc["COVbase"])/100)											#The number of base corresponding to 90% of a gene size
		df2.loc["Depth_per_coverage"] = (df2.loc["Depthbase"] - base_coverage)								#A negatif vallue it's mean that 90% of base of gene have not at less 3 reads.
		

		df2 = df2.iloc[-4:]   
		result = pd.concat([df1, df2]) 
		result = result.fillna(0)
		result = result.astype(int)

        
		result = result.to_csv(f"{outputfile}/{gene_name}")                                  #outpute file as csv in output directory and gene_name as file name
        #print(result)
		nb_gene+=1                                                                              			 #Run next gene file
 

if __name__=="__main__": 
    main() 
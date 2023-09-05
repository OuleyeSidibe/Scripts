#!usr/bin/python


#run script 
#script/table_stat_forGENES.py -c <coverage> -i inputdir/stat_byGENE -o outputdir


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

#csv table before filters
# outputfile_1 = $outputdir/stat_beforeFILTER.csv

#csv table after filters
# outputfile_2 = "/home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/coverage_by_SAMPLE/stat_afterFILTER.csv"
# outputfile_1 = os.path

# fileout_1 = open(outputfile_1, "w")
# fileout_2 = open(outputfile_2, "w")





def config_parameters():
	parser = ArgumentParser()
	parser.add_argument("-c", "--coverage",type = int, dest="coverage", help = "percentage of total coverage of the gene")
	parser.add_argument("-i", "--input", dest="inputdir", help = " directory of all genes stats")
	parser.add_argument("-o", "--output", dest="outputdir", help = " directory of tables for all genes stat before and after coverage and depth filter")
	args = parser.parse_args()
	return args.coverage, args.outputdir, args.inputdir



############################################
#Function for coverage and depth calculation
############################################


def main ():
	coverage, outputdir, inputdir = config_parameters()

    #Put inputfile genes into a list
	gene_files_list = glob.glob(f"{inputdir}/*")
    
	outputfile_1 = (f"{outputdir}/stat_beforeFILTER.csv")
	outputfile_2 = (f"{outputdir}/stat_afterFILTER.csv")
	
	fileout_1 = open(outputfile_1, "w")
	fileout_2 = open(outputfile_2, "w")


	for gene in gene_files_list:

		nb_gene=0
		gene_name = os.path.basename(gene)                                                      #gene name correspond to the basename of gene file
		#print(f"{gene}")                                                                        

		df = pd.read_csv(gene, sep =",", index_col=0)
		#print(df)

		for col in df.columns:
			#print(col)
			fileout_1.write(f'{gene_name}\t{col}\t{df.loc["mean", col]}\t{df.loc["std", col]}\t{df.loc["Coverage (%)", col]}\n')
			# fileout_1.write(f'{gene_name}\t{col}\t{df.loc["mean", col]}\t{df.loc["std", col]}\t{df.loc["Coverage (%)", col]}\t{df.loc["Depth_per_coverage", col]}\n')
				

			if df.loc["Depth_per_coverage", col] >= 0 and df.loc["Coverage (%)", col] >= coverage:
				# fileout_2.write(f'{gene_name}\t{col}\t{df.loc["mean", col]}\t{df.loc["std", col]}\t{df.loc["Coverage (%)", col]}\t{df.loc["Depth_per_coverage", col]}\n')
				fileout_2.write(f'{gene_name}\t{col}\t{df.loc["mean", col]}\t{df.loc["std", col]}\t{df.loc["Coverage (%)", col]}\n')
				
			nb_gene+=1
      

if __name__=="__main__": 
    main() 
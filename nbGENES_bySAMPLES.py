#!usr/bin/python3


#Library
import pandas as pd
import numpy as np
import os
import glob

inputdir ="/home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/coverage_by_SAMPLE/out_coveragebyGENE"
outputfile = "/home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/coverage_by_SAMPLE/out_nbGENES_bySAMPLES.csv"

files = glob.glob(f"{inputdir}/*")
#print(f"{files}")

result = pd.DataFrame(columns= ["gene_name","nb_of_samples"])
#print(result)
for gene in files :
    nb_gene =0
    df = pd.read_csv(gene, header=0, sep =",", index_col=0,  dtype= {'coordinate' : np.int32})
    #print(df)
    gene_name = os.path.basename(gene)
    #print(gene_name)
    nb_ofSAMPLES = df.columns
    nb_ofSAMPLES = len(nb_ofSAMPLES)
    #print(nb_ofSAMPLES)
    #print(gene_name)
    tmp = pd.DataFrame(data=[[gene_name, nb_ofSAMPLES]],  columns=["gene_name","nb_of_samples"])
   # print(tmp)
    result = pd.concat([result, tmp])
    nb_gene+=1
    #print(result)
   
    result.to_csv(f"{outputfile}", index = False)
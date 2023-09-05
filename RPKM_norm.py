import pandas as pd
import glob
import sys
import os
import numpy as np
from argparse import ArgumentParser	

#run script
#python3 /RPKM_norm.py -i directory -g out_nbGENES_bySAMPLES.csv(or list of gene_names) -c freq_alleliques



#Config arguments analyseur

def config_parameters():
    parser = ArgumentParser()
    parser.add_argument("-i", "--inputdir", dest="inputdir", default="", help="directory" )
    parser.add_argument("-g", "--genes", dest="gene_names", default="", help="list of genes_names")
    parser.add_argument("-c", "--freq_allelic", dest="counting_table", default="", help="table des fréquence de chaque gène")
    options =parser.parse_args()
    if len(sys.argv)<3:
        sys.exit("Warning:wrong number of arguments!!")
    return options.inputdir, options.gene_names, options.counting_table



def main():

    inputdir, gene_names, counting_table = config_parameters()

##Create empty df with samples, readsSample data , and gene_name 
    #Open redsSAMPLEdf_gene_namesfor reading and sample name as index_col
    df_reads = pd.read_csv(f"{inputdir}/nb_readsSAMPLE.csv", index_col=[0], sep=",")
    #Open out_nbGENES_bySAMPLES.csv for reading and gene_name as index_col
    #gene_names = pd.read_csv(gene_names, index_col=[0], sep=",")
    gene_names = pd.read_csv(gene_names, index_col=[0])
    # print(gene_names)
    #create list of gene_names
    gene_namesList = gene_names.index
    # print(gene_namesList)
    #Create DataFrame with gene_names as colnames
    df_gene_names= pd.DataFrame(columns=gene_namesList)
    #concate 2 df
    coverageRPK = pd.concat([df_reads,df_gene_names])
    # print(coverageRPK)
   
    # print(coverageRPK.head())
    
## calculate RPK (Read per Kilobas) for each gene using counting_table
    
    fichier_freq = glob.glob(f"{counting_table}/*")
    #for each sample of counting_table repository
    for freq in fichier_freq:
        #Define sample_name
        sample_name = os.path.basename(freq)
        #for each gene of sample_name repository
        for gene in gene_namesList:
            try:
                #Open freq file of gene for reading
                gene_counting_table=pd.read_csv(f"{counting_table}/{sample_name}/{gene}.freq", sep='\t', dtype={'coordinate': np.int32, 'A': np.int32, 'C': np.int32,'G': np.int32, 'T': np.int32}, index_col=[0])
                #Make sum of each colone
                gene_counting_table.loc["sum"] = gene_counting_table.sum(axis=0)
                #calcul number of mapped reads
                mapped_reads = sum(gene_counting_table.iloc[-1])/151
                RPK = mapped_reads/(len(gene_counting_table)/1000)
                # print(RPK)
                #insert valus of mapped read in df
                coverageRPK.loc[sample_name, gene] = RPK
            except:
                #if gene doesn't exist in thi sample, put 0
                # print(f"no {gene} for sample {sample_name}")
                coverageRPK.loc[sample_name, gene] = 0

    #Convert float type as int type      
    # print(coverageRPK)      
    
    coverageRPK = coverageRPK.astype(int)
    #save table as coverageRPK
    coverageRPK.to_csv(f"{inputdir}/coverageRPK22.csv")  

print('process done successfully!!') 

if __name__=="__main__": 
    main() 
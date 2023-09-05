import pandas as pd 
import os
import glob
import fnmatch

directory = "/home/osidibe/work/fluGenAvi_fugace/snakemake/snk_DESMAN_G3"
pattern_DESfile = "desman*"


gene_name_dir = glob.glob(f"{directory}/*")

for gene in gene_name_dir:
    genes_rep = os.path.basename(gene)
    with open (f"{directory}/{genes_rep}/allfit_bis.txt", "w") as bis_allfit:
        for DEsfile in os.listdir(gene):
            if fnmatch.fnmatch(DEsfile, pattern_DESfile):
                fit_dir = os.path.join(f"{directory}/{genes_rep}/{DEsfile}/fit.txt")           
                with open (fit_dir, "r") as fitfile:
                    bis_allfit.write(fitfile.read())
                    
            

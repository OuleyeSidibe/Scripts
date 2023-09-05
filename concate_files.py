import pandas as pd
import glob
import os
import fnmatch


#directory
directory = "/home/osidibe/work/fluGenAvi_fugace/snakemake/snk_DESMAN_G9"


#create function for concatenate stderr_file and stdout_file

def concat_files(outputfile, pattern_file):
#Open out_files_stderr from writting
    with open(outputfile, "w") as outfile:
        #For each file of logs folder
        for filename in os.listdir(logs):
            #take file if he has pattern
            if fnmatch.fnmatch(filename, pattern_file):
                #Open file for reading
                with open(os.path.join(logs, filename), "r") as infile:
                    #write file into outputfile (stderr_file)
                    outfile.write(infile.read())
                    # os.remove(f"{logs}/{filename}")
                    

#list genes folder
genes_rep = glob.glob(f"{directory}/*")

#for each gene folder 
for gene in genes_rep:
    gene_name = os.path.basename(gene)
    logs =(f"{directory}/{gene_name}/logs")

    #out_files and patterns
    stderr_file = (f"{directory}/{gene_name}/logs/out_stder.txt")
    stderr_pattern = "*.stderr"
    stdout_file = (f"{directory}/{gene_name}/logs/out_stdout.txt")
    stdout_pattern = "*.stdout"
    
    #concate stderr files
    concat_files(stderr_file, stderr_pattern)
    concat_files(stdout_file, stdout_pattern)
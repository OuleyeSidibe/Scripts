#!usr/bin


#Import librairies
import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser 

#Paths
outputfile = "/home/osidibe/work/fasta_to_freq"
fasta_path = "/home/osidibe/work/fasta_to_freq/vat_E__3_AF153312.fasta" 




#open fasta file
with open(fasta_path) as handle :
    #SimpleFasterParser : For each record a tuple of two strings is returned, the FASTA title line (without the leading ‘>’ character), and the sequence (with any whitespace removed).
    for sequence in SimpleFastaParser(handle):
        seq = sequence[1]
        # print(seq)
        length_seq = len(seq)
        #create a df with one column "coordinate" whole data equal to the length of fasta sequence
        table = pd.DataFrame(columns=["coordinate"], data = range(1,length_seq+1,1), index=None)

        #Create a empty df whose columns names equal to uniq bases in fasta sequence
        df = pd.DataFrame(columns = ["A", "C", "G","T"])

        #Concatenate two df
        table = pd.concat([table,df])
        #convert float64 into int64
        table["coordinate"] = table["coordinate"].astype(np.int64)
        #set coordinate as index of new df
        table = table.set_index("coordinate")
        #replace NAN by 0
        table =table.fillna(0)
        
    #loop for each line of df (index of df)
    for index in range (1, length_seq+1):
        #if base in fasta sequence equal to A , column A in table is equal to 50 for this line of df others equal to 0
        if seq[index-1] == "A":
            table.loc[index,"A"] =50
        #if base in fasta sequence equal to C , column C in table is equal to 50 for this line of df others equal to 0
        elif seq[index-1] == "C":
            table.loc[index,"C"] =50
        #if base in fasta sequence equal to G , column G in table is equal to 50 for this line of df others equal to 0
        elif seq[index-1] == "G":
            table.loc[index,"G"] =50
        #if base in fasta sequence equal to T , column T in table is equal to 50 for this line of df others equal to 0
        elif seq[index-1] == "T":
            table.loc[index,"T"] =50

    #save table in outputfile as file.freq      
    table.to_csv(f"{outputfile}/fre_allelic_vat_E.freq", sep = "\t")
        
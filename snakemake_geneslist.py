import pandas as pd

outputfile = "/home/osidibe/work/fluGenAvi_fugace/snakemake/genes_list.txt"
data = "/home/osidibe/work/fluGenAvi_fugace/analyse_preDESMAN/filters/coverage_by_SAMPLE/out_statbyGENE_5_90_90/stat_afterFilter_final.csv"
file = open(outputfile, "w")

genes = pd.read_csv(data, sep = "\t")
genes_list = genes.gene_name.unique()

for gene in genes_list:
    file.write(f"{gene}")

print(genes_list)
print(len(genes_list))

import os
import sys
import numpy as np
import pandas as pd

input_folder = sys.argv[1]
output_file = sys.argv[2]

# Import all STAR output files and append them in results_l list:
results_l = list()

for file in os.listdir("../{}/".format(input_folder)):
    if file.endswith("ReadsPerGene.out.tab"):
        results_l.append(pd.read_csv(
                  "../{}/".format(input_folder)+file, 
                  sep="\t", # tab as separater
                  skiprows=4, # Skip first 4 rows
                  names = ['gene_id', file.split("_ReadsPerGene")[0]], # name columns "gene_ID" and sample name
                  usecols = [0,1])) # Only use first two columns of STAR output. (Second column is gene counts)
        
# Set "gene_id" column as index for all dataframes in list (this is required for concatenating)
results_l = [df.set_index("gene_id") for df in results_l]

# Concatenate all results together in a single dataframe:
results_df = pd.concat(results_l, axis=1)

# Rename columns of results_df to match sample_id names in MIA_rnaseq_samples_meta_df:
results_df.columns = [x.split(".")[-1] for x in results_df.columns]

# Important dataframe by sample IDs:
results_df = results_df.sort_index(axis=1)

# Save dataframe to file to be used in edgeR pipeline
# reset index to get gene ids:
results_df.reset_index().to_csv("../{}.csv".format(output_file), index=False)
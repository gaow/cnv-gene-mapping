#!/usr/bin/env python3

import pandas as pd, numpy as np
import os


def main():
    cwd = os.path.expanduser("~/cnv-project/adsp_real_data")
    
    # Load data and do QC
    del_comm_cnv = pd.read_csv(f"{cwd}/del_comm_adsp_hq.gz", compression="gzip", sep="\t")
    del_comm_cnv = del_comm_cnv.drop(columns=["COUNTED","ALT"])
    del_comm_cnv.iloc[:,0] = del_comm_cnv.iloc[:,0].str.replace('chr','')

    n_cnv,n_cols = del_comm_cnv.shape

    percent_na_indv = [del_comm_cnv.iloc[:,i].isna().sum() / n_cnv for i in range(4,n_cols)]
    cutoff_indv = 0.2
    indv_to_keep = [index for index, element in enumerate(percent_na_indv) if element <= cutoff_indv]
    del_comm_cnv = del_comm_cnv.iloc[:,indv_to_keep]

    n_indvs = del_comm_cnv.shape[1] - 4
    percent_na_cnv = [del_comm_cnv.iloc[i,:].isna().sum() / (n_indvs) for i in range(n_cnv)] 
    cutoff_cnv = 0.4
    cnv_to_keep = [index for index, element in enumerate(percent_na_cnv) if element <= cutoff_cnv]
    del_comm_cnv = del_comm_cnv.iloc[cnv_to_keep,:]
    
    # Load reference genes and do global matching
    ref_gene = pd.read_table(f"{cwd}/refGene.clean.gz", compression = "gzip", sep = "\t", header = 0)
    cnv = del_comm_cnv.iloc[:,:4]
    cnv_ref_gene = cnv.merge(ref_gene, how='left', on='CHR')
    match_res = cnv_ref_gene[(cnv_ref_gene['CM'] >= cnv_ref_gene['start']) & (cnv_ref_gene['CM'] <= cnv_ref_gene['end'])
                            |(cnv_ref_gene['POS'] >= cnv_ref_gene['start']) & (cnv_ref_gene['POS'] <= cnv_ref_gene['end'])
                            |(cnv_ref_gene['CM'] <= cnv_ref_gene['start']) & (cnv_ref_gene['POS'] >= cnv_ref_gene['end'])
                            ]
    
    # Seperate reference genes to interrupted and normal
    genes_1 = match_res['gene'].unique()
    genes_0 = ref_gene[~ref_gene['gene'].isin(genes_1)]
    
    # Create a dictionary for building final X matrix dataframe
    gene_dict = {}
    for g in genes_1:
        gene_dict[g] = []
    
    # Iterate each individual and annotate genes
    global_match = match_res[['SNP','gene']]
    for i in range(del_comm_cnv.shape[1] - 4):
        person = del_comm_cnv.iloc[:,[1,(i+4)]]
        union = person.merge(global_match, how='outer', on='SNP')
        valid = union[union['gene'].isna() == False]
        
        # Iterate each interrupted gene in an individual and assign binary values
        for g in genes_1:
            if valid[valid['gene'] == g].iloc[:,1].sum() > 0:
                gene_dict[g].append(1)
            else:
                gene_dict[g].append(0)
    
    # Make a dataframe
    X = pd.DataFrame(gene_dict)
    X.to_csv(f"{cwd}/delX.gz", compression='gzip', sep = "\t", header = True, index = False)
    

if __name__ == "__main__":
    main()

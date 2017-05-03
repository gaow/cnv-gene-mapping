from simulation import *
from scipy import stats
import glob
def get_power(stats_table, causal_genes, p = 0.05):
    '''get power from each simulated dataset.
    First get overlapped genes from stats table and causal genes, and find corresponding p value
    Then power = #(p_value < 0.05) / #(p_value)'''
    overlap_causal_genes = [gene for gene in causal_genes if gene in stats_table["genes"]]
    overlap_p_value = [stats_table["p_value"][stats_table["genes"].index(gene)] 
                       for gene in overlap_causal_genes]
    pvalue_less_than_p = [x for x in overlap_p_value if x < p]
    return len(pvalue_less_than_p)/len(overlap_causal_genes)

def get_typeIerror(stats_table, causal_genes, p = 0.05):
    '''get type I error for each simulated dataset.
    First get the overlapped genes from noncausal genes and genes in stats table, find corresponding p value
    Then type I error = #(p_value < 0.05) / #(p_value)
    '''
    noncausal_overlap_genes = [gene for gene in stats_table["genes"] if gene not in causal_genes]
    noncausal_overlap_genes_pvalue = [stats_table["p_value"][stats_table["genes"].index(gene)]
                                      for gene in noncausal_overlap_genes]
    pvalue_less_than_p = [x for x in noncausal_overlap_genes_pvalue if x < p]
    return len(pvalue_less_than_p)/len(noncausal_overlap_genes_pvalue)

def test_contingency_table(gene_table, method = "Fisher", option = False): 
    if (method == "Fisher"):
        stats_table = [(pvalue(row["n_case_gene"], row["n_ctrl_gene"], row["n_case_nogene"], row["n_ctrl_nogene"]), 
                        row["gene_name"]) for idx, row in gene_table.iterrows()]
        p_value = [x[0].two_tail for x in stats_table]
        genes = [x[1] for x in stats_table]
        stats_table = {"genes": genes, "p_value": "p_value"}
    else:
        table = [( stats.chi2_contingency([[row["n_case_gene"], row["n_ctrl_gene"]], 
                                           [row["n_case_nogene"], row["n_ctrl_nogene"]]], 
                                          correction = option), row["gene_name"] ) 
                 for idx, row in gene_table.iterrows()]
        p_value = [x[0][1] for x in table]
        gene = [x[1] for x in table]
        stats_table = {"p_value": p_value, "genes": gene}
    return stats_table

def get_power_and_typeIerror(input_data, method_option = "chi2", correction_option = False, p_option = 0.05):
    '''use function "load_data" and "get_gene_table" from simulation.py, use simulated dataset as input data,
    and get stats table by using Fisher or chisquare test.
    Then get power and type I error by using functions above, input is stats table and causal genes'''
    sample_table = load_data(input_data)
    causal_genes = sample_table["debug"]["causal genes"]
    gene_table = get_gene_table(sample_table)
    stats_table = test_contingency_table(gene_table, method = method_option, option = correction_option)
    power = get_power(stats_table, causal_genes, p = p_option)
    typeI_error = get_typeIerror(stats_table, causal_genes, p = p_option)
    return {"power": power, "typeI_error": typeI_error, "debug": sample_table["debug"]["args"]}

def run_power_typeIerror(datasets):
    '''input data must be a list of simulated datasets'''
    res = {}
    i = 0
    for data in datasets:
        res_data = get_power_and_typeIerror(data)
        res["dataset_{}".format(i)] = res_data
        i += 1
    return res


from simulation import *
import numpy as np
from pandasql import sqldf

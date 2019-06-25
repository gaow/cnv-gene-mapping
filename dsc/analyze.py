import pandas as pd, numpy as np, seaborn as sns, matplotlib.pyplot as plt
from pandasql import sqldf
from scipy import stats
from fisher import pvalue
import pickle

def save_data(data, filename):
    pickle.dump(data, open(filename, "wb"))

def get_gene_table(gene_df):
    gene_table = {}
    for item in ["case", "ctrl"]:
        gene = pd.concat(gene_df[item])
        query = '''
        SELECT chrom, gene_name, count(gene_name)
        FROM gene
        GROUP BY chrom, gene_name
        ORDER BY count(gene_name) DESC
        '''
        gene_table[item] = sqldf(query)
    gene_table = pd.merge(gene_table["case"], gene_table["ctrl"], how = "outer", on = ["chrom", "gene_name"])
    gene_table["count(gene_name)_x"].fillna(0, inplace=True)
    gene_table["count(gene_name)_y"].fillna(0, inplace=True)
    gene_table = gene_table.rename(columns = {"count(gene_name)_x": "n_case_gene", "count(gene_name)_y": "n_ctrl_gene"})
    n_gene_case = sum(gene_table["n_case_gene"])
    n_gene_ctrl = sum(gene_table["n_ctrl_gene"])
    gene_table["n_case_nogene"] = n_gene_case - gene_table["n_case_gene"]
    gene_table["n_ctrl_nogene"] = n_gene_ctrl - gene_table["n_ctrl_gene"]
    gene_table = gene_table[["gene_name", "n_case_gene", "n_ctrl_gene", "n_case_nogene", "n_ctrl_nogene"]]
    return gene_table

def get_stats(gene_table, num, sort = 0):
    stats_table = [(pvalue(row["n_case_gene"], row["n_ctrl_gene"], row["n_case_nogene"], row["n_ctrl_nogene"]), row["gene_name"]) for idx, row in gene_table.iterrows()]
    p_value = [x[0].two_tail for x in stats_table]
    oddsratio_table = [(stats.fisher_exact([[row["n_case_gene"], row["n_ctrl_gene"]], [row["n_case_nogene"], row["n_ctrl_nogene"]]])[0], row["gene_name"]) for idx, row in gene_table.iterrows()]
    if sort != 0:
        stats_table = sorted(stats_table, reverse = True, key = lambda x: -np.log10(x[0].two_tail))
        oddsratio_table = sorted(oddsratio_table, reverse = True, key = lambda x: x[0] if np.isfinite(x[0]) else -x[0])
    logp_2side = [-np.log10(x[0].two_tail) for x in stats_table]
    top_logp_2side = [-np.log10(x[0].two_tail) for x in stats_table[:num]]
    logp_gene = [x[1] for x in stats_table]
    top_logp_gene = [x[1] for x in stats_table[:num]]
    OR_2side = [x[0] for x in oddsratio_table]
    OR_gene = [x[1] for x in oddsratio_table]
    stats_table = {"p_value": p_value, "top_logp_2side": top_logp_2side, "top_logp_gene": top_logp_gene, "OR_2side": OR_2side, "OR_gene": OR_gene}
    return stats_table

def get_stats_from_input(input_data, num, sort_data = 0):
    '''input data saved from run_simulate step: sample_dup and sample_del separately'''
    # input_data = load_data(input_data)
    sample_gene_table = get_gene_table(input_data)
    sample_stats_table = get_stats(sample_gene_table, num, sort = sort_data)
    return sample_stats_table

def run_stats(input_data, output_data, plt_path, num = 100, sort_option = 0):
    stats_table = get_stats_from_input(input_data, num, sort_data = sort_option)
    stats_table['debug'] = {'simulation_args': input_data['debug']}
    save_data(stats_table, output_data)
    sns.set(style = "whitegrid")
    sns.set_color_codes("pastel")
    f, ax = plt.subplots(figsize = (8, 25))
    p_df = pd.DataFrame({"-logp2side": stats_table["top_logp_2side"], "gene": stats_table["top_logp_gene"]})
    plot = sns.barplot(x = "-logp2side", y = "gene", data = p_df, color = "blue")
    plt.savefig(plt_path)
    return stats_table


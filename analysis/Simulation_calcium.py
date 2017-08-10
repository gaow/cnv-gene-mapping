import random
from random import randint
import pandas as pd
from pandasql import sqldf
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import stats
import math
from sklearn.mixture import GaussianMixture
from datetime import datetime
import pickle
from fisher import pvalue
def load_reference_gene(filename):
    '''Load reference gene database'''
    ref_gene = pd.read_table(filename, compression="gzip", sep="\t", 
                         header = None, usecols=(1,2,4,5,12), 
                         names = ["tx_name", "chrom", "tx_start", "tx_end", "gene_name"])
    return ref_gene.drop_duplicates(subset=("chrom", "tx_start", "tx_end"))

def load_pthwy_gene(filename, n_skiprow = 2):
    '''Load pathway genes. For example, calcium pathway genes.
       The input file must contain a column named "gene_name".
    '''
    pthwy_gene = pd.read_table(filename, skiprows = n_skiprow, header = None, names = ["gene_name"])
    return pthwy_gene

def load_cnv_data(filename):
    '''load cnv data'''
    cnvbed = {}
    dataset = None
    for line in open(filename).readlines():
        if not line.startswith("chr"):
            dataset = line.strip().split()[1].lstrip("name=")
            cnvbed[dataset] = {}
            continue
        line = line.strip().split()
        if not line[0] in cnvbed[dataset]:
            cnvbed[dataset][line[0]] = []
        cnvbed[dataset][line[0]].append((int(line[1]),int(line[2])))

    for dataset in cnvbed.keys():
        for chrom in cnvbed[dataset]:
            cnvbed[dataset][chrom].sort()

    cnvbed_df = {}
    for dataset in cnvbed.keys():
        cnvbed_df[dataset] = {"chrom":[], "cnv_start":[], "cnv_end":[]}
        for chrom in cnvbed[dataset]:
            start, end = tuple(zip(*cnvbed[dataset][chrom])) 
            cnvbed_df[dataset]["chrom"].extend([chrom] * len(start))
            cnvbed_df[dataset]["cnv_start"] += list(start)
            cnvbed_df[dataset]["cnv_end"] += list(end)
        cnvbed_df[dataset] = pd.DataFrame.from_dict(cnvbed_df[dataset]).drop_duplicates(
                                                subset=("chrom", "cnv_start", "cnv_end"))
    return cnvbed_df

def count_cnv_by_block(df, block_size):
    # check how many CNVs start in block_size = 100K blocks genome
    # this cell produces `block_counts`: ((chrom, block_position), count)
    def count_by_blocks(chrom):
        data = df.query('chrom == "{}"'.format(chrom))['cnv_start'].tolist()
        start_pos = min(data)
        end_pos = max(data)
        counts, bins = np.histogram(data, bins = int((end_pos - start_pos) / block_size) + 1, 
                                    range = (start_pos, end_pos))
        return counts, [int(x) for x in bins]
    block_counts = []
    for chrom in set(df['chrom']):
        counts, bins = count_by_blocks(chrom)
        # most blocks contain 0 CNV start. Add 0.5 to each block, so that CNV could start within any block.
        block_counts.extend([((chrom, x), y+0.5) for x, y in zip(bins, counts)])
    return block_counts

def fit_truncated_gaussian_mix(x, k = 10):
    x = x.extend([-i for i in x])
    clf = GaussianMixture(n_components=1, covariance_type='full')
    clf.fit(x)
    return clf

def sample_cnv_length(data, mean_num_cnv):
    return np.random.choice(data, np.random.poisson(mean_num_cnv))

def get_sample_blocks(block_counts, num_cnv):
    '''sample blocks from blocks across genome'''
    probability_distribution = np.array([x[1] for x in block_counts])
    sample_idx = np.random.choice(range(len(block_counts)), num_cnv, 
                                  p = probability_distribution / sum(probability_distribution))
    return sorted([block_counts[idx][0] for idx in sample_idx])

def assign_cnv_to_sample(sample_blocks, sample_len, block_size):
    samples = {'chrom': [], 'cnv_start': [], 'cnv_terminate': []}
    for x, y in zip(sample_blocks, sample_len):
        start_pos = randint(x[1], x[1] + block_size)
        samples['cnv_start'].append(start_pos)
        samples['cnv_terminate'].append(start_pos + int(y))
        samples['chrom'].append(x[0])
    return pd.DataFrame(samples)

def annotate_samples(samples, gene_df):
    query = """
        SELECT cnv.chrom, cnv.cnv_start, cnv.cnv_terminate, gene.tx_name, gene.gene_name
        FROM samples cnv LEFT JOIN gene_df gene
        WHERE cnv.chrom == gene.chrom 
        AND (
        (cnv.cnv_start >= gene.tx_start AND cnv.cnv_start <= gene.tx_end)
        OR
        (cnv.cnv_terminate >= gene.tx_start AND cnv.cnv_terminate <= gene.tx_end)
        OR
        (cnv.cnv_start <= gene.tx_start AND cnv.cnv_terminate >= gene.tx_end)
        )
        """
        # drop_duplicates(): make sure the case that CNV spread multiple txs but each gene to be counted only once
    return sqldf(query).drop_duplicates(subset=("chrom", "cnv_start", "cnv_terminate", "gene_name"))

def get_causal_genes(causal_genes, sample_genes):
    '''get causal genes for each simulated sample'''
    return [x for x in causal_genes if x in sample_genes]

def get_ccnv():
    '''get causal cnvs'''
    return None
    
def p_case(p, num_causal_genes_in_sample, sim_args):
    if num_causal_genes_in_sample == 0:
        return p
    baseline_odds = p / (1 - p)
    if sim_args["odds_ratio_params"]["shape"] is None or sim_args["odds_ratio_params"]["scale"] is None:
        odds_ratio = 1
    else:
        odds_ratio = np.prod([np.random.gamma(sim_args["odds_ratio_params"]['shape'], 
                                              sim_args["odds_ratio_params"]['scale']) 
                              for x in range(num_causal_genes_in_sample)])
    # obtain the power of fisher test by setting odds ratio to 1
    # odds_ratio = 1
    odds = baseline_odds * odds_ratio
    return odds / (1 + odds)

class Environment(dict):
    def __init__(self):
        parameters = {'block_size': 100000,
                      'avg_cnv_per_individual': 5,
                      'n_case': 2000,
                      'n_ctrl': 2000,
                       # set Gamma shape to be 3 instead of 5
                       # 'odds_ratio_params' : None # for H_0
                      'odds_ratio_params': {'shape': 5, 'scale': 1},
                      'prevalence': 0.005,
                      'n_causal_gene': 200,
                      'refgene_file': 'data/refGene.txt.gz',
                      'pthwy_gene_file': '../data/calciumgeneset.txt',
                      'cnv_file': 'data/ISC-r1.CNV.bed',
                      'case_dataset': 'delCases',
                      'ctrl_dataset': 'delControls',
                      'output': 'del_sample'
                     }
        self.update(parameters)
        ## select causal genes randomly, instead of the first 100 in enrichment analysis

def simulate(refgene, cnv_data, args, causal_genes):
    df = cnv_data.drop_duplicates(subset=("chrom", "cnv_start", "cnv_end"))
    block_counts = count_cnv_by_block(df, args['block_size'])
    cnv_length = cnv_data['cnv_end'] - cnv_data['cnv_start']
    status = 1
    case_data = []
    ctrl_data = []
    debug = {'p': [], 'niter': 0, 'time': [str(datetime.now()), None], 'args': dict(args), 'seed': args.seed, 
             'causal genes': causal_genes, 'number of causal genes': [], 'number of genes overlap CNV': [], 
             'simulated CNV length in case': [], 'simulated CNV length in ctrl': []}
    
    while(status):
        sample_len = sample_cnv_length(cnv_length, args['avg_cnv_per_individual'])
        sample_blocks = get_sample_blocks(block_counts, len(sample_len))
        samples = assign_cnv_to_sample(sample_blocks, sample_len, args['block_size'])
        samples = annotate_samples(samples, refgene)
        causal_genes_in_sample = get_causal_genes(causal_genes, samples['gene_name'].tolist())
        p = p_case(args['prevalence'], len(causal_genes_in_sample), args)
        # add the number of causal genes overlapped with simulated CNVs for each simulated sample
        debug['number of causal genes'].append(len(causal_genes_in_sample))
        # add the number of genes overlapped with simulated CNVs, both causal and non-causal genes
        debug['number of genes overlap CNV'].append(len( set(samples['gene_name'].tolist()) ))
        if random.random() < p and len(case_data) < args['n_case']:
            # sample data is a case
            case_data.append(samples)
            debug['p'].append(p)
            debug['simulated CNV length in case'].extend(sample_len)
        if random.random() > p and len(ctrl_data) < args['n_ctrl']:
            # sample data is a control
            ctrl_data.append(samples)
            debug['p'].append(p)
            debug['simulated CNV length in ctrl'].extend(sample_len)
        if len(case_data) == args['n_case'] and len(ctrl_data) == args['n_ctrl']:
            status = 0
        debug['niter'] += 1
    debug['time'][1] = str(datetime.now())
    return {'case': case_data, 'ctrl': ctrl_data, 'debug': debug}

def save_data(data, filename):
    pickle.dump(data, open(filename, "wb"))

def load_data(filename):
    return pickle.load(open(filename, "rb"))


def run_simulation(args, simulation_id = 0):
    np.random.seed(args.seed)
    ref_gene = load_reference_gene(args['refgene_file'])
    cnv_data = load_cnv_data(args['cnv_file'])
    sample_data = simulate(ref_gene, pd.concat([cnv_data[args['case_dataset']], cnv_data[args['ctrl_dataset']]]),
                          args, args.causal_genes)
    save_data(sample_data, '{}_{}.data.pkl'.format(args['output'], simulation_id))
    return sample_data


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
    gene_table = gene_table.rename(columns={"count(gene_name)_x":"n_case_gene", "count(gene_name)_y":"n_ctrl_gene"})
    n_gene_case = sum(gene_table["n_case_gene"])
    n_gene_ctrl = sum(gene_table["n_ctrl_gene"])
    gene_table["n_case_nogene"] = n_gene_case - gene_table["n_case_gene"]
    gene_table["n_ctrl_nogene"] = n_gene_ctrl - gene_table["n_ctrl_gene"]
    gene_table = gene_table[["gene_name", "n_case_gene", "n_ctrl_gene", "n_case_nogene", "n_ctrl_nogene"]]
    return gene_table



def get_stats(gene_table, sort = 0):
    # from website https://pypi.python.org/pypi/fisher/
    stats_table = [(pvalue(row["n_case_gene"], row["n_ctrl_gene"], row["n_case_nogene"], row["n_ctrl_nogene"]), 
                row["gene_name"]) for idx, row in gene_table.iterrows()]
    p_value = [x[0].two_tail for x in stats_table]
    oddsratio_table = [(stats.fisher_exact([[row["n_case_gene"], row["n_ctrl_gene"]], 
                                            [row["n_case_nogene"], row["n_ctrl_nogene"]]])[0], row["gene_name"]) 
                       for idx, row in gene_table.iterrows()]
    if not sort == 0:
        stats_table = sorted(stats_table, reverse=True, key = lambda x: -np.log10(x[0].two_tail))
        oddsratio_table = sorted(oddsratio_table, reverse=True, key=lambda x: x[0] if np.isfinite(x[0]) else -x[0])
    
    logp_2side = [-np.log10(x[0].two_tail) for x in stats_table]
    logp_gene = [x[1] for x in stats_table]
    OR_2side = [x[0] for x in oddsratio_table]
    OR_gene = [x[1] for x in oddsratio_table]
    stats_table = {"p_value": p_value, "logp_2side": logp_2side, "logp_gene": logp_gene, 
                   "OR_2side": OR_2side, "OR_gene": OR_gene}
    return stats_table

def get_stats_from_input(input_data, sort_data = 0):
    '''input data saved from run_simulate step: sample_dup and sample_del separately'''
    input_data = load_data(input_data)
    sample_gene_table = get_gene_table(input_data)
    sample_stats_table = get_stats(del_sample_gene, num=100, sort = sort_data)
    return sample_stats_table

def run_stats(input_data, output_data):
    stats_table = get_stats_from_input(input_data, sort_data=0)
    stats_table['debug'] = {'simulation_args': input_data['args']}
    save_data(stats_table, output_data)
    return stats_table



import plotly

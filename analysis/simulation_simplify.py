
# coding: utf-8

# In[5]:

import random
from random import randint
import pandas as pd
from pandasql import sqldf
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
import math
from collections import Counter
import scipy.special
from scipy.misc import logsumexp
from sklearn.mixture import GaussianMixture
from datetime import datetime
import pickle
from fisher import pvalue

get_ipython().magic('matplotlib inline')


def load_reference_gene(filename):
    '''Load reference gene database'''
    ref_gene = pd.read_table(filename, compression="gzip", sep="\t", 
                         header = None, usecols=(1,2,4,5,12), 
                         names = ["tx_name", "chrom", "tx_start", "tx_end", "gene_name"])
    return ref_gene.drop_duplicates(subset=("chrom", "tx_start", "tx_end"))

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
    #
    block_counts = []
    for chrom in set(df['chrom']):
        counts, bins = count_by_blocks(chrom)
        block_counts.extend([((chrom, x), y) for x, y in zip(bins, counts)])
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
        # drop_duplicates(): make sure the case that CNV spread multiple txs but one gene to be counted only once
    return sqldf(query).drop_duplicates(subset=("chrom", "cnv_start", "cnv_terminate", "gene_name"))

def get_causal_genes(causal_genes, sample_genes):
    '''get causal genes'''
    return [x for x in causal_genes if x in sample_genes]

def get_ccnv():
    '''get causal cnvs'''
    return None
    
def p_case(p, num_causal_genes_in_sample, sim_args):
    if num_causal_genes_in_sample == 0:
        return p
    baseline_odds = p / (1 - p)
    odds_ratio = np.prod([np.random.gamma(sim_args["odds_ratio_params"]['shape'], 
                                          sim_args["odds_ratio_params"]['scale']) 
                          for x in range(num_causal_genes_in_sample)])
    odds = baseline_odds * odds_ratio
    return odds / (1 + odds)

class Environment:
    def __init__(self):
        self.args = {'block_size': 100000,
                     'avg_cnv_per_individual': 5,
                     'n_case': 1,
                     'n_ctrl': 1,
                     # set Gamma shape to be 3 instead of 5
                     'odds_ratio_params': {'shape': 3, 'scale': 1},
                     'prevalence': 0.005
                    }
        
        self.causal_genes = {
            "causal_genes_del": ['RAB2B', 'CHD8', 'TOX4', 'SNORD8', 'METTL3', 'SALL2', 'SNORD9', 'SUPT16H', 
                                 'RPGRIP1', 'MIR3180-3', 'MIR3670-1', 'NOMO2', 'MIR6511A4', 'NPIPA8', 'FSIP2', 
                                 'FSIP2-AS1', 'LOC101927196', 'ATP6V1E1', 'BCL2L13', 'BID', 'CECR1', 'CECR2', 
                                 'CECR3', 'CECR5', 'CECR5-AS1', 'CECR6', 'CECR7', 'FLJ41941', 'GAB4', 'IL17RA', 
                                 'LINC00528', 'LOC100996342', 'LOC100996415', 'LOC101929372', 'LOC105379550', 
                                 'MICAL3', 'MIR3198-1', 'MIR648', 'PEX26', 'SLC25A18', 'TUBA8', 'USP18', 'FAM72C', 
                                 'FAM72D', 'LINC01138', 'NBPF8', 'PPIAL4G', 'MIR6770-2', 'MIR3179-1', 'MIR3180-2', 
                                 'MACROD2', 'FAM189A1', 'LOC100130111', 'MACROD2-AS1', 'LINC00623', 'LINC00869', 
                                 'PPIAL4C', 'MIR3680-2', 'ABCC6P1', 'HNRNPC', 'APBA2', 'C22orf39', 'CDC45', 'CLDN5',
                                 'CLTCL1', 'DGCR14', 'DGCR2', 'GNB1L', 'GOLGA6L7P', 'GOLGA8M', 'GP1BB', 'GSC2', 
                                 'HIRA', 'LINC00895', 'LINC01311', 'LOC100289656', 'MRPL40', 'NPIPB4', 'NSMCE3', 
                                 'PDCD6IPP2', 'SEPT5', 'SEPT5-GP1BB', 'SLC25A1', 'TBX1', 'TSSK2', 'UFD1L', 
                                 'WHAMMP2', 'HSFY1P1', 'PFN1P2', 'XKR3', 'OR4A47', 'C5orf17', 'CNTN4', 'EXOC4', 
                                 'LINC00540', 'LOC101927967', 'LRRC4C', 'SMG1P2', 'SPRY2', 'LOC100288162'],
            "causal_genes_dup": ['CHN2', 'ESYT2', 'KIF26B', 'RPGRIP1', 'C22orf39', 'CDAN1', 'CDC45', 'CLDN5', 
                                 'GNB1L', 'GP1BB', 'HIRA', 'LINC00895', 'MRPL40', 'SEPT5', 'SEPT5-GP1BB', 'STARD9', 
                                 'TBX1', 'TTBK2', 'UFD1L', 'LOC283177', 'AGAP6', 'COL18A1', 'COL18A1-AS1', 
                                 'FAM21EP', 'LOC440910', 'MIR6815', 'POTEE', 'SLC19A1', 'TIMM23B', 'VAV2', 
                                 'WASHC2A', 'PTPRN2', 'ANO2', 'CCNDBP1', 'COLEC12', 'EPB42', 'FAM118A', 'FAM160A1', 
                                 'FBLN1', 'KIAA0930', 'LINC01589', 'LOC100996325', 'LOC728613', 'LRRIQ3', 'MACROD2', 
                                 'MIR1249', 'NUP153', 'NUP50', 'NUP50-AS1', 'PRSS48', 'RAP1GAP2', 'RBM47', 'RIBC2', 
                                 'SH3D19', 'SMC1B', 'TGM5', 'TMEM62', 'UPK3A', 'VWF', 'DGKH', 'KIF13A', 'LINC01266', 
                                 'VWA8', 'VWA8-AS1', 'CRYM-AS1', 'SNX29P1', 'KCNJ12', 'KCNJ18', 'BNC1', 
                                 'LOC105375545', 'MIR128-2', 'AGAP7P', 'BTBD11', 'C3orf67', 'CHD8', 'COL18A1-AS2', 
                                 'FAM110C', 'GMDS', 'HIRIP3', 'INO80E', 'LINC01022', 'MARK3', 'MIR5707', 'MSMB', 
                                 'NCAPG2', 'NCOA4', 'PAK5', 'PLEKHB2', 'PWP1', 'RAB2B', 'RNF103-CHMP3', 'SNORD8', 
                                 'SNORD9', 'SUPT16H', 'TIMM23', 'ZCCHC14', 'C17orf51', 'COX20', 'DLG1', 'EFCAB2']
        }

    def print(self, do_not_print = False):
        if not do_not_print:
            print("Number of causal deletion genes {}".format(len(self.causal_genes)))
        self.a_new_one = 'hello'

env = Environment()
        
def simulate(refgene, cnv_data, args, causal_genes):
    df = cnv_data.drop_duplicates(subset=("chrom", "cnv_start", "cnv_end"))
    block_counts = count_cnv_by_block(df, args['block_size'])
    cnv_length = cnv_data['cnv_end'] - cnv_data['cnv_start']
    status = 1
    case_data = []
    ctrl_data = []
    debug = {'p': [], 'niter': 0, 'time': [str(datetime.now()), None]}
    
    while(status):
        sample_len = sample_cnv_length(cnv_length, args['avg_cnv_per_individual'])
        sample_blocks = get_sample_blocks(block_counts, len(sample_len))
        samples = assign_cnv_to_sample(sample_blocks, sample_len, args['block_size'])
        samples = annotate_samples(samples, refgene)
        causal_genes_in_sample = get_causal_genes(causal_genes, samples['gene_name'].tolist())
        p = p_case(args['prevalence'], len(causal_genes_in_sample), args)
        #debug['p'].append(p)
        if random.random() < p and len(case_data) < args['n_case']:
            # sample data is a case
            case_data.append(samples)
            debug['p'].append(p)
        if random.random() > p and len(ctrl_data) < args['n_ctrl']:
            ctrl_data.append(samples)
            debug['p'].append(p)
        if len(case_data) == args['n_case'] and len(ctrl_data) == args['n_ctrl']:
            status = 0
        debug['niter'] += 1
    debug['time'][1] = str(datetime.now())
    return {'case': case_data, 'ctrl': ctrl_data, 'debug': debug}

def save_data(data, filename):
    pickle.dump(data, open(filename, "wb"))

def load_data(filename):
    return pickle.load(open(filename, "rb"))


# In[6]:

get_ipython().system(' jupyter nbconvert --to script simulation_simplify.ipynb')


# In[ ]:




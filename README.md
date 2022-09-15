# Fine-mapping of genes using Copy Number Variation data
 
 Develop methods that leverage copy number variations (CNVs) for gene mapping in association studies
 
 https://gaow.github.io/cnv-gene-mapping
 
## Write-up documents

 * [Slides](https://www.overleaf.com/8687161qzfgsmmjwvzy) for 03/23/17 CNV collaboration meeting
 * [Manuscript](https://www.overleaf.com/project/5cc3658121e4e25e0c59df63)
 * [Other write-ups](https://github.com/gaow/cnv-gene-mapping/tree/master/writeup)

## Background and literature

* [A background on the SuSiE model, Wang et al 2020](https://github.com/gaow/cnv-gene-mapping/tree/master/reference/20200610_SuSiE.ipynb)
* [A background on VB and MCMC](https://github.com/gaow/cnv-gene-mapping/tree/master/reference/20200623_Approximate_Algorithms.ipynb)
* [Schizophrenia GWAS](https://github.com/gaow/cnv-gene-mapping/tree/master/reference/20200628_SCZ.ipynb)
* [A background on VB-based variable selection model, Carbonetto and Stephens 2012](https://github.com/gaow/cnv-gene-mapping/tree/master/reference/20200629_varbvs.ipynb)
* [An illustration to `pymc3` implementation of Bayesian variable selection via sparse logistic regression](https://github.com/gaow/cnv-gene-mapping/tree/master/reference/20200429_pymc3_example.ipynb)
    - [with a `pymc4` implementation](https://github.com/gaow/cnv-gene-mapping/blob/master/reference/20200519_PyMC4_explore.ipynb)

Please check out the [`reference`](https://github.com/gaow/cnv-gene-mapping/tree/master/reference) folder for more literature references.

## Preliminary analysis

* [Conventional CNV association study in SCZ data](https://github.com/gaow/cnv-gene-mapping/blob/master/analysis/20170216_Enrichment_analysis_of_CNV_in_schizophrenia.ipynb)
    * SCZ data:
        * sample size: ~1300
        * `.bed` file 
        * case-control data
        * no sex chromosomes
    * Reference gene panel
    * Enrichment analysis via _Fisher Exact Test_
* [Conventional CNV association study in simulated data as a diagnosis analysis](https://github.com/gaow/cnv-gene-mapping/blob/master/analysis/20170504_lfdr_pvalue_dist.ipynb)
    * Plot distribution of p-value from simulated data
* [Simulation study scheme, a prototype](https://github.com/gaow/cnv-gene-mapping/blob/master/analysis/20170526_Simulation.ipynb)
* [Assessment of inclusion of covariates (work in progress)](https://github.com/gaow/cnv-gene-mapping/blob/master/analysis/20200531_Confounder_effects.ipynb)

## Method prototypes

* [NUTS sampler for MCMC](https://github.com/gaow/cnv-gene-mapping/blob/master/prototype/Sparse_LR_NUTS.ipynb)
    * Give background of MCMC using spike-and-slab prior
    * Examples of using `pymc3`
* [Single effect model Bayesian logistic regression](https://github.com/gaow/cnv-gene-mapping/blob/master/prototype/LR_single_effect.Rmd)
    * Background of method 'SIER'
    * Code example of self-defined 'SIER'

## Simulation studies

* [Workflow process for running simulation, genome partition and various numerical methods](https://github.com/gaow/cnv-gene-mapping/blob/master/workflow/20190717_workflow.ipynb)
* Summary and visualization of simulation results
    * PIP calibration [analysis](https://github.com/gaow/cnv-gene-mapping/blob/master/analysis/20200519_PIP_Calibration.ipynb) and [results](https://github.com/gaow/cnv-gene-mapping/blob/master/analysis/20200519_PIP_Calibration_display.ipynb)
    * [ROC curve and calibrated mean-pip](https://github.com/gaow/cnv-gene-mapping/blob/master/analysis/20200603_FDR_curve_calibration.ipynb)
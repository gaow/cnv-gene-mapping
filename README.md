# cnv-gene-mapping
 Develop methods that leverage copy number variations (CNVs) for gene mapping
 
 https://gaow.github.io/cnv-gene-mapping
 
 ## Links
 * [Slides](https://www.overleaf.com/8687161qzfgsmmjwvzy) for 03/23/17 CNV collaboration meeting

---
Below are useful notebooks for building software page.

## `analysis`

* `20170216_Enrichment_analysis_of_CNV_in_schizophrenia`
    * Real data:
        * sample size: ~1300
        * `.bed` file 
        * case-control data
        * no sex chromosomes
    * Reference gene panel
    * Enrichment analysis via _Fisher Exact Test_

* `20170504_lfdr_pvalue_dist`
    * Plot distribution of p-value from simulated data

* `20170526_Simulation`
    * Scheme of simulating CNV data manually 

* `20200610_SuSiE`
    * Background of SuSiE

* `20200623_Approximate_Algorithms`
    * Background of MCMC

* `20200628_SCZ`
    * Background of SCZ

* `20200629_varbvs`
    * Background of `varbvs`

## `prototype`

* `Sparse_LR_NUTS`
    * Give background of MCMC using spike-and-slab prior
    * Examples of using `pymc3`

* `LR_single_effect.Rmd`
    * Background of method 'SIER'
    * Code example of self-defined 'SIER'

## `reference` 
* `varbvs-jss.pdf`
    * Introduction to software usage of `varbvs`

* `methodology1.pdf`
    * Assessing risk genes of SCZ using rare number of CNV data

## `workflow`

* `20190717_workflow.ipynb`
    * workflow process for running simulation, genome partition and various numerical methods

* `20200429_pymc3_example.ipynb`
    * General usage and examples of `pymc3`

* `20200518_Manuscripts_Methods.ipynb`
    * Overview of the project's structure including simulation process and numerical methods used

* `20200603_FDR_curve_calibration.ipynb`
    * PIP calibration figures are displayed in `20200519_PIP_Calibration_display.ipynb`
    * plot ROC curve and calibrated mean-pip figure`

* `20200531_Confounder_effects.ipynb`
    * Include gender covariate in the model, but has limited information



# ADSP real data 

**Backgrounds**:

* ADSP: Alzheimer's Disease Sequencing Project

**Logistics**:

Files in `analysis` folder are running analysis for ADSP data in `data` folder and return outputs in `outputs` folder. 

## Descriptions

### `analysis`

* `data_explore`:
    * `cnv_data_report.ipynb`: exploring data quality and how it fits to our method SuSiE-CNV-hybrid. 
    * `get_hist.sh`: Get the histogram for number of genes in a block
    
* `get_genotype_matrix`:
    * `get_X_matrix.ipynb`: Develop a pipeline for intake ADSP CNV data and calculate genotype(X) matrix.
    * `getXMatrix.py`: A complete version of obtaining X matrix to be running on HPC.
    * `Xmatrix(draft).ipynb`: Please ignore.
    
* `workflow.ipynb`: SoS workflow file, which is the same as one `20190717_workflow.ipynb`.

### `data`:

* `all.apoe.pheno.txt`: The file contains sample IDs for each individual and their corresponding phenotypes.
* `del_comm_adsp_hq.gz`: The high-quality ADSP CNV data whose rows are CNVs across genome and columns are CNV information and CNV values of each individual.
* `refGene.clean.gz`: Cleaned hg19 reference gene from GCSC.

### `outputs`:

* `del_X.gz`: The obtained genotype matrix from file `/analysis/get_genotype_matrix/getXMatrix.py`
* `del_y.gz`: The phenotype matrix which consists of two columns `sample.id` and `AD` in `all.apoe.pheno.txt`.
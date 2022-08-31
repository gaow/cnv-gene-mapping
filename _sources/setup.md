# Softwares and packages setup

## [Source](https://github.com/minqiao/cnv-gene-mapping)
```{note}
We only provide installation guide using conda environment to manage packages since it is tested well during whole simulation process.
```

```{warning}
Using `pip` or other package managers have uncertain consequences and errors.
```
## conda installation

### Install conda packages via `.yml` file
```
conda env create -f cnv-env.yml
```
_Note_: the default name of the environment is `cnv`. Users can change it by editing the name.

### Or install conda packages manually
```{note}
* If you have your own conda installed, you can skip this step. 
* We recommned Miniconda due to its "light weight".
```
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
#### Create a conda environment
If you're using _Intel_ chip computer, do
```
conda create -c conda-forge -n cnv pymc3=3.11.2 theano-pymc=1.1.2 mkl=2022.1.0 mkl-service
```
Otherwise, do
```
conda create -c conda-forge -n cnv pymc3 theano-pymc
```

[__PyMC3__](https://docs.pymc.io/en/v3/index.html) is a probabilistic programming package for Python and it allows users to fit Bayesian models using different numerical methods. In our simulation framework, we use Markov chain Monte Carlo (MCMC) method.
#### Install SoS 
```
conda install sos sos-pbs -c conda-forge
conda install sos-notebook jupyterlab-sos sos-papermill -c conda-forge
conda install sos-r sos-python sos-bash -c conda-forge
```
[__SoS__](https://vatlab.github.io/sos-docs/) has full name of "Script of Script" and it is a computational environment for the development and execution of scripts in multiple languages in one Jupyter Notebook:
* `sos` os SoS workflow engine with its command line interface.
* `sos-pbs` is PBS task engine for submitting jobs to Torch, Slurm, IBM LSF etc.
* `sos-notebook` is core SoS notebook module.
* `jupyterlab-sos` is a JupyterLab extension for SoS Notebook.
* `sos-papermill` is a Papermill extension for running SoS notebooks in command line.
* `sos-r`,`sos-python`,`sos-bash` are required language modules for CNV simulation. 

#### Install PLINK
```
conda install -c bioconda plink
```
[__PLINK__](https://www.cog-genomics.org/plink/) is whole genome association analysis toolset developed by Christopher Chang, the Purcell Lab and others for analyzing genotype/phenotype data purely.
(_Note_: The above PLINK downloaded from bioconda channel has version 1.9.)

### Install required R packages
Firstly, we install packages that can be directly installed from conda-forge:
```
conda install -c conda-forge r-devtools 
conda install -c conda-forge r-biocmanager 
conda install -c conda-forge r-data.table
conda install -c conda-forge r-microbenchmark
conda install -c conda-forge r-cowplot
conda install -c conda-forge r-igraph
conda install -c conda-forge r-expm
```
Then, entering R console to install rest of packages:
```
R
BiocManager::install('DNAcopy')
install.packages('saasCNV')
install.packages('varbvs')
install.packages('susieR')
```
The key packages are `varbvs` and `susieR` where they involve comparison with other methods, for example, PyMC3, to conclude power of finding causal genes. Other than those packages, the rest of packages are prerequesites for `susieR`. 

## Docker image installation



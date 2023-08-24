#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=adsp_analysis
#SBATCH --output=adsp_susie_cnv.out
#SBATCH --error=adsp_susie_cnv.err

module load python/cpython-3.8.5
module load R/4.2.0
wkdir="/home/bohjiang/cnv-project/adsp_real_data"
outputs="/home/bohjiang/cnv-project/adsp_real_data/outputs"

cd /home/bohjiang/cnv-project/adsp_real_data

# Whole genome analysis using VARBVS
#sos run workflow.ipynb varbvs_wg \
#    --name varbvs_wg \
#    --cwd $outputs \
#    --genotype-file $wkdir/delX.gz \
#    --phenotype-file $wkdir/dely.gz

# Apply SuSiE-CNV-Hybrid
sos run workflow.ipynb susie_cnv_hybrid \
    --name susie_cnv_hybrid \
    --cwd $outputs \
    --genotype-file $wkdir/delX.gz \
    --phenotype-file $wkdir/dely.gz \
    --hyperparam-file $outputs/delX.varbvs_wg.hyperparams \
    --varbvs-wg-pip $outputs/delX.varbvs_wg.pip

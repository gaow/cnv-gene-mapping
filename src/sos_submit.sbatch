#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1000
#SBATCH --job-name=sos-submitter
#SBATCH --output=sos-log-%J.out
#SBATCH --error=sos-log-%J.err

cd /home/gaow/GIT/github/cnv-gene-mapping
JOB_OPT="-q midway2 -c midway2.yml"

for dataset in simu_dele_30_0528 simu_dele_30_0529; do

sos run workflow/20190717_workflow.ipynb susie \
    --name susie --cwd output/$dataset \
    --genotype-file output/$dataset/deletion.X.colnames_b30.$dataset.X.gz \
    --phenotype-file output/$dataset/deletion.X.colnames_b30.$dataset.y.gz \
    --varbvs_wg_pip output/$dataset/deletion.X.colnames_b30.$dataset.varbvs_wg.pip \
    --hyperparam_file output/$dataset/deletion.X.colnames_b30.$dataset.varbvs_wg.hyperparams \
    -s build --job-size 100 $JOB_OPT &> susie_$dataset.log

sos run workflow/20190717_workflow.ipynb varbvs \
    --name varbvs --cwd output/$dataset \
    --genotype-file output/$dataset/deletion.X.colnames_b30.$dataset.X.gz \
    --phenotype-file output/$dataset/deletion.X.colnames_b30.$dataset.y.gz \
    -s build --job-size 100 $JOB_OPT &> varbvs_$dataset.log
    
sos run workflow/20190717_workflow.ipynb sier \
    --name sier --cwd output/$dataset \
    --genotype-file output/$dataset/deletion.X.colnames_b30.$dataset.X.gz \
    --phenotype-file output/$dataset/deletion.X.colnames_b30.$dataset.y.gz \
    --varbvs_wg_pip output/$dataset/deletion.X.colnames_b30.$dataset.varbvs_wg.pip \
    --hyperparam_file output/$dataset/deletion.X.colnames_b30.$dataset.varbvs_wg.hyperparams \
    -s build --job-size 1 $JOB_OPT &> sier_$dataset.log
    
sos run workflow/20190717_workflow.ipynb hybrid \
    --name hybrid --cwd output/$dataset \
    --genotype-file output/$dataset/deletion.X.colnames_b30.$dataset.X.gz \
    --phenotype-file output/$dataset/deletion.X.colnames_b30.$dataset.y.gz \
    --varbvs_wg_pip output/$dataset/deletion.X.colnames_b30.$dataset.varbvs_wg.pip \
    --hyperparam_file output/$dataset/deletion.X.colnames_b30.$dataset.varbvs_wg.hyperparams \
    -s build --job-size 1 $JOB_OPT &> hybrid_$dataset.log
    
sos run workflow/20190717_workflow.ipynb mcmc \
    --name mcmc --cwd output/$dataset \
    --genotype-file output/$dataset/deletion.X.colnames_b30.$dataset.X.gz \
    --phenotype-file output/$dataset/deletion.X.colnames_b30.$dataset.y.gz \
    --varbvs_wg_pip output/$dataset/deletion.X.colnames_b30.$dataset.varbvs_wg.pip \
    --hyperparam_file output/$dataset/deletion.X.colnames_b30.$dataset.varbvs_wg.hyperparams \
    -s build --job-size 1 $JOB_OPT &> mcmc_$dataset.log

done

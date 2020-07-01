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

for dataset in deldup; do

sos run workflow/20190717_workflow.ipynb mcmc \
    --name mcmc --cwd output/$dataset \
    --genotype-file data/del_dup.X.colnames.gz \
    --phenotype-file data/del_dup.y.no_colrownames.gz \
    --varbvs_wg_pip output/$dataset/del_dup.X.varbvs_wg.pip \
    --hyperparam_file output/$dataset/del_dup.X.varbvs_wg.hyperparams \
    -s build --job-size 10 $JOB_OPT &> mcmc_$dataset.log

sos run workflow/20190717_workflow.ipynb sier \
    --name sier --cwd output/$dataset \
    --genotype-file data/del_dup.X.colnames.gz \
    --phenotype-file data/del_dup.y.no_colrownames.gz \
    --varbvs_wg_pip output/$dataset/del_dup.X.varbvs_wg.pip \
    --hyperparam_file output/$dataset/del_dup.X.varbvs_wg.hyperparams \
    -s build --job-size 6 $JOB_OPT &> sier_$dataset.log
    
sos run workflow/20190717_workflow.ipynb varbvs \
    --name varbvs --cwd output/$dataset \
    --genotype-file data/del_dup.X.colnames.gz \
    --phenotype-file data/del_dup.y.no_colrownames.gz \
    -s build --job-size 30 $JOB_OPT &> varbvs_$dataset.log
    
done
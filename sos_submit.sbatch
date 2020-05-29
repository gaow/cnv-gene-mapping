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

wall_time="2h"
sos run workflow/20190717_workflow.ipynb sier     --name sier --cwd output/simu_dele_30_0528     --genotype-file output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.X_genename.gz     --phenotype-file output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.y.gz     --hyperparam_file output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.X_genename.varbvs_wg.hyperparams     --varbvs_wg_pip output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.X_genename.varbvs_wg.genename.pip --job_size 5 -s build  -q midway2 -c midway2.yml &> sier.log
sos run workflow/20190717_workflow.ipynb hybrid     --name hybrid --cwd output/simu_dele_30_0528     --genotype-file output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.X_genename.gz     --phenotype-file output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.y.gz     --hyperparam_file output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.X_genename.varbvs_wg.hyperparams     --varbvs_wg_pip output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.X_genename.varbvs_wg.genename.pip --mcmc_walltime $wall_time --job_size 5 -s build  -q midway2 -c midway2.yml &> hybrid.log
sos run workflow/20190717_workflow.ipynb mcmc     --name mcmc --cwd output/simu_dele_30_0528     --genotype-file output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.X_genename.gz     --phenotype-file output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.y.gz     --hyperparam_file output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.X_genename.varbvs_wg.hyperparams     --varbvs_wg_pip output/simu_dele_30_0528/deletion.X.colnames_b30.simu_dele_30_0528.samples.X_genename.varbvs_wg.genename.pip --mcmc_walltime $wall_time --job_size 2 -s build  -q midway2 -c midway2.yml &> mcmc.log

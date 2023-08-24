#!/bin/bash

outputs="/home/bohjiang/cnv-project/adsp_real_data/outputs"
wkdir="/home/bohjiang/cnv-project/adsp_real_data"

cd /home/bohjiang/cnv-project/adsp_real_data/analysis

sos run workflow.ipynb get_hist \
    --name real_data_hist \
    --cwd $outputs \
    --genotype-file /home/bohjiang/cnv-project/cnv-gene-mapping-public/data/deletion.X.gz
    #--genotype-file $wkdir/delX.gz 

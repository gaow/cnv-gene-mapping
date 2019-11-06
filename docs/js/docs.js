var analysisArray = ['20171030_Obtain_beta', '20171010_swcnv_TORUS', '20171005_SCZ_CNV', 'ng.3725', '20170906_toy_multiCNVs_multigenes_TORUS', '20170906_Overview_TORUS_Estimate_Alpha', '20170905_Results_toy_example_multi_causal_gene', '20170817_Overview_of_DAP_on_calcium_pathway_simulation', '20170806_Overview_of_dap_on_simulation', '20170801_Overview_TORUS_DAP', '20170726_Dap_TORUS_on_simulation_for_Calcium_Pathway', '20170726_CalciumSignalingPathway', '20170721_swcnv', '20170706_R_varbvs_toy', '20170705_toy_example_multi_outline', '20170627_toy_multigenes', '20170623_toy_example_outline', '20170526_Simulation', '20170525_R_varbvs', '20170504_lfdr_pvalue_dist', '20170504_lfdr_plot', '20170503_Power_typeIerror', '20170425_feather', '20170413_Plot_enrichment_analysis_CNV', '20170216_Enrichment_analysis_of_CNV_in_schizophrenia', 'test.pipeline']
var analysisDict = {'Enrichment-analysis-of-CNVs-in-schizophrenia-1': '20170216_Enrichment_analysis_of_CNV_in_schizophrenia', 'Plots-for-enrichment-analysis-of-CNV-in-Schizophrenia-1': '20170413_Plot_enrichment_analysis_CNV', 'Obtain-.feather-file-from-.pkl-to-analyze-in-R-1': '20170425_feather', 'Obtain-power-and-type-I-error-for-each-simulated-dataset-(.pkl)-1': '20170503_Power_typeIerror', 'Obtain-LFDR-using-.RDS-and-p-value-using-simulated-dataset-(.pkl)-1': '20170504_lfdr_plot', 'Plots-of-p-value-and-LFDR-for-all-genes-in-simulated-dataset-1': '20170504_lfdr_pvalue_dist', 'R:-use-varbvs-for-all-datasets-from-simulation-and-save-as-.RDS-1': '20170525_R_varbvs', 'Simulation-of-Exome-wide-CNV-data-for-case-control-samples-1': '20170526_Simulation', 'Toy-example-of-gene-configuration-overlapped-with-CNV-1': '20170623_toy_example_outline', 'Toy-example-of-multiple-gene-overlapped-with-CNVs-1': '20170627_toy_multigenes', 'Toy-example-of-multi-gene-configuration-overlapped-with-multi-CNV-1': '20170705_toy_example_multi_outline', 'R:-use-package-varbvs-and-feather-to-file-to-obtain-.RDS-file-1': '20170706_R_varbvs_toy', 'Implement-of-table-1-and-Figure-3-in-CNV-in-schizophrenia-in-Sweden-1': '20170721_swcnv', 'Implement-calcium-signaling-pathway-1': '20170726_CalciumSignalingPathway', 'DAP-on-simulation-results-for-calcium-pathway-1': '20170726_Dap_TORUS_on_simulation_for_Calcium_Pathway', 'Overview-of-DAP-and-TORUS-1': '20170801_Overview_TORUS_DAP', 'DAP-on-simulation-results-1': '20170806_Overview_of_dap_on_simulation', 'Overview-of-DAP-on-Calcium-Pathway-Simulation-1': '20170817_Overview_of_DAP_on_calcium_pathway_simulation', 'Results-of-toy-example-for-multi-causal-gene-in-a-region-1': '20170905_Results_toy_example_multi_causal_gene', 'Use-TORUS-to-estimate-for-toy-of-multi-CNV-causal-gene-example-1': '20170906_Overview_TORUS_Estimate_Alpha', 'Toy-example-for-multi-genes-and-multi-CNVs-in-a-region-1': '20170906_toy_multiCNVs_multigenes_TORUS', 'Use-data-from-the-paper-of-Contribution-of-CNV-to-Schizophrenia-(http:-www.nature.com-ng-journal-v49-n1-full-ng.3725.html)-for-analysis-1': '20171005_SCZ_CNV', 'Implement-TORUS-by-using-CNV-data-of-Schizophrenia-in-Sweden-1': '20171010_swcnv_TORUS', 'Integration-of-beta-1': '20171030_Obtain_beta'}
var analysisArrayMap = {'20170216_Enrichment_analysis_of_CNV_in_schizophrenia': 'Enrichment analysis of  ... chizophrenia', '20170413_Plot_enrichment_analysis_CNV': 'Plots for enrichment an ... chizophrenia', '20170425_feather': 'Obtain .feather file fr ... analyze in R', '20170503_Power_typeIerror': 'Obtain power and type I ... taset (.pkl)', '20170504_lfdr_plot': 'Obtain LFDR using .RDS  ... taset (.pkl)', '20170504_lfdr_pvalue_dist': 'Plots of p-value and LF ... ated dataset', '20170525_R_varbvs': 'R: use varbvs for all d ... save as .RDS', '20170526_Simulation': 'Simulation of Exome-wid ... trol samples', '20170623_toy_example_outline': 'Toy example of gene con ... ped with CNV', '20170627_toy_multigenes': 'Toy example of multiple ... ed with CNVs', '20170705_toy_example_multi_outline': 'Toy example of multi-ge ... th multi-CNV', '20170706_R_varbvs_toy': 'R: use package varbvs a ... in .RDS file', '20170721_swcnv': 'Implement of table 1 an ... ia in Sweden', '20170726_CalciumSignalingPathway': 'Implement calcium signaling pathway', '20170726_Dap_TORUS_on_simulation_for_Calcium_Pathway': 'DAP on simulation resul ... cium pathway', '20170801_Overview_TORUS_DAP': 'Overview of DAP and TORUS', '20170806_Overview_of_dap_on_simulation': 'DAP on simulation results', '20170817_Overview_of_DAP_on_calcium_pathway_simulation': 'Overview of DAP on Calc ... y Simulation', '20170905_Results_toy_example_multi_causal_gene': 'Results of toy example  ... in a region', '20170906_Overview_TORUS_Estimate_Alpha': 'Use TORUS to estimate α ... gene example', '20170906_toy_multiCNVs_multigenes_TORUS': 'Toy example for multi-g ... in a region', '20171005_SCZ_CNV': 'Use data from the paper ... for analysis', '20171010_swcnv_TORUS': 'Implement TORUS by usin ... ia in Sweden', '20171030_Obtain_beta': 'Integration of $\\\\beta$'}
var prototypeArray = ['toy_example_test', '20170822_Fisher_test_by_block', '20170820_Pickle_Block_Matrix', '20170710_dap_on_toy', '20170606_toy_4gene', '20170605_toy_3gene', '20170504_R_useRDS', '20170113_Toy_from_saasCNV', '20160817_VCF_Toy']
var prototypeDict = {'Exploratory-of-a-toy-CNV-dataset-in-VCF-format-1': '20160817_VCF_Toy', 'Visualize-genome-wide-SCNA-Profile-1': '20170113_Toy_from_saasCNV', 'R:-read-.RDS-file-and-draw-histogram-of-PIP-1': '20170504_R_useRDS', 'Toy-example-of-3-gene-configuration-overlapped-with-one-CNV-1': '20170605_toy_3gene', 'Toy-example-of-4-gene-configuration-overlapped-with-one-CNV-1': '20170606_toy_4gene', 'Analysis-of-toy-data-via-DAP-1': '20170710_dap_on_toy', 'Transfer-old-version-pickle-to-new-version-pickle-in-debug-1': '20170820_Pickle_Block_Matrix', 'import-pandas-as-pd-1': '20170822_Fisher_test_by_block', 'from-simulation-import-1': 'toy_example_test'}
var prototypeArrayMap = {'20160817_VCF_Toy': 'Exploratory of a toy CN ... n VCF format', '20170113_Toy_from_saasCNV': 'Visualize genome-wide SCNA Profile', '20170504_R_useRDS': 'R: read .RDS file and d ... ram of "PIP"', '20170605_toy_3gene': 'Toy example of 3-gene c ... with one CNV', '20170606_toy_4gene': 'Toy example of 4-gene c ... with one CNV', '20170710_dap_on_toy': 'Analysis of toy data via DAP', '20170820_Pickle_Block_Matrix': 'Transfer old version pi ... e in "debug"', '20170822_Fisher_test_by_block': 'import pandas as pd', 'toy_example_test': 'from simulation import *'}
var dscArray = ['20190627_Clean_RefGene', '20190628_Xmatrix', '20190628_Xmatrix_python_kernel', '20190710_dele_Genome_window', '20190712_dup_genome_window', '20190716_minimal_working_example_gene_block_matrix', '20190717_workflow', '20190717_workflow_R_kernel_test', '20190717_workflow_python_kernel_test', '20190819_workflow_report', '20190827_varbvs_whole_genome', '20190828_R_test', '20190908_proposal_manuscript', '20190909_Boom_Spike_Slab', '20190917_logitBvs', '20190918_logitBvs_report', '20191005_PyMC3_SpikeSlab', '20191008_PyMC3_report', '20191009_varbvs_gene_block', '20191012_ASHG_poster', '20191022_logistic_single_effect', '20191023_logistic_report', '20191031_Data_Incubator', '20191102_Bike', '20191103']
var dscDict = {'Obtain-clean-reference-gene-1': '20190627_Clean_RefGene', 'Obtain-X-matrix-1': '20190628_Xmatrix', 'import-pandas-as-pd,-numpy-as-np-1': '20191031_Data_Incubator', 'Obtain-Genome-Window-(40-100-genes-per-window)-1': '20190710_dele_Genome_window', 'Copy-model-simulation-and-analysis-workflow-1': '20190717_workflow', 'library(susieR)-1': '20190717_workflow_R_kernel_test', 'Simulation,-analyze-and-plots-1': '20190819_workflow_report', 'Simulation-and-analysis-result-overview-1': '20190827_varbvs_whole_genome', 'Evaluate-copy-model-simulation-1': '20190828_R_test', 'Integrated-fine-mapping-of-non-coding-disease-variants-with-functional-genomics-data-1': '20190908_proposal_manuscript', 'Use-BoomSpikeSlab-1': '20190909_Boom_Spike_Slab', 'library(pogit)-1': '20190917_logitBvs', 'The-results-of-varbvs-result-and-logitBvs-1': '20190918_logitBvs_report', 'Logistic-regression-spike-slab-prior-using-PyMC3-1': '20191005_PyMC3_SpikeSlab', 'PyMC3-report:-Logistic-regression-spike-slab-prior-1': '20191008_PyMC3_report', 'Use-varbvs-result-to-compare-1': '20191009_varbvs_gene_block', 'Introduction-(Use-dot-for-numbering)-1': '20191012_ASHG_poster', 'Single-effect-using-logistic-1': '20191022_logistic_single_effect', 'Single-effect-CNV-gene-pattern-report-using-new-logistic-and-other-methods-1': '20191023_logistic_report', 'import-pandas-as-pd,-numpy-as-np,-matplotlib.pyplot-as-plt-1': '20191103'}
var dscArrayMap = {'20190627_Clean_RefGene': 'Obtain clean reference gene', '20190628_Xmatrix': 'Obtain X matrix', '20190628_Xmatrix_python_kernel': 'import pandas as pd, numpy as np', '20190710_dele_Genome_window': 'Obtain Genome Window (4 ... per window)', '20190712_dup_genome_window': 'import pandas as pd, numpy as np', '20190716_minimal_working_example_gene_block_matrix': 'import pandas as pd, numpy as np', '20190717_workflow': 'Copy model simulation a ... sis workflow', '20190717_workflow_R_kernel_test': 'library(susieR)', '20190717_workflow_python_kernel_test': 'import pandas as pd, numpy as np', '20190819_workflow_report': 'Simulation, analyze and plots', '20190827_varbvs_whole_genome': 'Simulation and analysis ... ult overview', '20190828_R_test': 'Evaluate copy model simulation', '20190908_proposal_manuscript': 'Integrated fine-mapping ... enomics data', '20190909_Boom_Spike_Slab': 'Use BoomSpikeSlab', '20190917_logitBvs': 'library(pogit)', '20190918_logitBvs_report': 'The results of varbvs r ... and logitBvs', '20191005_PyMC3_SpikeSlab': 'Logistic regression spi ... using PyMC3', '20191008_PyMC3_report': 'PyMC3 report: Logistic  ... e slab prior', '20191009_varbvs_gene_block': 'Use varbvs result to compare', '20191012_ASHG_poster': 'Introduction (Use dot for numbering)', '20191022_logistic_single_effect': 'Single effect using logistic', '20191023_logistic_report': 'Single effect CNV-gene  ... ther methods', '20191031_Data_Incubator': 'import pandas as pd, numpy as np', '20191102_Bike': 'import pandas as pd, nu ... yplot as plt', '20191103': 'import pandas as pd, nu ... yplot as plt'}
var setupArray = ['plink', 'cnv.', 'index']
var setupDict = {'plink-1': 'plink', 'install-saasCNV-1': 'setup'}
var setupArrayMap = {'plink': 'plink', 'setup': 'install saasCNV'}
var writeupArray = ['20190603_ASHG_Abstract', 'CNV_meeting_2016', 'CNV_meeting_2017', 'Dataset_List', 'Glossary', 'Learning_Resources', 'Literature_Notes', 'Useful_Info', 'dbGaP_Proposal']
var writeupDict = {'CNV-Abstract-1': '20190603_ASHG_Abstract', 'Project-meeting-20161208-1': 'CNV_meeting_2016', 'CNV-meeting-notes-starting-from-June-2017-1': 'CNV_meeting_2017', 'CNV-datasets-1': 'Dataset_List', 'Glossary-1': 'Glossary', 'Learning-resources-1': 'Learning_Resources', 'Literature-notes-1': 'Literature_Notes', 'Useful-Info-1': 'Useful_Info', 'The-application-below-follows-from-templates-on-this-site-(https:-github.com-aaronwolen-gibbs-brain-project-wiki-Other-dbGap-applications).-1': 'dbGaP_Proposal'}
var writeupArrayMap = {'20190603_ASHG_Abstract': 'CNV Abstract', 'CNV_meeting_2016': 'Project meeting 20161208', 'CNV_meeting_2017': 'CNV meeting notes start ... om June 2017', 'Dataset_List': 'CNV datasets', 'Glossary': 'Glossary', 'Learning_Resources': 'Learning resources', 'Literature_Notes': 'Literature notes', 'Useful_Info': 'Useful Info', 'dbGaP_Proposal': 'The application below f ... plications).'}
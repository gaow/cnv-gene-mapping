# pipeline variables
# $block_size: number of SNPs in a block
# $avg_cnv_per_individual: number of CNV in each individual
# $n_case: number of cases
# $n_ctrl: number of controls
# $odds_ratio_params: odds ratio parameters, Gamma distribution, shape and scale
# $prevalence: ASD prevalence
# $n_causal_gene: number of causal genes

# simulation steps
simulation: simulation_functions.py + Python(data = run_simulation(seed, ref_gene_fn, pathway_gene_fn, CNV_fn, indel, block_size, prevalence, avg_cnv_per_individual, odds_ratio_params, n_case, n_ctrl, id)))
    seed: 999
    ref_gene_fn: ${ref_gene}
    pathway_gene_fn: ${pathway_gene}
    CNV_fn: "../data/ISC-r1.CNV.bed"
    indel: "del"
    block_size: 20000
    prevalence: 0.005
    avg_cnv_per_individual: 5
    odds_ratio_params: {'shape': 5, 'scale': 1}
    n_case: 2000
    n_ctrl: 2000
    id: 0
    $simu_res: data

# analyze
fisher_test: analyze.py + Python(res = run_stats(x, plt_path, num, sort_option))
    num: 100
    sort_option: 1
    plt_path: file(png)
    # plt_path: None
    x: $simu_res
    $stats: res

# score
# score: Python(...)

DSC:
    define:
        simulate: simulation
        analyze: fisher_test
#        score: 
    run: simulate * analyze
    global:
        ref_gene: "../data/refGene.txt.gz"
	pathway_gene: "../data/calciumgeneset.txt"

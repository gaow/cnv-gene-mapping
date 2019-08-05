simulate: simulation_y.py + Python(x1, x2 = simulate_y(genotype_file, prevalence, shape, scale, ctrl_case_ratio, seed))
    genotype_file: "deletion_geneblock.sample.1.gz"
    prevalence: 0.05
    shape: 3
    scale: 1
    ctrl_case_ratio: 1
    seed: 999
    $case: x1
    $ctrl: x2

# analyze
fisher_test: analyze.py + Python(res = run_stats(x, plt_path, num, sort_option))
    num: 100
    sort_option: 1
    plt_path: file(png)
    # plt_path: None
    x: $simu_res
    $stats: res

susie: R(y = c(rep(1,nrow(case), rep(0, nrow(ctrl)))); 
	 res = susieR::susie(rbind(case,ctrl), y, L = L, scaled_prior_variance = pve))
    case: $case
    ctrl: $ctrl
    L: 10
    pve: 0.005
    $result: res

DSC:
    define:
        analyze: susie
    run: simulate * analyze

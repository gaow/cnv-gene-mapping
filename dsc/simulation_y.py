def simulate_y(genotype_file, prevalence, shape, scale, ctrl_case_ratio = 1, seed = 999):
  import pandas as pd, numpy as np
  np.random.seed(seed)
  data = pd.read_csv(genotype_file, compression = "gzip", sep = "\t", header = None)
  beta0 = np.log(prevalence/(1-prevalence))
  beta1s = np.log(np.random.gamma(shape, scale, data.shape[1])) # ORs follow gamma(5,1)
  logit_y = np.matmul(data.values, beta1s) + beta0
  ys_p = np.exp(logit_y) / (1+np.exp(logit_y))
  ys = np.random.binomial(1, ys_p)
  case_index = np.ravel(np.where(ys == 1))
  print(case_index[:10])
  ctrl_index = sorted(np.random.choice(np.ravel(np.where(ys == 0)), int(len(case_index) * ctrl_case_ratio)))
  genotype_case = data.iloc[case_index,:]
  genotype_ctrl = data.iloc[ctrl_index,:]
  return genotype_case, genotype_ctrl

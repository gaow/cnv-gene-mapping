# Overview

## CNV
Copy Number Variation (CNV) is a structural variation, for example, large genomic insertion or deletion events, which are a type of structural variation of an organism’s chromosome. The length of CNVs varies to a great extent, often spread from 50 base pairs to kilo- or even mega-bases.

## Goal
The goal of this project is to develop a statistical framework that leverage genome-wide CNVs for mapping susceptibility genes in schizophrenia (SCZ).

## Challenge
The challenge is that CNVs often span multiple genes, and it is difficult to distinguish susceptible or causal gene(s) from other genes in the same CNV event.

## Idea
We developed a new approach that exploits large-scale genome-wide CNV data in case-control studies to map susceptibility genes. It is inspired by statistical fine-mapping of causal variants in linkage-disequilibrium blocks from GWAS. Unlike existing approaches that directly test for association signals between CNV and disorders, our method seeks to identify true susceptibility genes in CNV events in a rigorous statistical framework. Genome-wide CNV events are first divided into disjoint CNV-gene genomic regions or blocks to ensure no CNV events span more than one block, i.e. no CNVs in common between different regions. For genes located in a certain block, we test for their associations with SCZ while accounting for correlations between genes induced by CNV in the same block. We accomplish this by 3 existing and 2 newly developed approaches as follows.

First we use 3 existing methods as follows.

1). Variational Bayesian variable selection (varbvs) [varbvs]: a Bayesian variable selection methods for the analysis of large-scale data sets built on Bayesian models for variable selection in regression and variational approximation techniques.

2). Sum of Single Effect (SuSiE): a recently developed Bayesian variable selection method, which select a small number of putative risk genes among multiple correlated genes that best explain the data. SuSiE estimates posterior inclusion probabilities (PIPs) of all putative risk genes as well as 95% credible sets (i.e. the set of genes that cover all risk genes with high probability). 

3). Fisher's p-value

Then we compare the results with 2 newly developed approaches as follows.

1). Markov Chain Monte Carlo (MCMC): We use Monte Carlo strategy to sample from the posterior distribution. 

2). Single effect logistic model


## Table of Contents
```{tableofcontents}
```

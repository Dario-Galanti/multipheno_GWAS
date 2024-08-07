# gemmaGWAS
Scripts for running multi-phenotype GWAS using mixed models with [GEMMA](https://github.com/genetics-statistics/GEMMA).<br/>
Additional functionality includes pleasant visualization, optional phenotype transformations and enrichment of variants neighbouring a-priori candidate genes.


SCRIPTS DESCRIPTION: <br/>

[gemma_input.sh](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/gemmaGWAS/gemma_input.sh)<br/>
Prepare PLINK input files for GEMMA.

[gemmaGWA_multipleno.sh](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/gemmaGWAS/gemmaGWA_multipleno.sh)<br/>
Running in parallel for each phenotype (i) [GEMMA](https://github.com/genetics-statistics/GEMMA), (ii) manhattan and qqplot visualization, (iii) enrichment of a-priori candidates as in [Atwell et al. 2010](https://www.nature.com/articles/nature08800) and (iv) [enrichment_plot.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/gemmaGWAS/enrichment_plot.R) for visualization of the enrichment.

[qqman_GEMMA.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/gemmaGWAS/qqman_GEMMA.R)<br/>
Making manhattan and qqplots from [GEMMA](https://github.com/genetics-statistics/GEMMA) results.

[enrichment_plot.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/gemmaGWAS/enrichment_plot.R)<br/>
Plotting the enrichment of a-priori candidates for increasing -log(p) thresholds.

[enrich_sig_genome_rotation.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/gemmaGWAS/enrich_sig_genome_rotation.R)<br/>
Calculating significance of the enrichment of variants closeby a priori candidate genes for a single -log(p) threshold (usually bonferroni). It implements a genome rotation scheme described at the beginning of the script or for example by [Brachi et al. 2010](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000940).


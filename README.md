# multipheno_GWAS
Scripts for running multi-phenotype GWAS using either mixed models with the R package "rrBLUP", or alternatively GEMMA in [gemmaGWAS](https://github.com/Dario-Galanti/multipheno_GWAS/tree/main/gemmaGWAS).<br/>
Additional functionality includes pleasant visualization, optional phenotype transformations and enrichment of variants neighbouring a-priori candidate genes.


SCRIPTS DESCRIPTION: <br/>
[vcfGWAformat_IBSmx.sh](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/vcfGWAformat_IBSmx.sh)<br/>
This script takes an imputed vcf file and does all the following processes which are useful before running GWAS:<br/>
- It uses [PLINK](https://www.cog-genomics.org/plink/) to prune variants in strong LD and produce an Isolation By State matrix.<br/>
- It applies a stricter MAF filtering and converts the vcf to PLINK format for GWAS (this part is reference genome specific!!!)<br/>
- It calculates the number of indipendent GWAS tests after LD pruning according to [Sobota et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4334751/).
- Optionally, it outputs a zygosity report for all individuals.

[multiphenoGWAS.sh](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/multiphenoGWAS.sh)<br/>
Wrapper script running [gwas_rr_kinship_multipheno.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/gwas_rr_kinship_multipheno.R) for multi-phenotypes GWAS, enrichment of a-priori candidates based on the method established by [Atwell et al. 2010](https://www.nature.com/articles/nature08800) and [enrichment_plot.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/enrichment_plot.R) for visualization of the enrichment.

[gwas_rr_kinship_multipheno.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/gwas_rr_kinship_multipheno.R)<br/>
Script for running multi-phenotype GWAS using mixed models with the R package ["rrBLUP"](https://cran.r-project.org/web/packages/rrBLUP/index.html). Transformation can be applied to phenotypes deviating heavily from normality and visually nice manhattan and qqplots are produced using the R package ["qqman"](https://cran.r-project.org/web/packages/qqman/index.html). NAs in the multi-phenotypes file are allowed.

[enrichment_plot.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/enrichment_plot.R)<br/>
Script for plotting the enrichment of a-priori candidates.

PUBLICATIONS: <br/>
[Genetic and environmental drivers of large-scale epigenetic variation in Thlaspi arvense](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010452)

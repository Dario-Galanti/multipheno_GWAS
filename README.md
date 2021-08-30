# multipheno_GWAS
Scripts for running multi-phenotype GWAS using mixed models with the R package "rrBLUP".<br/>
Additional functionality includes high quality visualization, phenotype transformations and enrichment of variants neighbouring a-priori candidate genes.


SCRIPTS DESCRIPTION: <br/>
[vcfGWAformat_IBSmx.sh](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/vcfGWAformat_IBSmx.sh)<br/>
This script takes an imputed vcf file and does all the following processes which are useful before running GWAS:<br/>
- It uses [PLINK](https://www.cog-genomics.org/plink/) to prune variants in strong LD and produce an Isolation By State matrix.<br/>
- It applies a stricter MAF filtering and converts the vcf to PLINK format for GWAS (this part is reference genome specific!!!)<br/>
- It calculates the number of indipendent GWAS tests after LD pruning according to [Sobota et al. 2015[(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4334751/).
- Optionally, it outputs a zygosity report for all individuals.

[multiphenoGWAS.sh](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/multiphenoGWAS.sh)<br/>
Wrapper script running [gwas_rr_kinship_multipheno.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/sderio_gwas_multipheno_BinAC.sh) for multi-phenotypes GWAS, enrichment of a-priori candidates based on the method established by [Atwell et al. 2010](https://www.nature.com/articles/nature08800) and [enrichment_plot.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/enrichment_plot.R) for visualization of the enrichment.

[gwas_rr_kinship_multipheno.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/gwas_rr_kinship_multipheno.R)<br/>
Script for running multi-phenotype GWAS using mixed models with the R package ["rrBLUP"](https://cran.r-project.org/web/packages/rrBLUP/index.html). Transformation can be applied to phenotypes deviating heavily from normality and visually nice manhattan and qqplots are produced using the R package ["qqman"](https://cran.r-project.org/web/packages/qqman/index.html). NAs in the multi-phenotypes file are allowed.

[enrichment_plot.R](https://github.com/Dario-Galanti/multipheno_GWAS/blob/main/enrichment_plot.R)<br/>
Script for visualization of the enrichment of a-priori candidates.

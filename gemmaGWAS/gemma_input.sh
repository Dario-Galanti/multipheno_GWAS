## Author: Dario Galanti (July 2022)
## Aim: Prepare input files for running GEMMA
## Documentation: https://www.xzlab.org/software/GEMMAmanual.pdf
## Tutorial: http://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf

##NB: GEMMA runs with two kinds of input files:
## - PLINK -> Easy to obtain and more symplistic
## - BIMBAM -> It stores dosage information or posterior mean genotypes (Is a value
## between 0 to 2 that can be interpreted as the minor allele dosage: 0 is homozygous for the
## major allele, 1 is a heterozygote, and 2 is a homozygote for the minor allele.)
## For conversion check -> https://github.com/genetics-statistics/GEMMA/issues/103

## DEFINE SOME TOOLS
Rscript=~/miniconda3/envs/R/bin/Rscript
vcftools=~/miniconda3/envs/vcftools/bin/vcftools
bgzip=~/miniconda3/envs/bwa/bin/bgzip
plink=~/miniconda3/envs/vcftools/bin/plink

## DEFINE INPUT FILES
#geno_file=input_v4v5lift/Ta_v5_vrts_10NA_biall_4MAF_imp_withUSspls_withref.raw		# 4MAF filt after imputation (1856975). More accurate!
#geno_file=input_v4v5lift/Ta_v5_vrts_10NA_biall_4MAF_imp_GWAspls_noCommRegs.raw		# 4MAF filt after imputation (141552). More accurate!
#geno_file=input_v4v5lift/Ta_v5_vrts_10NA_biall_4MAF_imp_GWAspls_withref.raw		# 4MAF filt after imputation (1422472). More accurate!
input=/beegfs/work/bbmdg01/Ta_v4v5_lifted/withUSspls/Ta_v5_vrts_10NA_biall_1MAF_imp_withUSspls_withref
output=/beegfs/work/bbmdg01/GEMMA/plink_SNPs/Ta_v5_vrts_10NA_biall_1MAF_imp_withUSspls_withref
kinship=input_v4v5lift/Ta_v5_vrts_10NA_biall_1MAF_imp_8prun_withUSspls_withref.kinship	# OPTIONAL PARAM. K estimated from 1MAF, 0.8 pruned SNPs (more accurate!)

mkdir -p $outdir

## PLINK SNP INPUT
$plink --file $input --allow-extra-chr --make-bed --out $output


## PREPARE PHENO FILE
## We need to paste a phenotype on the 6th column of the plink .fam file
template=Ta_v5_vrts_10NA_biall_4MAF_imp_withUSspls_withref.fam
pheno_input=SupFamTIPs_withUS_GWA_log_1out.txt
col=9
pheno_file=SupFamTIPs_withUS_GWA_log_1out_RIR.fam
rm $pheno_file
re='^[0-9]+'	# To test if my variable is a number
while read l; do
 spl=$(echo $l | cut -d" " -f1)
 val=$(grep -P $spl"\t" $pheno_input | cut -f $col)
 if ! [[ $val =~ $re ]];then val="NA";fi
 newline=$(echo $l | cut -d" " -f1-5)" "$val
 echo $newline >> $pheno_file
done < $template



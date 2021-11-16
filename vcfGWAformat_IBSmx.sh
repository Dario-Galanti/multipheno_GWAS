#PBS -l nodes=1:ppn=8 #Nodes and cores
#PBS -l walltime=20:00:00
#PBS -l mem=20gb 
#PBS -S /bin/bash
#PBS -N GWASformat #Job name
#PBS -j oe
#PBS -q short
#PBS -o /beegfs/work/bbmdg01/logs/GWASformat_IBDmx.out
#PBS -e /beegfs/work/bbmdg01/logs/GWASformat_IBDmx.err

## IMPORTANT: Strict MAF filtering (0.05) should be applied AFTER imputation with BEAGLE to increase imputation accuracy as found in the paper below!
## "Impact of pre-imputation SNP-filtering on genotype imputation results"

## IMPORTANT: Define number of independent tests in GWAS: this is done by SNP pruning with Plink according to the study
## "Addressing population-specific multiple testing burdens in genetic association studies"
## Independent tests = Variants passing 0.1NAs, 0.05MAF, biallelic, pruning (--indep-pairwise 100 5 0.3)

### Aim: This script takes an imputed vcf file with a mild MAF filtering and performs the following analysis and steps:
###			1) From 1MAF_imputed input -> LD pruning (--indep-pairwise 50 5 0.8) IBS matrix extraction and heatmap + dendrogram plotting
###			2) From 1MAF_imputed input -> strict 5MAF filtering and formatting for GWAS
###			3) From 5MAF_imputed -> LD pruning to calculate number of independent tests (--indep-pairwise 100 5 0.3) see note above
###			4) From 1MAF_imputed input -> strict 5MAF filtering, LD pruning (--indep-pairwise 50 5 0.95) and formatting for GWAS with DMRs (to reduce computational load)
### Author: Dario Galanti Mar 2021 (updated)
### Input: imputed vcf file with mild 1MAF filtering (0.01)
### Run: qsub -q short vcfGWAformat_IBSmx.sh

## Define tools
work="/beegfs/work/bbmdg01"
java=~/miniconda3/envs/vcftools/bin/java
beaglejar=~/miniconda3/envs/vcftools/share/beagle-5.1_24Aug19.3e8-0/beagle.jar
bgzip=~/miniconda3/envs/bwa/bin/bgzip
plink=~/miniconda3/envs/vcftools/bin/plink
Rscript=~/miniconda3/envs/R/bin/Rscript
IBSplot_script=${work}/ibs_plot.R

## IMPORTANT: Define parameters for strict MAF filtering and pruning
MAF=0.04				#DEFINE strict Minor Allele Frequency for GWAS
pruneLD_IBS=0.8			#DEFINE LD for pruning SNP to output the IBS matrix
pruneLD_dmrGWA=0.95		#DEFINE LD for pruning SNPs to use for DMR GWAS

## Define input and outputs
wDir=${work}/Ta_v3_10NA_1MAF
#cd $wDir
fin=${wDir}/Ta_v3_vrts_1MAF_imputed_GWAspls_withref.vcf.gz
fin_base=${wDir}/Ta_v3_vrts_1MAF_imputed_GWAspls_withref
ibs_base=${wDir}/Ta_v3_vrts_1MAF_imp_8prun_GWAspls_withref
gwa_base=${wDir}/Ta_v3_vrts_4MAF_imp_GWAspls_withref
ind_tests=${wDir}/4MAF_3prun_ind_tests
dmrGWA_base=${wDir}/Ta_v3_vrts_4MAF_imp_98prun_GWAspls_withref
zigosity_report=${wDir}/zigosity_report_1MAF_imp_8prun_v3.txt
filter_report=${wDir}/Report_vcfGWAformat_v4.txt

## OPTIONAL: ADD REFERENCE GENOTYPE TO THE VCF FILE
#Ref=MN106A
#zcat $fout | awk -v Ref="$Ref" '{if(substr($1,1,2)=="##"){print}else if(substr($1,1,2)=="#C"){OFS="\t";print $0,Ref}else{OFS="\t";print $0,"0|0"}}' | $bgzip -c > $vcf_withref
## OPTIONAL: INDEX VCF FILE
#~/miniconda3/envs/bwa/bin/tabix -p vcf $fin

## IMPORTANT!!! NB: FORMATTING FOR GWAS REQUIRES SOME RELEVANT FORMATTING DESCRIBED BELOW.
## 1) Substitute "_" with "-" in sample names or plink crashes
## 2) Add variant name -> "chr_pos"
## 3) IMPORTANT!!! Plink is very very bad at handling non-model species with Scaffold-level assemblies.
## 3) We use --allow-extra-chr, but this causes eg. Chr1 --> 1 but Scaffold_1 --> Scaffold_1.
## 3) SOLUTION: AFTER running plink we change all scaffold names in the map file with the following rationale:
## 3) SOLUTION: Scaffold_n --> n+7 (eg. Scaffold_1 --> 8) in order not to re-use Chr numbers

## Substitute "_" with "-" in sample names or plink crashes
zcat $fin | awk 'OFS="\t"{if(!/\#/){$3=$1"_"$2;print} else {gsub("_","-");print}}' > ${wDir}/temp_ready.vcf
$plink --vcf ${wDir}/temp_ready.vcf --allow-extra-chr --recode --out $fin_base

### 1a) IBS matrix LD pruning (50 SNPs window, remove a SNP if LD>0.8, shift window by 10SNPs and repeat process)
$plink --file $fin_base --allow-extra-chr --indep-pairwise 50 5 $pruneLD_IBS --out tmp1
$plink --file $fin_base --allow-extra-chr --extract tmp1.prune.in --recode --out $ibs_base
rm tmp1.prune.in tmp1.prune.out

### 1b) Extract IBS matrix from vcf and plot matrix heatmap and dendrogram
$plink --file $ibs_base --allow-extra-chr --distance square ibs flat-missing --out $ibs_base
$Rscript --vanilla $IBSplot_script $ibs_base.mibs $ibs_base.mibs.id

### 2) Strict MAF filtering (0.05), format for GWAS
$plink --file $fin_base --allow-extra-chr --maf $MAF --recode --out $gwa_base
$plink --file $gwa_base --allow-extra-chr --recode A --out $gwa_base			# Make raw file for GWAS
awk 'OFS="\t"{if($1~/^Scaffold/){gsub(/Scaffold_/,"",$1);$1=($1+7)};print}' ${gwa_base}.map > ${gwa_base}_chrint.map

### 3) LD pruning to calculate number of independent tests (--indep-pairwise 100 5 0.3)
$plink --file $gwa_base --allow-extra-chr --indep-pairwise 100 5 0.3 --out $ind_tests
rm ${ind_tests}.prune.out

### 4) Strict MAF filtering (0.05), format for GWAS with DMRs
$plink --file $fin_base --allow-extra-chr --maf $MAF --recode --out $dmrGWA_base
$plink --file $gwa_base --allow-extra-chr --indep-pairwise 50 5 $pruneLD_dmrGWA --out ${wDir}/tmp2
$plink --file $gwa_base --allow-extra-chr --extract ${wDir}/tmp2.prune.in --recode --out $dmrGWA_base
rm ${wDir}/tmp2.*
$plink --file $dmrGWA_base --allow-extra-chr --recode A --out $dmrGWA_base			# Make raw file for GWAS
awk 'OFS="\t"{if($1~/^Scaffold/){gsub(/Scaffold_/,"",$1);$1=($1+7)};print}' ${dmrGWA_base}.map > ${dmrGWA_base}_chrint.map

### 5) Summary
echo "Report from vcfGWAformat_IBSmx.sh, this script starts from already imputed 10NA 1MAF filtered variants" > $filter_report
echo $(wc -l ${fin_base}.map | cut -d" " -f1) "variants after 0.1NA, biallelic, 1MAF filtering" >> $filter_report
echo $(wc -l ${ibs_base}.map | cut -d" " -f1) "variants after 0.1NA, biallelic, 1MAF filtering and pruning (--indep-pairwise 50 5 0.8). These were used for IBS matrix" >> $filter_report
echo $(wc -l ${gwa_base}.map | cut -d" " -f1) "variants after 0.1NA, biallelic," $MAF"MAF filtering. These variants are ready for GWAS" >> $filter_report
echo $(wc -l ${ind_tests}.prune.in | cut -d" " -f1) "variants after 0.1NA, biallelic," $MAF"MAF filter and pruning (--indep-pairwise 100 5 0.3). This is the number of independent tests in GWAS" >> $filter_report
echo $(wc -l ${dmrGWA_base}.map | cut -d" " -f1) "variants after 0.1NA, biallelic," $MAF"MAF filter and pruning (--indep-pairwise 50 5" $pruneLD_dmrGWA"). These pruned variants are ready for GWAS with DMRs" >> $filter_report

### 6) OPTIONAL: Count homozygous, hetero and... sites (using pruned SNPs nused for the IBS matrix)
$plink --file $ibs_base --allow-extra-chr --recode A --out $ibs_base
awk 'OFS="\t"{if(NR==1){print "sample","homo_ref","het","homo_alt"}
 else{href=0;het=0;halt=0;for(i=7;i<=NF;i++){if($i==0){href++}else if($i==1){het++}else{halt++}};print $1,href,het,halt}}' ${ibs_base}.raw > $zigosity_report
 
### 7) Cleanup
rm ${wDir}/*.log ${wDir}/*.nosex ${wDir}/temp_ready.vcf


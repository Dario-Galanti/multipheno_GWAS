#PBS -l nodes=1:ppn=4 #Nodes and cores
#PBS -l walltime=40:00:00
#PBS -l mem=40Gb 
#PBS -S /bin/bash
#PBS -N GWAS #Job name
#PBS -j oe
#PBS -q short
#PBS -o /beegfs/work/bbmdg01/GWAS/logs/aphid_v3_GWAS.out
#PBS -e /beegfs/work/bbmdg01/GWAS/logs/aphid_v3_GWAS.err

## Author: Dario Galanti (Jan 2021)
## Aim: Run GWAS on multiple phenotypes with "gwas_rr_kinship_multipheno.R" on BinAC cluster (PBS queuing system). Optionally also run enrichment of a-priori candidate variants
## Input 1): multi-phenotypes txt file. Headers are "id" "whatever" "phenotype1" "phenotype2" ... NAs are allowed. NB: "id" names have to be the same as the kinship matrix
## Input 2): output directory name (individual phenotype directories will be created within this directory)
## Input 3) OPTIONAL: kinship: genetic similarity matrix (kinship or IBS). Sample IDs should correspond to phenotypes.txt. If not provided kinship matrix will be calculated
## Input 4) OPTIONAL: flanked_candidates: Bed file containing flanked (we used 20kb) regions of a-priori candidates
## Run: bash sderio_gwas_multipheno_BinAC.sh multiphenotypes.txt outdir
## Run BinAC: qsub -q short -F "pheno_v3/Multi_pheno_GC2020.txt GC2020_multitrait" sderio_gwas_multipheno_BinAC.sh
## Dependencies: conda activate R (Important to use R version downloaded on conda "R")

## SETUP DIR AND DEFINE ALL TOOLS AND SCRIPTS
cd /beegfs/work/bbmdg01/GWAS
Rscript=~/miniconda3/envs/R/bin/Rscript
GWAS_script=gwas_rr_kinship_multipheno.R
# For enrichment
bedtools=~/miniconda3/envs/bedtools/bin/bedtools
flanked_candidates=/beegfs/work/bbmdg01/GWAS/candidates/flanked_Methylation_Ta_genes_final.bed			# Bed file containing flanked (we used 20kb) regions of a-priori candidates
enrich_script=/beegfs/work/bbmdg01/GWAS/candidates/enrichment_plot.R

## DEFINE INPUT FILES
geno_file="input_v3/Ta_v3_vrts_4MAF_imp_GWAspls_withref.raw"		# 5MAF filt after imputation (1422472). More accurate!
snps="input_v3/Ta_v3_vrts_4MAF_imp_GWAspls_withref_chrint.map"
pheno_file=$1
kinship="input_v3/Ta_v3_vrts_1MAF_imp_8prun_GWAspls_withref.kinship"	# OPTIONAL PARAMETER. This K was estimated from 1MAF, 0.8 pruned SNPs (more accurate!!!)
outdir=$2

mkdir -p $outdir

## RUN GWAS
$Rscript --vanilla $GWAS_script genotype_file=$geno_file snp_map=$snps phenotype_file=$pheno_file outdir=$outdir kinship=$kinship		# Kinship


## OPTIONAL 1): ENRICHMENT OF A-PRIORI CANDIDATES
for f in ${outdir}/*/Ta_v3_vrts_*GWAS_rrBLUP.results
do
 enrich_bed=$(dirname $f)/$(basename $f | cut -d"_" -f7-9)_enrichanalysis.bed
 enrichment=$(dirname $f)/$(dirname $f | rev | cut -d/ -f1 | rev)_enrichment_analysis.txt
 tail -n+2 $f | awk -F, 'OFS="\t"{split($1,loc,"_");if($1 ~ /^Scaffold/){print loc[1]"_"loc[2],($3-1),$3,$4}else{print loc[1],($3-1),$3,$4}}' | $bedtools intersect -a stdin -b $flanked_candidates -wa -c > $enrich_bed
 maxi=$(awk 'BEFORE{m=0} {if($4>m){m=$4}} END {print int(m*10)}' $enrich_bed)
 exp_freq=$(awk '{if($5>0){cand++}} END {print cand/NR}' $enrich_bed)
 echo -e -logP"\t"Cand_sig_vrts"\t"Sig_vrts"\t"Observed_freq"\t"Expected_freq"\t"Enrichment > $enrichment
 for ((p=0;p<=maxi;p+=5));
 do
  awk -v logp=$p -v exp_fq=$exp_freq 'BEGIN{logp=(logp/10)}{if($4>=logp){sig++;if($5>0){cand++}}}
   END{OFS="\t";obs_fq=(cand/sig);print logp,cand,sig,obs_fq,exp_fq,(obs_fq/exp_fq)}' $enrich_bed >> $enrichment
 done
 rm $enrich_bed
 $Rscript --vanilla $enrich_script $enrichment
done


## OPTIONAL 2): FORMAT RESULTS TO DOWNLOAD AND VISUALIZE IN IGV (https://software.broadinstitute.org/software/igv/GWAS)
## I remove all SNPs with -log(p)<1 to make the results lighter and easier to download and I format them in the correct way to use them in IGV
for f in ${outdir}/*/Ta_v3_vrts_*GWAS_rrBLUP.results;
do
 newdir=${outdir}/${outdir}_4download/$(echo $f | rev | cut -d"/" -f2 | rev)
 mkdir -p ${newdir}
 IGVfout=${newdir}/$(basename $f | cut -d"_" -f7-9)_IGVresults.gwas
 echo -e "CHR\tBP\tSNP\tP" > $IGVfout
 tail -n+2 $f | awk -F"," 'OFS="\t"{if($4>1){split($1,loc,"_");p=(10^(-$4));if($1 ~ /^Scaffold/){print loc[1]"_"loc[2],loc[3],$1,p}else{print loc[1],loc[2],$1,p}}}' >> $IGVfout
 plot1=$(dirname $f)/FDR_*manhattan_rrBLUP.png
 plot2=$(dirname $f)/FDR_*qqplot_rrBLUP.png
 plot3=$(dirname $f)/*enrich_FDR.pdf
 cp $plot1 ${newdir}
 cp $plot2 ${newdir}
 cp $plot3 ${newdir}
done





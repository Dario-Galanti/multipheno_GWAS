#!/bin/bash
## Author: Dario Galanti (July 2022)
## Aim: Run GWAS with GEMMA with plink input files
## Input: 1) PLINK format input files for GEMMA (.bed, .bim and .fam)
## Input: 2a) template phenotype file in .fam format and 2b) real phenotype file with headers "id" "empty" "pheno1" "pheno2" ...
## Input: 3) kinship matrix
## Documentation: https://www.xzlab.org/software/GEMMAmanual.pdf
## Tutorial: http://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf
## Run: bash gemmaGWA_multipheno.sh phenotype_file.txt output_directory
## Run: bash gemmaGWA_multipheno.sh /beegfs/work/bbmdg01/GWAS/pheno_v5/SupFamTIPs_withUS_1out_noAM.txt SupFamTIPs_1out_noAM

##NB: GEMMA runs with two kinds of input files:
## - PLINK -> Easy to obtain and more symplistic
## - BIMBAM -> It stores dosage information or posterior mean genotypes (Is a value between 0 to 2 that can be interpreted as the
## minor allele dosage: 0 is homozygous for the major allele, 1 is a heterozygote, and 2 is a homozygote for the minor allele.)

## PROCESS AND PRINCIPLES:
## GEMMA cannot be run with multiple threads, so even though they have a model that runs multiple phenotypes, it is very slow
## So we run the single phenotype lmm and we parallelise submitting each phenotype as a different job to the queue
## This requires having multiple copies of the input files, because a single file basename has to be provided for plink.
## So we create multiple copies of the plink input files where only the phenotype .fam file differs 


## DEFINE TOOLS
Rscript=~/miniconda3/envs/R4/bin/Rscript
gemma=~/miniconda3/envs/GEMMA/bin/gemma
bedtools=~/miniconda3/envs/bedtools/bin/bedtools
work=/beegfs/work/bbmdg01/GEMMA
#qqman=${work}/qqman_GEMMA.R
qqman=${work}/qqman_GEMMA_grey.R
enrich_script=/beegfs/work/bbmdg01/GEMMA/enrichment_plot.R
cd /beegfs/work/bbmdg01/GEMMA		# This allows to use relative paths

#flanked_candidates=/beegfs/work/bbmdg01/GWAS/candidates/flanked_Methylation_TaV5_genes_final.bed		# Methylation
flanked_candidates=/beegfs/work/bbmdg01/GWAS/candidates/DefResp_AtOrtho_genes_BSMT1_20kbFlk.bed	# Defense response (from At Ortho with GO:0006952) + BSMT1 genes

## DEFINE INPUT
## NB: Using plink format has the stupid problem that the phenotype file cannot be defined alone but needs to have the same file prefix of the SNP file
## So we use a workaround: We define the input file and copy it to same dir and filename prefix of the SNP file
pheno_file=$1																		# /beegfs/work/bbmdg01/GWAS/pheno_v5/SupFamTIPs_withUS_GWA_log_1out.txt
pheno_fam=${work}/plink_pheno/$(basename $1 .txt).fam								# Will be created and stored in the "plink_pheno" dir
plink_base=${work}/plink_SNPs/Ta_v5_vrts_10NA_biall_4MAF_imp_withUSspls_withref		# 4MAF filt after imputation (1856975). More accurate!
#plink_base=${work}/plink_SNPs/Ta_v3v5_vrts_10NA_4MAF_bial_imp_GWAspls_withref		# 4MAF filt after imputation. More accurate!
plink_bed=${plink_base}.bed
plink_bim=${plink_base}.bim
plink_fam=${plink_base}.fam
## DEFINE K
kinship=${work}/IBS/Ta_v5_vrts_10NA_biall_1MAF_imp_8prun_withUSspls_withref.kinship
k_nohead=${work}/IBS/$(basename $kinship .kinship)_nohead.kinship
## DEFINE OUTPUT
outdir=${work}/$2
download_dir=${outdir}/${2}_4download
mkdir -p $download_dir
report=${outdir}/gemmaGWA_report.txt

## 1) CHECK AND FORMAT INPUT
## We use the first 5 fields of the original .fam file and paste all the phenotypes to it, so that gemma finds the correct plink input
pheno_arr=($(head -1 $pheno_file | cut -f3- | tr -d '\r'))
num_pheno=${#pheno_arr[@]}
rm $pheno_fam							# Clean start (FIX SO IT DOESN'T THROW AN ERROR!!!!!!!!!!!!!!!!!!!)
while read l; do
 spl=$(echo $l | cut -d" " -f1)
 # grep values from phenotype file and replace empty fields with NAs (unix newline is \n, while DOS/windows is \r\n)
 val=$(grep -P $spl"\t" $pheno_file | cut -f3- | tr -d '\r' | awk -F"\t" '{OFS=" ";for(i=1;i<=NF;i++){if($i==""){$i="NA"}};print}')
 if [[ ${#val} == 0 ]];then val=$(seq $num_pheno | sed "c NA" | tr "\n" " ");fi
 newline=$(echo $l | cut -d" " -f1-5)" "$val
 echo $newline >> $pheno_fam
done < $plink_fam
tail -n+2 $kinship > $k_nohead		#NB: I should implement a check for sample names: they should be the same as in the plink files!!!!!!!!!

## PRINT SOME INFO
echo plink input: $plink_base > $report
echo phenotype file: $pheno_file >> $report
echo phenotypes: ${pheno_arr[@]} >> $report
echo kinship: $kinship >> $report

## LOOP THROUGH PHENOTYPES AND SUBMIT ONE JOB FOR EACH
mkdir -p ${work}/work/${2}
mkdir -p ${work}/logs/${2}
n=0
for phe in ${pheno_arr[@]};
do
	n=$((n+1))
	c=$((n+5))
	jobName=${work}/work/${2}/gemma.${phe}.sh
	(
	echo "#PBS -l nodes=1:ppn=2 #Nodes and cores"
	echo "#PBS -l walltime=10:00:00"
	echo "#PBS -l mem=8Gb"
	echo "#PBS -S /bin/bash"
	echo "#PBS -N gemma.${phe}"
	echo "#PBS -j oe"
	echo "#PBS -q short"
	echo "#PBS -o ${work}/logs/${2}/gemma.${phe}.out"
	echo "#PBS -e ${work}/logs/${2}/gemma.${phe}.err"
	echo ""
	
	## Define and prepare input and output
	echo "ph_dir=${outdir}/${phe}"
	echo "indir=\${ph_dir}/input"
	echo "plink_phe_base=\${ph_dir}/input/$(basename $plink_fam .fam)"
	echo "results=\${ph_dir}/${phe}.assoc.txt"
	echo "mkdir -p \$indir"
	echo "cd \$ph_dir"
	echo "cp $plink_bed \${indir}/"
	echo "cp $plink_bim \${indir}/"
	echo "cut -d ' ' -f1-5,${c} $pheno_fam > \${indir}/$(basename $plink_fam)"
	echo ""

	## 2) RUN GWAS
	echo "$gemma -bfile \$plink_phe_base -k $k_nohead -lmm 1 -o $phe"
	echo "mv output/* ."													# Move GEMMA output to $ph_dir
	echo "rm -r input output"
	
	## 3) VISUALIZATION (manhattan and qqplot)
	echo "$Rscript --vanilla $qqman \$results"
	echo ""
	
	## 4) ENRICHMENT OF A-PRIORI CANDIDATES
	echo "enrich_bed=Enrichment_${phe}.bed"
	echo "enrichment=Enrichment_${phe}.txt"
	echo "awk 'NR>1{OFS=\"\t\";logp=-log(\$12)/log(10);print \$1,(\$3-1),\$3,logp}' \$results | $bedtools intersect -a stdin -b $flanked_candidates -wa -c > \$enrich_bed"
	echo "maxi=\$(awk 'BEFORE{m=0} {if(\$4>m){m=\$4}} END {print int(m*10)}' \$enrich_bed)"
	echo "exp_freq=\$(awk '{if(\$5>0){cand++}} END {print cand/NR}' \$enrich_bed)"
	echo "echo -e -logP'\t'Cand_sig_vrts'\t'Sig_vrts'\t'Observed_freq'\t'Expected_freq'\t'Enrichment > \$enrichment"
	echo "for ((p=0;p<=maxi;p+=5));"										# Iterate through -log(p) by 0.5 increasing steps and calculate enrichment
	echo "do"
	echo " awk -v logp=\$p -v exp_fq=\$exp_freq 'BEGIN{logp=(logp/10)}{if(\$4>=logp){sig++;if(\$5>0){cand++}}} END{OFS=\"\t\";obs_fq=(cand/sig);print logp,cand,sig,obs_fq,exp_fq,(obs_fq/exp_fq)}' \$enrich_bed >> \$enrichment"
	echo "done"
	echo "rm \$enrich_bed"
	echo "$Rscript --vanilla $enrich_script \$enrichment"
	echo ""
	
	## 5): FORMAT RESULTS FOR IGV (https://software.broadinstitute.org/software/igv/GWAS) AND PREPARE FOLDER FOR DOWNLOAD
	## I remove all SNPs with -log(p)<1 to make the results lighter and easier to download and I format them in the correct way to use them in IGV
	echo "newdir=${download_dir}/${phe}"
	echo "mkdir \$newdir"
	echo "IGVfout=\${newdir}/${phe}_IGVresults.gwas"
	echo "echo -e CHR'\t'BP'\t'SNP'\t'P > \$IGVfout"
	echo "tail -n+2 \$results | awk 'OFS=\"\t\"{if(\$12<0.1){print \$1,\$3,\$2,\$12}}' >> \$IGVfout"
	echo "cp ./*.png \${newdir}"
	echo "cp ./*.pdf \${newdir}"
	
	) > $jobName
		chmod +x $jobName
		echo "Bash file $jobName created"
		qsub -q short ${jobName}
done
exit


## Author: Dario Galanti 
## Aim: R script to carry out multi-phenotype GWAS with the package rrBLUP, using mixed models
## kinship matrix is used to account for population structure. It can be provided as "kinship=K_matrix.kinship" or it will be calculated from the genotype file.
## Input geno: Plink .raw and .map files NB: Map file should have Chrom/Scaff as integers!!
## Input pheno: phenotype txt file with headers: id <TAB> pop <TAB> phenotype1 <TAB> phenotype2 ... NAs are allowed
## Run: Rscript --vanilla gwas_rr_kinship_multipheno.R genotype_file=path_to_genotypes snp_map=path_to_map phenotype_file=path_to_phenotypes 
## Run: Rscript --vanilla gwas_rr_kinship_multipheno.R genotype_file=Ta_v2_gwas_ready.raw snp_map=Ta_v2_gwas_ready.map phenotype_file=pheno/grpix_samples_newnames.txt

## Install dependencies:
## NB: If running without a kinship matrix, the script needs the "gMatrix" package, which has to be installed from source file (https://cran.r-project.org/src/contrib/Archive/gmatrix/)

## IMPORTANT: PHENOTYPE TRANSFORMATION !!!!!
## If phenotype is deviating strongly from normality (shapiro.test pvalue < 0.0001), a transformation should be applied.
## BY DEFAULT NO TRANSFORMATION IS APPLIED, BUT TWO OPTIONS CAN BE ACTIVATED:
## OPTION 1: Inverse Normal Transformation (Brute force transformation, the phenotype will become perfectly normal)
## OPTION 2: "Gaussianise" function can be applied instead. Better option even though the resulting phenotype might still not be completely normal
## To activate the options, see the "TRANSFORMATION" part of the script, around line 130

## OPTIONAL: You can define number of independent tests. If not define will be calculated = n° variants tested
## This can be done by SNP pruning with Plink according to the study "Addressing population-specific multiple testing burdens in genetic association studies". 
## Independent tests = pruning (--indep-pairwise 100 5 0.3)
#tests <- 143346 # 4MAF. Pruning (--indep-pairwise 100 5 0.3)

library("qqman")
library("dplyr")
library("rrBLUP")
library("gMatrix")     # Only necessary if no genetic distance matrix is provided
library("data.table")
library("moments")     # For investigating skewness
library("LambertW")    # OPTIONAL: For transforming the phenotypes with Gaussianize function

print("GWAS using the rrBLUP package")

###################################
## read arguments from command line
###################################
allowed_parameters = c(
  'genotype_file',
  'snp_map',
  'phenotype_file',
  'kinship',
  'outdir'
)

args <- commandArgs(trailingOnly = TRUE)

### 1) FORMATTING INPUT PARAMETERS AND ASSIGNING THEM TO VARIABLES
for (p in args){
  pieces = strsplit(p, '=')[[1]]
  #sanity check for something=somethingElse
  if (length(pieces) != 2){stop(paste('badly formatted parameter:', p))}
  if (pieces[1] %in% allowed_parameters)  {
    assign(pieces[1], pieces[2])
    next
  }
  stop(paste('bad parameter:', pieces[1]))
}

print(paste("genotype file name:",genotype_file))
print(paste("SNP map:",snp_map))
print(paste("phenotype file name:",phenotype_file))
if (exists("kinship")) {print(paste("kinship:",kinship))}
print(paste("outdir:",outdir))

dataset <- basename(genotype_file)
dataset <- substr(dataset, 1, nchar(dataset)-4)

report <- paste(outdir,"/GWASreport.txt",sep="") #Only for tracking transformations, for everything else check the logs

### 2) READING GENOTYPE DATA
print("now reading in the data ...")
## genotypes
X <- fread(genotype_file, header = TRUE)
print(paste(nrow(X),"records read from the genotype file",sep=" "))
SNP_INFO <- fread(snp_map)
names(SNP_INFO) <- c("Chr","SNP","cM","Pos")
SNP_INFO$cM <- NULL

X_IDs <- (X[,c(1:6)])
X <- as.matrix(X[,-c(1:6)])
colnames(X) <- gsub("\\_[A-Z]{1,}$","",colnames(X)) #delete SNP/INDEL extention from variant name (probably useless)
rownames(X) <- X_IDs$IID

SNP_INFO <- bind_cols(SNP_INFO,as.data.frame(t(X)))
SNP_INFO <- as.data.frame(SNP_INFO)

print(paste(nrow(SNP_INFO),"SNPs read from the map file",sep=" "))

## Check whether .raw and .map files have the same number of markers
if ((ncol(X)) != nrow(SNP_INFO)) {
  stop("!! N. of SNPs in the map file not equal to the number of genotyped SNPs in the genotype file")
} else print("N. of SNPs in the map and genotype files is the same: this is correct!!")

### 3) CALCULATING KINSHIP MATRIX (if not already provided!!)
if (exists("kinship")) {
  print(paste("Using provided kinship matrix:", kinship))
  K <- as.matrix(fread(kinship))
  rownames(K) <- colnames(K) #Because the kinshim matrix is printed out without rownames
} else {
  print("No kinship matrix was provided, so we calculate it!!")
  K <- gVanRaden.2(X) #Requires to install introduction_to_gwas/software/gMatrix_0.2.tar.gz
  print("writing out the kinship matrix ...")
  fname = paste(outdir,"/",dataset,".kinship",sep="")
  write.table(K, file=fname, quote = FALSE, row.names = FALSE)
  print("producing the heatmap kinship matrix ...")
  pdf(paste(outdir,"/",dataset,"_kinship_heatmap",".pdf",sep=""), width=30, height=30) # These dims fit ~200 sample names
  heatmap(K,col=rev(heat.colors(75)))
  dev.off()
}

### 4) READING PHENOTYPES
phenotypes <- fread(phenotype_file, data.table=FALSE)
if(names(phenotypes)[1] != "id"){stop("id column header has to be \"id\" and it is not!!!")}

### 5) ITERATE THROUGH PHENOTYPES AND RUN GWAS FOR EACH
for (p in 3:ncol(phenotypes)){
  # Read and format phenotype
  trait <- colnames(phenotypes)[p]
  pheno <- phenotypes[!is.na(as.numeric(phenotypes[,p])),c(1,p)]
  # Intersect phenotype individuals with genotype individuals
  pheno <- pheno[pheno$id %in% X_IDs$IID,]
  print(paste("Running", trait, "with", nrow(pheno), "individuals", sep=" "))
  # Intersect genotype individuals and K matrix individuals with phenotype individuals
  vec <- colnames(K) %in% pheno$id
  Kinship <- K[vec,vec]
  SNP_info <- SNP_INFO[,c(TRUE,TRUE,TRUE,vec)]
  
## TRANSFORMATION!!! (IMPORTANT STEP!!!!)
## We use Gaussianize function from LambertW to transform the phenotype if deviating from normality
## We hope this works better than INT, as the latter reduces power in my dataset!!!
  norm <- shapiro.test(as.numeric(pheno[,2]))
  if(norm$p.value < 0.0001){
    write(paste(names(pheno)[2],"deviates heavily from normality: shapiro.p=",norm$p.value, "Consider applying transformation"), file=report, append=TRUE)
    
## OPTION 1: Inverse Normal Transformation (Remove #s to apply)
#    write(paste(trait,"deviates from normality: shapiro.p=",norm$p.value, "INT was applied"),file=report, append=TRUE)
#    pheno$INTpheno <- qnorm((rank(pheno[,2],na.last="keep")-0.5)/sum(!is.na(pheno[,2])))
#    names(pheno)[3] <- names(pheno)[2]
#    pheno[,2] <- NULL
    
## OPTION 2: Gaussianize (Remove #s to apply)
#    sk <- skewness(as.numeric(pheno[,2]))  # Inspect skwness to decide on the Gaussianize "type" parameter
#    if( sk > -1 && sk < 1){
#      if( sk > -0.5 && sk < 0.5){
#        write(paste(names(pheno)[2],"not skewed"),file=report, append=TRUE);type="h"} else {
#          write(paste(names(pheno)[2],"mildly skewed"),file=report, append=TRUE);type="s"}} else {
#            write(paste(names(pheno)[2],"heavily skewed"),file=report, append=TRUE);type="hh"}
#    write(paste(trait,"deviates heavily from normality: shapiro.p=",norm$p.value, "Gaussianize with type=",type,"was applied"),file=report, append=TRUE)
#    pheno$GAUSSpheno <- Gaussianize(data=as.numeric(pheno[,2]), type=type)
#    names(pheno)[3] <- paste(names(pheno)[2],"_gaus",sep="")
#    pheno[,2] <- NULL
    
  } else {write(paste(names(pheno)[2],"is close to normal distribution, no transformation necessary"), file=report, append=TRUE)}

  ## RUN GWAS FOR EACH PHENOTYPE (https://cran.r-project.org/web/packages/rrBLUP/rrBLUP.pdf)
  model1 <- GWAS(
    pheno = pheno,
    geno = SNP_info,
    # n.PC = 4, #Choose between this, nothing or K. PCs are used as fixed factors
    # fixed = #array of phenotype_file headers containing fixed factors to correct the phenotype
    K = Kinship, #Control for pop str with K
    #n.core = 4, #Somehow this makes it slower instead of faster
    plot = FALSE
  )
  names(model1)[length(model1)] <- trait
  gwasResults <- model1[,c("SNP","Chr","Pos",trait)]
  names(gwasResults) <- c("SNP","CHR","BP","P")
  gwasResults <- gwasResults[order( gwasResults$CHR, gwasResults$BP),]
  ## SAVE RESULTS
  dir.create(paste(outdir,"/",trait, sep=""))
  fout <- paste(outdir,"/",trait,"/",dataset,"_",trait,"_GWAS_rrBLUP.results", sep="")
  fwrite(x = gwasResults, file = fout)
  
  ## MANHATTAN PLOT
  topval <- max(gwasResults$P)
  gwasResults$P = 10^((-gwasResults$P)) #rrblup results already contain -log(p)
  SNPs <- nrow(gwasResults)
  if (exists("tests")==FALSE) {tests <- SNPs}
  Bonferroni <- -log10(0.05/tests)
  pvec <- sort(gwasResults$P)
  FDRtr <- 6.0
  for (i in 1:length(pvec)) {
    P <- pvec[i]
    cv <- (i/SNPs)*0.20
    if(P < cv) { FDRtr <- -log10(P)
    }
  }
  print(paste("Bonferroni treshold for",tests,"indipendent tests: ", Bonferroni))
  print(paste("FDR 0.20 treshold for",SNPs,"SNPs:",FDRtr))
  
  ## MAKING MANHATTAN
  ytop <- max(Bonferroni, topval)
  ## "cex" regulates dot size; "cex.axes" axes numbers; "cex.lab" axes lables
  #man.pal <- c("blue4","orange3") # OLD
  man.pal <- c("darkblue","darkcyan") # NEW
  png(paste(outdir,"/",trait,"/FDR_0.20_",trait,"_manhattan_rrBLUP.png",sep=""), width=1400, height=600)
  par(mar = c(5, 5.2, 4, 2) + 0.1) # Add space for the y axes lable
  manhattan(gwasResults, cex=1.6, cex.axis=1.7, cex.lab=1.9, col=man.pal, main=trait,
            suggestiveline = FDRtr, genomewideline = Bonferroni, 
            ylim = c(0, ceiling(ytop)))
  dev.off()
  
  ## MAKING QQ-PLOT
  png(paste(outdir,"/",trait,"/FDR_0.20_",trait,"_qqplot_rrBLUP.png",sep=""), width = 500, height = 500)
  par(mar = c(5, 5.2, 4, 2) + 0.1) # Add space for the y axes lable
  qq(gwasResults$P, cex = 1.2, cex.axis = 1.6, cex.lab=1.7, col="darkblue")
  dev.off()
}


print("#########")
print("## END GWAS ##")
print("#########")


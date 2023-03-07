## Author: Dario Galanti
## Aim: Make manhattan and qqplot from GEMMA output file
## INPUT: gwas output from GEMMA
## RUN: Rscript --vanilla qqman_GEMMA.R gemma_output/phenotype.assoc.txt
## RUN: Rscript --vanilla qqman_GEMMA.R SupFamTIPs_log_1out_RIR/SupFamTIPs_log_1out_RIR.assoc.txt

library(qqman)
library(fastman)  # For plot with larger lables
library(data.table)

## IMPORTANT: Define number of independent tests. This is done by SNP pruning with Plink according to the study
## "Addressing population-specific multiple testing burdens in genetic association studies"
## Independent tests = Variants passing 0.1NAs, 0.05MAF, biallelic, pruning (--indep-pairwise 100 5 0.4)
#tests <- 148404 # 4MAF. Pruning (--indep-pairwise 100 5 0.3) v4v5 lifted withUSspls (withref)

args <- commandArgs(trailingOnly = TRUE)

# Results <- fread("SupFamTIPs_log_1out_RIR.assoc.txt",header=T, data.table = F, sep="\t")
Results <- fread(args[1], header=T, data.table = F, sep="\t")
file <- basename(args[1])
file <- substr(file, 1, nchar(file)-10)
outdir <- dirname(args[1])

## FORMAT FOR QQMAN
Results <- Results[c(2,1,3,12)]
colnames(Results) <- c("SNP","CHR","BP","P")
Results$CHR <- as.integer(gsub("Scaffold_","",Results$CHR))

# MANHATTAN PLOT
topval <- -log10(min(Results$P)) #save top value for manhattan ylim scaling
#Results$P=10^((-Results$P)) #rrblup results already contain -log(p), convert to p
"##CALCULATE BONFERRONI THRESHOLD"
SNPs <- nrow(Results)
if (exists("tests")==FALSE) {tests <- SNPs}
Bonferroni <- -log10(0.05/tests)

"##CALCULATE FDR TRESHOLD"
pvec <- sort(Results$P)
FDRtr <- 6.0
for (i in 1:length(pvec)) {
  P <- pvec[i]
  cv <- (i/SNPs)*0.2
  if(P < cv) { FDRtr <- -log10(P)
  }
}
print(paste("Bonferroni treshold for",tests,"indipendent tests: ", Bonferroni))
print(paste("FDR 0.20 treshold for",SNPs,"SNPs:",FDRtr))

##MAKING MANHATTAN
ytop <- max(Bonferroni, topval)
man.pal <- c("grey50","grey20") # grey
top.pal <- c("darkorange2", "orangered2")
#top.pal <- c("royalblue2", "blue")   # blues


# Remove very high pvalues for nicer plotting
Results <- Results[Results$P<0.85,]

## FASTMAN WITH LARGER LABLES
png(paste(outdir,"/FDR_0.2_",file,"_manhattan_GEMMA.png",sep=""), width = 1400, height = 600)
par(mar = c(5, 5.4, 4, 2) + 0.1) # Add space for the y axes label. c(bottom, left, top, right) margins with default c(5,4,4,2)+0.1
par(mgp=c(3.1,1.1,0)) #Add space between axes lables and plot. c(axis.title, axis.label, axis.line) with default c(3,1,0)
fastman(Results, cex=1.6, cex.axis=2.1, cex.lab=2.1, main=file, maxP=NULL,
          genomewideline = Bonferroni,
          #suggestiveline = FDRtr,
          #col = man.pal,
          colAbovePval=T, annotatePval=Bonferroni, col2=man.pal, col=top.pal, annotationCol="white",
          ylim = c(0, ceiling(ytop)))
dev.off()

## QQPLOT
png(paste(outdir,"/",file,"_qqplot_GEMMA.png",sep=""), width = 500, height = 500)
par(mar = c(5, 5.2, 4, 2) + 0.1) # Add space for the y axes lable
qq(Results$P, cex = 1.4, cex.axis = 1.7, cex.lab=1.8, col="darkblue")
dev.off()

print("#########")
print("## END QQMAN ##")
print

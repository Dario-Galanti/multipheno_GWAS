## Author: Dario Galanti April 2021
## Aim: Plot a-priori candidate SNP enrichment
## Input: txt file with enrichment analysis for increasing -log(p) values
## Run: Rscript --vanilla enrichment_plot.R GBM_avg_more5meth/CpG_avg_GBM_more5avmeth/CpG_avg_GBM_more5avmeth_enrichment_analysis.txt

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

file <- basename(args[1])
file <- substr(file, 1, nchar(file)-17) # Remove extra stuff from file name
outdir <- dirname(args[1])

## 1) READ DATA AND TIDY
enrichment <- read.table(args[1], header=TRUE, sep="\t")
#enrichment <- read.table("CpG_avg_GBM_more5avmeth_enrichment_analysis.txt", header=TRUE, sep="\t")

## Remove lines with observed frequency = 0
enrichment <- enrichment[!(enrichment$Observed_freq == 0),]

## 2) CALCULATE FDR as in "Atwell et al. 2010" (https://static-content.springer.com/esm/art%3A10.1038%2Fnature08800/MediaObjects/41586_2010_BFnature08800_MOESM250_ESM.pdf)
Cand_vrts <- enrichment[1,2]
Tot_vrts <- enrichment[1,3]
Non.cand <- (Tot_vrts - Cand_vrts)
enrichment$Non.cand_sig_vrts <- (enrichment$Sig_vrts - enrichment$Cand_sig_vrts)
enrichment$x <- enrichment$Cand_sig_vrts/Cand_vrts
enrichment$y <- enrichment$Non.cand_sig_vrts/Non.cand
enrichment$FDR <- (1 -(enrichment$x - enrichment$y)/enrichment$x)
enrichment$FDR <- ifelse(enrichment$FDR > 1, 1, enrichment$FDR) # Limit FDR to max 1. Cannot be higher by definition

## 3a)  PRINT ENRICHMENT AND FDR PLOT
max_enrich <- ceiling(max(enrichment$Enrichment))

#png(paste(outdir,"/",file,"_FDR.png",sep=""), width = 550, height = 450) # In png
pdf(paste(outdir,"/",file,"_FDR.pdf",sep=""), width = 7, height = 6)   # In pdf
ggplot(enrichment, aes(x=X.logP)) +
  geom_line(aes(y=Enrichment), size=1.2, color="dodgerblue4") +
  geom_line(aes(y=FDR*max_enrich), size=1.2, color="firebrick") +
  scale_x_continuous(name = "-log(p)") +
  scale_y_continuous(
    name = "Enrichment",
    sec.axis = sec_axis( trans=~./max_enrich, name="FDR")
  ) +
  #theme_minimal() +
  theme_light() +
  theme(
    axis.title.x = element_text(hjust = 0.5, size=21),
    axis.text.x = element_text(size=18),
    axis.title.y = element_text(color = "dodgerblue4", hjust = 0.5, size=21),
    axis.text.y = element_text(color="dodgerblue4", size=18),
    axis.title.y.right = element_text(color = "firebrick", hjust = 0.5, vjust = 1.5, size=21),
    axis.text.y.right = element_text(color="firebrick", size=18)
  )
dev.off()








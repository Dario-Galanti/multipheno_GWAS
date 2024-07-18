## Author: Dario Galanti July 2024
## Aim: Perform genome rotation scheme to calculate significance of a given enrichment of a priori candidate genes from GWAS results.
## Process: Starting from a file with SNPs, pvalues and a-priori candidate status as (0/1), randomly rotate pvalues and candidate status within each chromosome and calculate enrichment for 1M permutations. This gives the null distribution of enrichments. Test if the observed enrichment is significantly different from the null distribution.
## Input: bed file with: Chr Start End pvalue Candidate_status(0/1)
## Run: Rscript --vanilla enrich_sig_genome_rotation.R ErysiRes_enrichanalysis.bed enrich_sig.txt $iterations $logs_file

## NB: I am randomly rotating pvalues and candidate status within each chromosome
## but for now I am only rotating and not randomly inverting. That can be easily implemented.

library(data.table)
library(dplyr) #for bind_rows()
library(foreach)   # for parallelization
library(doParallel)  # for parallelization
library(parallelly)   # for availableCores() - to detect available cores on HPC

### 0) Define input and output
args <- commandArgs(trailingOnly = TRUE)
fin <- args[1]  # fin <- "ErysiRes_enrichanalysis.bed"
df <- fread(fin, header = F, data.table = F)
fout <- args[2]  # fout <- "enrich_sig_paral_R.txt"
logP <- as.numeric(args[3])  # Sig. threshold al -log(p) | logP <- 6.424212
iterations <- as.integer(args[4])  # iterations <- 100000
logs <- args[5]


## 1) Split data frame by scaffold
X <- split(df, df$V1)

## OPTIONAL: REMOVE SCAFFOLDS WITH <3 VARIANTS
scaff_keep <- lapply(X, function(df) nrow(df) >= 3) # Logical vector of scaffolds to keep
X <- X[unlist(scaff_keep)]    # Subset the list
df <- dplyr::bind_rows(X)

### CALCULATE EXPECTED FREQ AND ENRICHMENT
exp_freq <- length(which(df$V5>0))/nrow(df)
sig_df <- df[df$V4>logP,]
observed_freq <- length(which(sig_df$V5>0))/nrow(sig_df)     # Identical
#observed_freq <- sum((sig_df$V5>0),na.rm=TRUE)/nrow(sig_df) # Identical
enrichment <- observed_freq/exp_freq


## FUNCTION TO ROTATE THE PVALUES AND A-PRIORI CANDIDATE STATUS OF EACH CHROMOSOME
# Function to rotate the vector by a random number
rotate_vector <- function(vec) {
  n <- length(vec)   # Get the length of the vector
  rand <- sample(0:(n-1), 1)
  if (rand == 0) return(vec)
  c(vec[(n - rand + 1):n], vec[1:(n - rand)])
}
# Test the rotate_vector function
#vec <- 1:10
#rotated_vec <- rotate_vector(vec)


### SETUP FOR PARALLELIZATION
#cores <- availableCores() -1  # For laptop (1 cpu free not to crush laptop)
cores <- availableCores()   # For server
print(paste("detecting",cores,"cores", sep=" "))
cl <- makeCluster(cores)
registerDoParallel(cl) #register the parallel backend with the foreach package

## PRINT INFO TO OUTPUT FILE
## Checking if things are correct before the script is finished
cat(paste("using", cores, "cpus", sep=" "), sep="\n", file=fout)
cat(paste("running",iterations,"iterations", sep=" "), sep="\n", file=fout, append=T)
cat(paste("input file:", fin, sep=" "), sep="\n", file=fout, append=T)
cat(paste("Sig. threshold:", logP, sep=" "), sep="\n", file=fout, append=T)
cat(paste("Exp. freq:", exp_freq, sep=" "), sep="\n", file=fout, append=T)
cat(paste("Obs. freq:", observed_freq, sep=" "), sep="\n", file=fout, append=T)
cat(paste("Enrichment:", enrichment, sep=" "), sep="\n", file=fout, append=T)

### PARALLEL PART: Iterate through dataframes and rotate pvals and states
null_distr <- foreach(i = 1:iterations, .combine = c) %dopar% {
  iter_X <- X
  if (i%%1000 == 0){cat(paste(i,"iterations done",sep=" "),sep="\n",file=logs,append=T)} #check progress
  ### Iterate through scaffolds df and rotate pvals and states
  for (scaf in 1:length(X)){
    iter_X[[scaf]]$V4 <- rotate_vector(X[[scaf]]$V4)
    iter_X[[scaf]]$V5 <- rotate_vector(X[[scaf]]$V5)
  }
  ### Re-join scaffold dataframes
  rotated_df <- dplyr::bind_rows(iter_X)
  
  ### CALCULATE ROTATION ENRICHMENT
  sig_df <- rotated_df[rotated_df$V4>logP,]
  obs_freq <- length(which(sig_df$V5>0))/nrow(sig_df)
  enrich <- obs_freq/exp_freq
  
  # Return the enrichment
  enrich
}

### stop cluster
stopCluster(cl)

### CALCULATE P VALUE OF THE ENRICHMENT
pvalue <- mean(null_distr >= enrichment)

### SAVE P VALUE ON OUTPUT FILE
cat(paste("Enrichment pvalue:", pvalue, sep=" "), sep = "\n", file=fout, append=T)

### IF WILLING TO PRINT NULL DISTRIBUTION
#fwrite(as.data.frame(null_distr), "enrich_null_dist.txt", sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE)




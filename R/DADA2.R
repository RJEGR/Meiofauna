# Running data based on 388 libarries
# Ricardo Gomez-Reyes
# https://github.com/RJEGR/Small-RNASeq-data-analysis/blob/master/DOWNSTREAM_BKP/DIFFEXP/DADA2_iso_miRNA_LAB.R

rm(list=ls())

# source('https://raw.githubusercontent.com/RJEGR/BM/7b306de3fedc539640e475dd157a67db9fdb8a47/filteringAndTriming_learnError_test.R')

args = commandArgs(trailingOnly=TRUE)


# ==============
## Checking and Load packages ----
# ==============
.cran_packages <- c("ggplot2", "GGally", "reshape2", "dplyr","ggExtra")
.bioc_packages <- c("dada2")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# +++++++++++++++

# ================================
## DADA2 Parameters in the config:
# ================================
 
date <- format(Sys.time(), "%Y%m%d")

out_prefix <- "multirun" # Output prefix (for plots and files)
run <- paste0("multirun_", date,"_18S")

# ================
# Outputs in `pwd`:
# ================

path <- "~/MEIOFAUNA/raw-seqs/"

setwd(path)

out_path <- file.path(path, run)

# create it manually if widows

# system(command = paste0("mkdir -p ", out_path), intern = F)


# outdir

system(paste0('mkdir -p ', path, '/outputs'))

outpath <- paste0(path, "/outputs")

# Input dataset

str(fastq.gz <- sort(list.files(path, pattern=".gz$")))

str(fqFs <- sort(fastq.gz[grep('R1', fastq.gz)]))
str(fqRs <- sort(fastq.gz[grep('R2', fastq.gz)]))

# sampleo <- c(1, 85, 100)
# 
# fqFs <- fqFs[sampleo]
# fqRs <- fqRs[sampleo]


# source("https://raw.githubusercontent.com/RJEGR/metagenomics/master/plotQualityProfile.R")

# plotQP(c(fqFs[1], fqRs[1]))

if(length(fqFs) != length(fqRs)) stop("Forward and reverse files do not match.")

# File parsing

system(paste0('mkdir -p ', path, 'filtered'))

filtpathF <- file.path(path, "filtered")
filtpathR <- file.path(path, "filtered")

filt <- sort(list.files(path, pattern="gz"))

filtFs <- sort(filt[grep('R1', filt)])
filtRs <- sort(filt[grep('R2', filt)])

# filtFs <- filtFs[sampleo]
# filtRs <- filtRs[sampleo]

str(snF <- sapply(strsplit(filtFs, "[_]"), `[`, 1))
str(snR <- sapply(strsplit(filtRs, "[_]"), `[`, 1))

snF[!snF %in% snR]

if(!identical(snF, snR)) stop("Forward and reverse files do not match.")

names(filtFs) <- snF
names(filtRs) <- snR


packageVersion("dada2")

# Filtrar de manera homogenea todas las bibliotecas como estandar 

truncLen <- c(150,140)

filterAndTrim <- filterAndTrim(fwd=file.path(path, fqFs), filt=file.path(filtpathF, fqFs),
                               rev=file.path(path, fqRs), filt.rev=file.path(filtpathR, fqRs),
                               truncLen = truncLen, trimRigh = c(100,100),
                               maxEE=2, truncQ=10, maxN=0, rm.phix=TRUE,
                               compress=TRUE, verbose=TRUE, multithread = F)

filtFs <- c(file.path(filtpathF, fqFs)[1], file.path(filtpathR, fqRs)[1])

# plotQP(filtFs)

pct_trim <- 1 - filterAndTrim[,2]/filterAndTrim[,1]

pct_trim

# Learn Error ----
# MAX_CONSIST <- 40
# 
# threads <- NULL
# 
# errF <- learn_Err(filtFs)
# 
# plotErrors(errF, nominalQ=TRUE)

write.table(filterAndTrim, 
            file = paste0(outpath, "/", out_prefix, "_filterAndTrim.tsv"), 
            sep="\t", row.names = F, col.names = T)

# Infer Sequence Variants

sample.groups <- unique(sapply(strsplit(fqFs, "[-]"), `[`, 1))
mergers <- vector("list", length(sample.groups))
names(mergers) <- sample.groups


# Learn error and , derep and run dada

# Learn error and , derep and run dada

for(sam in sample.groups) {
  run <- sapply(strsplit(fqFs, "[-]"), `[`, 1)
  which_samples <- which(run %in% sam)
  
  cat("Processing Run:", sam, "\n")
  
  Fs <- file.path(filtpathF, fqFs)[which_samples]
  Rs <- file.path(filtpathR, fqRs)[which_samples]
  
  if(length(Fs) != length(Rs)) stop("Forward and reverse files do not match.")
  
  samples <- sapply(strsplit(basename(Fs), "[_]"), `[`, 1)
  cat("Processing Files:",  samples, "\n")
  
  cat("learn Errors in Run: ", sam, "\n")
  
  errF <- learnErrors(Fs, nbases=1e8, multithread = threads)
  # Learn reverse error rates
  errR <- learnErrors(Rs, nbases=1e8, multithread = threads)
  
  cat("Dereplicate files: ", samples, "\n")
  
  derepF <- derepFastq(Fs)
  ddF <- dada(derepF, err=errF, multithread = threads)
  derepR <- derepFastq(Rs)
  
  cat("Infering amplicon sequence variance using the error rate learned from run:", sam, "\n")
  
  ddR <- dada(derepR, err=errR, multithread = threads)
  
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  
  propagateCol <- names(ddR)
  
  # nrow(mergers <- mergePairs(dadaFs, derepFs, 
  #                            minOverlap = 12, maxMismatch = 1,
  #                            verbose=TRUE, returnRejects = TRUE, 
  #                            propagateCol=propagateCol))
  
  mergers[[sam]] <- merger
}

library(tidyverse)

write_rds(mergers, file = "mergers.rds")

rm(derepF); rm(derepR)

# names(mergers) <- paste("G", names(mergers), sep="")

# Construct sequence table  -----
# Load first seqtab
# seqtab <- makeSequenceTable(mergers[[1]])
# 
# # Add other seqtabs
# for(sam in 2:length(names(mergers))){
#   add_seqtab <- makeSequenceTable(mergers[[sam]])
#   seqtab <- mergeSequenceTables(seqtab, add_seqtab)
# }

seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, "seqtab.rds")

# and remove chimeras ----

# We filter high-quality reads by requiring an exact match to the proximal primer and the presence of the distal primer.


minOverlap <- 12

seqtab <- collapseNoMismatch(seqtab, 
                             minOverlap = minOverlap, 
                             orderBy = "abundance", verbose=TRUE)

seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "consensus", verbose = T)


# insilico amplicon trimming of ASVs
targetLength <- seq(100,180)
seqtab.nochim.targetLength <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% targetLength]


# save ----
# revisar nombres y renombrar, despues procesar combine_features e incluir fasta para filtrar!

save_seqtab <- function(seqtab) {
  # ============
  # Save Fasta
  # ============ 
  
  asv_seqs <- colnames(seqtab)
  n_asv <- dim(seqtab)[2]
  asv_headers <- vector(n_asv, mode="character")
  
  for (i in 1:n_asv) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
  }
  
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  
  write(asv_fasta,
        file=paste0(outpath, "/",
                    out_prefix, "_ASVs.fasta"))
  
  # ============
  # Save count table: 
  # ============
  
  
  
  
  asv_tab <- t(seqtab)
  
  id_samples <- colnames(asv_tab)
  camp <- sapply(strsplit(id_samples, "[-]"), `[`, 2)
  ship <- sapply(strsplit(id_samples, "[-]"), `[`, 3)
  
  row.names(asv_tab) <- sub(">", "", asv_headers)
  colnames(asv_tab) <- paste(ship, camp, sep = '.')
  
  write.table(asv_tab,
              file= paste0(outpath, "/",
                           out_prefix, "_ASVs_count.table"),
              sep="\t",
              row.names = TRUE,
              col.names = TRUE,
              quote=FALSE
  )
}

save_seqtab(seqtab.nochim.targetLength)


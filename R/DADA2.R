# Running data based on 388 libarries
# Ricardo Gomez-Reyes
# https://github.com/RJEGR/Small-RNASeq-data-analysis/blob/master/DOWNSTREAM_BKP/DIFFEXP/DADA2_iso_miRNA_LAB.R
# https://github.com/RJEGR/dada2/blob/de2ddd8a2f67607b694e81403f3f2d1c9d6126a3/multirun_18S.R

# Step 2 - Trimming (amplicon-) primers from both ends of paired-end reads


# ./trim_galore -q 10 --paired X7-TS_TCTGTCGT-TAGTTGCG_L002_R1_001.fastq.gz X7-TS_TCTGTCGT-TAGTTGCG_L002_R2_001.fastq.gz./trim_galore -q 10 --paired X7-TS_TCTGTCGT-TAGTTGCG_L002_R1_001.fastq.gz X7-TS_TCTGTCGT-TAGTTGCG_L002_R2_001.fastq.gz
# cutadapt -a TTGTACACACCGCCC...GTAGGTGAACCTGCRGAAGG -A CCTTCYGCAGGTTCACCTAC...GGGCGGTGTGTACAA --discard-untrimmed -o Nombre_archivo_F -p Nombre_archivo_R input_file_F input_file_R
# for f in $(ls *R1* | grep fastq); do ./trim_galore -q 10 --paired ${f%_R*}_R1_001.fastq.gz ${f%_R*}_R2_001.fastq.gz; done

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

packageVersion("dada2")

# devtools::install_github("benjjneb/dada2", ref="v1.26")

# remove.packages("dada2")

# install.packages("E:/Users/Israel V/Downloads/dada2-1.4.1.zip",
#                 repos = NULL,
#                 type = "source",
#                 dependencies = c("Depends", "Suggests","Imports"))

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

path <- "C:/Users/Israel V/Documents/MEIOFAUNA//raw-seqs/"

# path <- "D:/raw-seqs/"

setwd(path)

out_path <- file.path(path, run)

# create it manually if widows

# system(command = paste0("mkdir -p ", out_path), intern = F)


# outdir

system(paste0('mkdir -p ','/outputs'))

outpath <- paste0(path, "/outputs")

# Input dataset

str(fastq.gz <- sort(list.files(path, pattern=".gz$")))

str(fqFs <- sort(fastq.gz[grep('_R1_', fastq.gz)]))
str(fqRs <- sort(fastq.gz[grep('_R2_', fastq.gz)]))

# sampleo <- c(1, 85, 100);fqFs <- fqFs[sampleo];fqRs <- fqRs[sampleo]


# source("https://raw.githubusercontent.com/RJEGR/metagenomics/master/plotQualityProfile.R")

# plotQP(c(fqFs[1], fqRs[1]))

if(length(fqFs) != length(fqRs)) stop("Forward and reverse files do not match.")

# File parsing

filtpathF <- file.path(path, "filtered")
filtpathR <- file.path(path, "filtered")

filt <- sort(list.files(path, pattern="gz"))

filtFs <- sort(filt[grep('R1', filt)])
filtRs <- sort(filt[grep('R2', filt)])

# filtFs <- filtFs[sampleo]
# filtRs <- filtRs[sampleo]

str(snF <- sapply(strsplit(filtFs, "[_]"), `[`, 1))
str(snR <- sapply(strsplit(filtRs, "[_]"), `[`, 1))

snR[!snR %in% snF]

if(!identical(snF, snR)) stop("Forward and reverse files do not match.")

names(filtFs) <- snF
names(filtRs) <- snR

# metagroups

table(sapply(strsplit(filtFs, "-"), `[`, 1))


# Filtrar de manera homogenea todas las bibliotecas como estandar 

# truncLen <- c(150,140)

truncLen <- c(0,0)

filterAndTrim <- filterAndTrim(fwd=file.path(path, fqFs), filt=file.path(filtpathF, fqFs),
                               rev=file.path(path, fqRs), filt.rev=file.path(filtpathR, fqRs),
                               truncLen = truncLen, 
                               # trimRigh = c(100,100), 
                               maxEE=1, truncQ=2, maxN=0, rm.phix=TRUE,
                               compress=TRUE, verbose=TRUE, multithread = F)

file_path_filt <- paste0(filtpathF, "/", out_prefix, "_filterAndTrim.tsv")

head(filterAndTrim <- as_tibble(filterAndTrim, rownames = "Sample"))

filterAndTrim <- filterAndTrim %>% mutate(Sample = sapply(strsplit(Sample, "[_]"), `[`, 1))

write.table(filterAndTrim, 
            file = file_path_filt, 
            sep="\t", row.names = F, col.names = T)

# sample.groups <- sapply(strsplit(filtFs, "-"), `[`, 1)

# filtFs <- c(file.path(filtpathF, fqFs)[1], file.path(filtpathR, fqRs)[1])
# plotQP(filtFs)
# 
# read.delim(paste0(filtpathF, "/", out_prefix, "_filterAndTrim.tsv"))

pct_trim <- filterAndTrim$reads.out/filterAndTrim$reads.in

hist(pct_trim)

# Learn Error ----
# MAX_CONSIST <- 40
# 

threads <- F

# 
# errF <- learn_Err(filtFs)
# 
# plotErrors(errF, nominalQ=TRUE)



# Infer Sequence Variants

sample.groups <- unique(sapply(strsplit(fqFs, "[-]"), `[`, 1))
mergers <- vector("list", length(sample.groups))
names(mergers) <- sample.groups


# Learn error individually (by batch) and , derep, then run dada 
# upgrade filtered-pass files

filt <- sort(list.files(filtpathF, pattern="gz"))

str(fqFs <- sort(filt[grep('R1', filt)]))
str(fqRs <- sort(filt[grep('R2', filt)]))


table(sapply(strsplit(fqFs, "-"), `[`, 1))


for(sam in sample.groups) {
  
  run <- sapply(strsplit(fqFs, "[-]"), `[`, 1)
  
  which_samples <- which(run %in% sam)
  
  cat("Processing batch:", sam, "\n")
  
  Fs <- file.path(filtpathF, fqFs)[which_samples]
  Rs <- file.path(filtpathR, fqRs)[which_samples]
  
  if(length(Fs) != length(Rs)) stop("Forward and reverse files do not match.")
  
  samples <- sapply(strsplit(basename(Fs), "[_]"), `[`, 1)
  
  cat("Processing Files:",  samples, "\n")
  
  cat("learn Errors in Run: ", sam, "\n")
  
  errF <- learnErrors(Fs, nbases=1e8, multithread = F, randomize = T)
  errR <- learnErrors(Rs, nbases=1e8, multithread = F, randomize = T)
  
  cat("Dereplicate files: ", samples, "\n")
  
  derepF <- derepFastq(Fs)
  derepR <- derepFastq(Rs)
  
  cat("Infering amplicon sequence variance using the error rate learned from run:", sam, "\n")
  
  ddF <- dada(derepF, err=errF, multithread = F)
  
  ddR <- dada(derepR, err=errR, multithread = F)
  
  # propagateCol <- names(ddR)
  
  nrow(merger <- mergePairs(ddF, derepF, ddR, derepR
                             # minOverlap = 10, maxMismatch = 1,
                             # verbose=TRUE, 
                             # propagateCol=propagateCol
                            ))
  
  mergers[[sam]] <- merger
}


write_rds(list(ddF, ddR), file = "dadaError.rds")

class(mergers)

write_rds(mergers, file = "mergers.rds")

library(tidyverse)

# mergers <- read_rds("mergers.rds")




rm(derepF); rm(derepR)

# names(mergers) <- paste("G", names(mergers), sep="")

# Construct sequence table  -----
# this is for big-data processing step
# Load first seqtab

str(seqtab <- makeSequenceTable(mergers[[1]]))

# # Add other seqtabs

sample.groups <- names(mergers)
add_seqtab <- vector("list", length(sample.groups))
names(add_seqtab) <- sample.groups


for(sam in 2:length(names(mergers))){
  cat("makeSequenceTable :", names(mergers)[sam], "\n")
  
  add_seqtab <- makeSequenceTable(mergers[[sam]])

  seqtab <- mergeSequenceTables(seqtab, add_seqtab)
}

# upgrade

# 
# 
# for(sam in 1:length(names(mergers))){
#   
#   cat("makeSequenceTable :", names(mergers)[sam], "\n")
#   
#   add_seqtab[[sam]] <- makeSequenceTable(mergers[[sam]])
#   
#   # seqtab <- mergeSequenceTables(add_seqtab)
# }
# 
# seqtab.all <- dada2::mergeSequenceTables(tables = add_seqtab)
# 
# sqtb.all <- dada2::mergeSequenceTables(table1 = add_seqtab[[1]], table2 = add_seqtab[[2]])

write_rds(seqtab, file = "seqtab.rds")

# and remove chimeras ----

# We filter high-quality reads by requiring an exact match to the proximal primer and the presence of the distal primer.
dim(seqtab)


minOverlap <- 10

seqtab.collapsed <- collapseNoMismatch(seqtab, 
                             minOverlap = minOverlap, 
                             orderBy = "abundance", verbose=TRUE)

# Output 10160 collapsed sequences out of 11276 input sequences.

write_rds(seqtab.collapsed, file = "seqtab.collapsed.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab.collapsed, method = "consensus", verbose = T)

# Identified 491 bimeras out of 10160 input sequences.

dim(seqtab.nochim)

write_rds(seqtab.nochim, file = "seqtab.nochim.rds")


# insilico amplicon trimming of ASVs

targetLength <- seq(100,200)

seqtab.nochim.targetLength <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% targetLength]

dim(seqtab.nochim.targetLength)

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
        file=paste0(out_prefix, "_ASVs.fasta"))
  
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
              file= paste0(out_prefix, "_ASVs_count.table"),
              sep="\t",
              row.names = TRUE,
              col.names = TRUE,
              quote=FALSE
  )
}

save_seqtab(seqtab.nochim)


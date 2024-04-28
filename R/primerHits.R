
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


# From https://benjjneb.github.io/dada2/ITS_workflow.html

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
    RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}


primerHits <- function(primer, fn) {
  require(ShortRead)
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


# test primer detection

raw_path <- "~/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/"

# clean_path <- "~/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/filtered/"

path <- raw_path
  
setwd(path)

str(fastq.gz <- sort(list.files(path, pattern=".fastq.gz$", full.names = T)))

str(fqFs <- sort(fastq.gz[grep('_R1_', fastq.gz)]))
str(fqRs <- sort(fastq.gz[grep('_R2_', fastq.gz)]))

# lapply(fqFs, function(x) system(paste0("md5sum ", x)))

# system(paste0("md5sum ", fqFs))
# system(paste0("md5sum ", fqRs))

# About the amplicon 

# A universal-specific forward/eukaryotic-specific reverse primer combination: 

# 1389F/1510R (Amaral-Zettler et al 2009)

# Roche adapters (shown after ^) and the “barcodes” or 5-base keys (shown as X's) for distinguishing between samples 


# 18S V9 - Euk_1389f (Universal):   AmpliconPCR Forward Primer
# 5p GCCTCCCTCGCGCCATCAGXXXXXTTGTACACACCGCCC 3p
# 5p GCCTCCCTCGCGCCATCAGXXXXX^ p-f-primer

# 18S V9 - Euk_1510r: Amplicon PCR Reverse Primer
# 5p GCCTTGCCAGCCCGCTCAGCCTTCYGCAGGTTCACCTAC 3p
# 5p GCCTTGCCAGCCCGCTCAG^ p-r-primer

# Illumina primers designed for 

# 5p TTGTACACACCGCCC--- 3p Euk_1389f (Amaral-Zettler et al 2009)
# 5p --GTACACACCGCCCGTC	3p LinkerPrimerSequence from illumina

# 5p ----CCTTCYGCAGGTTCACCTAC 3p Euk_1510r (Amaral-Zettler et al 2009)
# 5p TGATCCTTCTGCAGGTTCACCTAC	3p ReversePrimer from illumina

# Easy way to draw adapter/primer 
# zgrep --color=always "AGATCGGAAGAGC" *.fastq.gz

# FWD <- "TTGTACACACCGCCC"  
# REV <- "CCTTCYGCAGGTTCACCTAC" 

FWD <- "GTACACACCGCCCGTC"  
REV <- "TGATCCTTCTGCAGGTTCACCTAC" 

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


# 
# Ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. 
# Because The percentage of base calls at each position for which an N was called is > 0

fnFs.filtN <- file.path(path, "filtN", basename(fqFs))
fnRs.filtN <- file.path(path, "filtN", basename(fqRs))

# filterAndTrim(fqFs, fnFs.filtN, fqRs, fnRs.filtN, maxN = 0, truncQ = 0, compress = TRUE, verbose=TRUE)

fqFs_hit <- fnFs.filtN[[10]]
fqRs_hit <- fnRs.filtN[[10]]


rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fqFs_hit), 
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fqRs_hit), 
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fqFs_hit), 
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fqRs_hit))


# Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0    4882
# REV.ForwardReads       0          0       0    5889
# REV.ReverseReads       1          0       0       0


# remove primers

path.cut <- file.path(path, "filtN/cutadapt")

if(!dir.exists(path.cut)) dir.create(path.cut)


fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt

finding_adapter_params <- c("-e", 0.1, "-O", 1, "-q", 10)

for(i in seq_along(fqFs)) {
  system2(command = "/Users/cigom/opt/anaconda3/bin/cutadapt", 
    args = c("-j", 2, finding_adapter_params,
      R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
      "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
      fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# test primer zero

fqFs_hit <- fnFs.cut[[10]]
fqRs_hit <- fnRs.cut[[10]]
# 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fqFs_hit),
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fqRs_hit),
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fqFs_hit),
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fqRs_hit))


# Additional primer trimming
# (Illumina TruSeq, Sanger iPCR; auto-detected w/ trimgalone)

TruSeq <- "AGATCGGAAGAGC" 
TruSeq.orients <- allOrients(TruSeq)

fqFs_hit <- fnFs.cut[[10]]
fqRs_hit <- fnRs.cut[[10]]

rbind(
  TruSeq.ForwardReads = sapply(TruSeq.orients, primerHits, fn = fqFs_hit),
  TruSeq.ReverseReads = sapply(TruSeq.orients, primerHits, fn = fqRs_hit))


path.Illum <- file.path(path, "filtN/cutadapt/Illumina")

if(!dir.exists(path.Illum)) dir.create(path.Illum)


fnFs.Illum <- file.path(path.Illum, basename(fnFs.filtN))
fnRs.Illum <- file.path(path.Illum, basename(fnRs.filtN))

TruSeq.RC <- dada2:::rc(TruSeq)

R1.flags <- paste("-g", TruSeq, "-a", TruSeq.RC) 

R2.flags <- paste("-G", TruSeq, "-A", TruSeq.RC) 

for(i in seq_along(fqFs)) {
  system2(command = "/Users/cigom/opt/anaconda3/bin/cutadapt", 
    args = c("-j", 2, finding_adapter_params,
      R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
      "-o", fnFs.Illum[i], "-p", fnRs.Illum[i], # output files
      fnFs.cut[i], fnRs.cut[i])) # input files
}



fqFs_hit <- fnFs.Illum[[10]]
fqRs_hit <- fnRs.Illum[[10]]

rbind(
  TruSeq.ForwardReads = sapply(TruSeq.orients, primerHits, fn = fqFs_hit),
  TruSeq.ReverseReads = sapply(TruSeq.orients, primerHits, fn = fqRs_hit))


# datviz
# system2(command = "mkdir", args = c("-f", "fastqc"))

# system2(command = "fastqc", args = c("*.fastq.gz","-o", "fastqc"))

# system2(command = "multiqc", args = c("fastqc/*zip","-o", "multiqc"))


# Run dada2denoisePaired.R




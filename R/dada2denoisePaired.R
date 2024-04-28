
# After primerHit.R was done
# Run dada2denoisePaired.R

# Filtrar de manera homogenea todas las bibliotecas como estandar 
# esto permite hacer el merge de cada lote mas adelante.


work_path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp//filtN/cutadapt/Illumina"

fastq.filt <- list.files(work_path, pattern = ".fastq.gz" ,full.names = T)

str(fnFs.Illum <- sort(fastq.filt[grep('_R1_', fastq.filt)]))
str(fnRs.Illum <- sort(fastq.filt[grep('_R2_', fastq.filt)]))

fwd <- fnFs.Illum #[1:5] # fnFs.cut
rev <- fnRs.Illum # [1:5] # fnRs.cut


filtered_path <- file.path(dirname(fwd[1]), "filterAndTrim") # file.path(path, "cutadapt")

filt.fwd <- file.path(filtered_path, basename(fwd))

filt.rev <- file.path(filtered_path, basename(rev))

filterAndTrim <- filterAndTrim(
  fwd = fwd, filt = filt.fwd,
  rev= rev, filt.rev = filt.rev,
  maxEE = Inf, truncQ = 2, maxN = 0, rm.phix=TRUE,
  compress = TRUE, verbose=TRUE, multithread = F)

# Infer Sequence Variants

DenoisePaired <- function(filtFs, filtRs) {
  
  if(length(filtFs) != length(filtRs)) stop("Forward and reverse files do not match.")
  
  S <- filtFs
  
  B <- unique(sapply(strsplit(basename(S), "[-]"), `[`, 1))
  
  cat("learn Errors in the batch: ", B, "\n")
  
  cat("Using : ", length(S), " samples \n")
  
  errF <- learnErrors(filtFs, nbases=1e8, multithread = T, randomize = T)
  
  errR <- learnErrors(filtRs, nbases=1e8, multithread = T, randomize = T)
  
  ERR <- list("errF" = errF, "errR" = errR)
  
  readr::write_rds(ERR, file = file.path(dirname(S[1]), paste0(B, "_learnErrors", ".RData")))
  
  cat("Dereplicate files: ", S, "\n")
  
  derepF <- derepFastq(filtFs, verbose = T)
  
  derepR <- derepFastq(filtRs, verbose = T)
  
  cat("Infering amplicon sequence variance using the error rate learned in this batch\n")
  
  dadaFs <- dada(derepF, err=errF, multithread = TRUE)
  
  dadaRs <- dada(derepR, err=errR, multithread = TRUE)
  
  
  merger <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 12, maxMismatch = 0) 
  
  return(merger)
  
  # seqtab <- makeSequenceTable(mergers)
  # 
  # seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
}

fastq.filt <- list.files(filtered_path, pattern = ".fastq.gz" ,full.names = T)

str(filt.fwd <- sort(fastq.filt[grep('_R1_', fastq.filt)]))
str(filt.rev <- sort(fastq.filt[grep('_R2_', fastq.filt)]))

filt.fwd[!basename(filt.fwd) %in% basename(fwd)]

# merger_out <- DenoisePaired(filt.fwd, filt.rev)
# str(seqtab <- makeSequenceTable(merger_out))


sample.groups <- unique(sapply(strsplit(basename(filt.fwd), "[-]"), `[`, 1))
mergers <- vector("list", length(sample.groups))
names(mergers) <- sample.groups


batch.f <- split(filt.fwd, sapply(strsplit(basename(filt.fwd), "[-]"), `[`, 1))
batch.r <- split(filt.rev, sapply(strsplit(basename(filt.rev), "[-]"), `[`, 1))


t(lapply(batch.f, length))
t(lapply(batch.r, length))

for(i in seq_along(batch.f)) {
  
  mergers[[i]] <- DenoisePaired(batch.f[[i]], batch.r[[i]])
  
}

# merger_out <- DenoisePaired(batch.f$X1, batch.r$X1)
# str(seqtab <- makeSequenceTable(merger_out))


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


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


dim(seqtab.nochim)

table(nchar(getSequences(seqtab.nochim)))

targetLength <- seq(90,150)

seqtab.nochim.target <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% targetLength]

table(nchar(getSequences(seqtab.nochim.target)))

hist(nchar(getSequences(seqtab.nochim.target)))

readr::write_rds()

outrds <- list("filterAndTrim" = filterAndTrim, 
  "mergers" = mergers, "seqtab" = seqtab, 
  "seqtab.nochim" = seqtab.nochim, 
  "seqtab.nochim.target" = seqtab.nochim.target)

readr::write_rds(outrds, file = file.path(filtered_path, paste0("dada2denoisedPaired_outputs", ".RData")))



# Track

getN <- function(x) sum(getUniques(x))

# sapply(mergers[[1]], getN) %>% dplyr::bind_rows() %>% pivot_longer(names(.))

merge_out <- vector("list", length(mergers))

names(merge_out) <- names(mergers)

for (i in 1:length(mergers)) {
  merge_out[[i]] <- sapply(mergers[[i]], getN) %>% dplyr::bind_rows() %>% 
    pivot_longer(names(.), names_to = "Sample", values_to = "merged") %>%
    mutate(Sample = sapply(strsplit(Sample, "[_]"), `[`, 1))
  
}

do.call(bind_rows, merge_out) -> merge_out 

track <- data.frame(tabled = rowSums(seqtab), nonchim = rowSums(seqtab.nochim), targetLen = rowSums(seqtab.nochim.target)) %>% 
  dplyr::as_tibble(rownames = "Sample") %>%
  mutate(Sample = sapply(strsplit(Sample, "[_]"), `[`, 1)) 

view(track)

save_seqtab <- function(seqtab, out_path = getwd()) {
  
  # ============
  # Save Fasta
  # ============ 
  
  require(digest)
  
  
  asv_seqs <- dada2::getSequences(seqtab)
  
  asv_headers <- sapply(asv_seqs, digest, algo="md5")
  
  asv_headers <- paste(">", as.vector(asv_headers), sep = "")
  
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  
  write(asv_fasta, file = file.path(out_path, "dna_sequences.fa"))
  
  # ============
  # Save count table: 
  # ============
  
  # Include a sanity check for md5 ids
  
  asv_tab <- t(seqtab)
  
  asv_tab <- data.frame(asv_tab)
  
  if(!identical(rownames(asv_tab), asv_seqs)) stop("md5 names do not match.")
  
  # id_samples <- colnames(asv_tab)
  
  LIBRARY_ID <- sapply(strsplit(colnames(asv_tab), "_"), `[`, 1)
  LIBRARY_ID <- sub("[.]", "-", LIBRARY_ID)
  
  row.names(asv_tab) <- sub(">", "", asv_headers)
  
  colnames(asv_tab) <-  LIBRARY_ID
  
  # write_tsv(data.frame(asv_tab), file= paste0(out_prefix, "_ASVs_count.table"))
  
  write.table(asv_tab,
    file= file.path(out_path, "feature_table.tsv"),
    sep="\t",
    row.names = TRUE,
    col.names = TRUE,
    quote=FALSE
  )
}


save_seqtab(seqtab.nochim.target, out_path = filtered_path)

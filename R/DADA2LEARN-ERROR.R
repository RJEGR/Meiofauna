# PLOT individual error Learned from batch
# as https://github.com/RJEGR/dada2/blob/de2ddd8a2f67607b694e81403f3f2d1c9d6126a3/learnError_dada_test.R#L4

# write_rds(list(ddF, ddR), file = "dadaError.rds")


path <- "C:/Users/Israel V/Documents/MEIOFAUNA//raw-seqs/"

setwd(path)

fastq.gz <- sort(list.files(path, pattern="fq.gz"))


fqFs <- sort(fastq.gz[grep('R1', fastq.gz)])
fqRs <- sort(fastq.gz[grep('R2', fastq.gz)])


library(tidyverse)

# Plot error ----
# plotErrors()

transdf_q <- function (dq, 
                       nti = c("A", "C", "G", "T"), 
                       ntj = c("A", "C", "G", "T"), 
                       obs = TRUE, err_out = TRUE, 
                       err_in = FALSE, nominalQ = FALSE) 
{
  
  require(dada2);require(reshape2)
  
  ACGT <- c("A", "C", "G", "T")
  
  if (!(all(nti %in% ACGT) && all(ntj %in% ACGT)) || any(duplicated(nti)) || 
      any(duplicated(ntj))) {
    stop("nti and ntj must be nucleotide(s): A/C/G/T.")
  }
  dq <- getErrors(dq, detailed = TRUE, enforce = FALSE)
  if (!is.null(dq$trans)) {
    if (ncol(dq$trans) <= 1) {
      stop("plotErrors only supported when using quality scores in the error model (i.e. USE_QUALS=TRUE).")
    }
    transdf = melt(dq$trans, factorsAsStrings = TRUE)
    colnames(transdf) <- c("Transition", "Qual", "count")
  }
  else if (!is.null(dq$err_out)) {
    if (ncol(dq$err_out) <= 1) {
      stop("plotErrors only supported when using quality scores in the error model (i.e. USE_QUALS=TRUE).")
    }
    transdf = melt(dq$err_out, factorsAsStrings = TRUE)
    colnames(transdf) <- c("Transition", "Qual", "est")
  }
  else {
    stop("Non-null observed and/or estimated error rates (dq$trans or dq$err_out) must be provided.")
  }
  transdf$from <- substr(transdf$Transition, 1, 1)
  transdf$to <- substr(transdf$Transition, 3, 3)
  if (!is.null(dq$trans)) {
    tot.count <- tapply(transdf$count, list(transdf$from, 
                                            transdf$Qual), sum)
    transdf$tot <- mapply(function(x, y) tot.count[x, y], 
                          transdf$from, as.character(transdf$Qual))
    transdf$Observed <- transdf$count/transdf$tot
  }
  else {
    transdf$Observed <- NA
    obs <- FALSE
  }
  if (!is.null(dq$err_out)) {
    transdf$Estimated <- mapply(function(x, y) dq$err_out[x, 
                                                          y], transdf$Transition, as.character(transdf$Qual))
  }
  else {
    transdf$Estimated <- NA
    err_out <- FALSE
  }
  if (!is.null(dq$err_in)) {
    ei <- dq$err_in
    if (is.list(ei)) 
      ei <- ei[[1]]
    transdf$Input <- mapply(function(x, y) ei[x, y], transdf$Transition, 
                            as.character(transdf$Qual))
  }
  else {
    transdf$Input <- NA
    err_in <- FALSE
  }
  transdf$Nominal <- (1/3) * 10^-(transdf$Qual/10)
  transdf$Nominal[transdf$Transition %in% c("A2A", "C2C", 
                                            "G2G", "T2T")] <- 1 - 10^-(transdf$Qual[transdf$Transition %in% 
                                                                                      c("A2A", "C2C", "G2G", "T2T")]/10)
  
  
  return(transdf)
}


nti = c("A", "C", "G", "T")
ntj = c("A", "C", "G", "T")

merge_errFs <- lapply(errFs, transdf_q) %>% do.call(rbind, .)
merge_errRs <- lapply(errRs, transdf_q) %>% do.call(rbind, .)

sample.groups <- unique(sapply(strsplit(fqFs, "[-]"), `[`, 1))
mergers <- vector("list", length(sample.groups))
names(mergers) <- sample.groups

for (i in 1:n_asv) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

transdf_q(ddF[[1]]$trans)

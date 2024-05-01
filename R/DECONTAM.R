

# A defined set of “negative control” samples in which sequencing was performed on blanks without any biological sample added. Extraction controls are preferred, and in amplicon sequencing the negative controls should also be carried through the PCR step, as each step in the workflow has the potential to introduce new contaminants.

# Therefore, omit this step
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/"

library(phyloseq)
library(tidyverse)

ps <- read_rds(paste0(wd, '/phyloseq.rds'))

library(decontam)

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

df$Sample_or_Control <- ifelse(df$Cruise %in% "XIXIMI-7", "Ctrl", "Sampl")
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

# In our phyloseq object, "quant_reading" is the sample variable that holds the concentration information:

contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")

head(contamdf.freq)

# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

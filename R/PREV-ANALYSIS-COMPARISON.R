
# LOAD ABUNDANCE AND SUM READS PER BATCH

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/Documents/MEIOFAUNA_PAPER/INPUTS/"

subdir <- paste0(wd, "CURATED_DB_DIR")
# library(Biostrings)
# library(rstatix)

library(tidyverse)

# tax_f <- read_tsv(list.files(path = subdir, pattern = 'taxonomy.tsv$', full.names = T))

ab_f <-  read_tsv(list.files(path = wd, pattern = 'table_100_80', full.names = T), skip = 1)

# COUNT SAMPLES

table(sapply(strsplit(names(ab_f)[-1], "-"), `[`, 1))

sum(table(sapply(strsplit(names(ab_f)[-1], "-"), `[`, 1)))


ab_f %>% pivot_longer(-`Feature ID`) %>%
  mutate(name = sapply(strsplit(name, "-"), `[`, 1)) %>%
  count(name)

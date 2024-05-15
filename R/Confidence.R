
# Contrast confidence values

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/CLASSIFICATION/classify-sklearn_dir/"

# 1

tax <- list.files(path = wd, pattern = "taxonomy.tsv", full.names = T)

tax <- read_tsv(tax) %>% dplyr::rename('DBRef' = 'Confidence')

# 2


SSURef <- read_tsv(list.files(path = "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/CLASSIFICATION/AGAINST_DB/SSURef_dir", pattern = "taxonomy.tsv", full.names = T)) %>%  
  dplyr::rename('SSURef' = 'Confidence')

# 3

PR2Ref <- read_tsv(list.files(path = "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/CLASSIFICATION/AGAINST_DB/pr2_dir/", pattern = "taxonomy.tsv", full.names = T)) %>%
  dplyr::rename('PR2Ref' = 'Confidence')


df <- cbind(tax, SSURef = SSURef$SSURef, PR2Ref = PR2Ref$PR2Ref)

head(df)

library(rstatix)

rstatix::cor_mat(df, DBRef,SSURef,PR2Ref)


df <- df %>% pivot_longer(cols = c("DBRef", "SSURef", "PR2Ref"), values_to = "Confidence", names_to = "Dataset")

df %>% 
  mutate(Dataset = factor(Dataset, levels = c("DBRef", "SSURef", "PR2Ref"))) %>%
  ggplot(aes(Confidence, color = Dataset)) + 
  geom_histogram(aes(fill = Dataset)) + facet_grid(~ Dataset) +
  # ggplot2::stat_ecdf(size = 1.2) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  scale_fill_grey() +
  scale_color_grey()

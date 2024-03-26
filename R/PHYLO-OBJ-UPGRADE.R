# READ TAX AND PROCESS
# READ COUNT AND MERGE TO TAX
# PREPARE PHYLO OBJECT

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# wd <- "C:/Users/Israel V/Documents/MEIOFAUNA/raw-seqs/CLASSIFIER-STEP/sklearn-classify_dir/"

# wd <- "~/Downloads/RESULTS/classify-consensus-blast_dir/"
wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/MULTIRUN-ANALYSIS/LIBRARY/classify-consensus-blast_dir"

# 1

tax <- list.files(path = wd, pattern = "taxonomy.tsv", full.names = T)

tax <- read_tsv(tax)

# hist(tax$Confidence)
hist(tax$Consensus)

into <- c("k", "p", "c", "o","f", "g", "s")

mutate_ranks <- c("none" = NA_character_, "os__" = "", "[a-z]__" = "", 
                                            "Incertae_Sedis" = NA_character_, 
                                            "uncultured" = NA_character_, 
                                            "Unassigned"=NA_character_)
tax <- tax %>% 
  separate(Taxon, sep = ";", into = into) %>%
  mutate_all(list(~ str_replace_all(., mutate_ranks))) %>%
  mutate_all(function(x) {na_if(x,"")}) %>%
  data.frame(row.names = .$`Feature ID`)
  


# 2) BIND W/ ABUNDANCE

# wd <- "C:/Users/Israel V/Documents/MEIOFAUNA/raw-seqs/"

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/MULTIRUN-ANALYSIS/LIBRARY/"


ab_f <- list.files(path = wd, pattern = 'multirun_ASVs_count.table', full.names = T)

ab <- read.delim(ab_f, sep = "\t") # %>% as_tibble(rownames = "Feature ID")

names(ab) %>% as_tibble()

DB <- read_tsv(list.files(path = wd, pattern = "mapping-file-corregido.tsv", full.names = T)) %>%
  dplyr::rename("LIBRARY_ID" = "#SampleID") %>% mutate(LIBRARY_ID = gsub("^XI", "X", LIBRARY_ID))

all_sam <- gsub("[.]", "-", names(ab))



all_sam[!all_sam %in% DB$LIBRARY_ID]

data.frame("LIBRARY_ID" = all_sam  ) %>% as_tibble() %>%
  anti_join(DB, by = "LIBRARY_ID") 



MTD <-  read_tsv(list.files(path = wd, pattern = "mapping-file-corregido.tsv", full.names = T)) %>% 
  select(`#SampleID`, Depth, Region) %>% 
  mutate(Region = factor(Region, levels = c("Yucatan", "NW Shelf", "NW Slope", "Deep-sea"))) %>%
  dplyr::rename("LIBRARY_ID" = "#SampleID") %>%
  mutate(LIBRARY_ID = gsub("-", ".", LIBRARY_ID)) %>%
  filter(LIBRARY_ID %in% names(ab)) %>%
  data.frame(row.names = .$LIBRARY_ID)


identical(names(ab),rownames(MTD))

identical(rownames(ab), rownames(tax))

library(phyloseq)

phyloseq = phyloseq(otu_table(ab, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')) 
                    # sample_data(MTD)
)

saveRDS(phyloseq, paste0(wd, "phyloseq.rds"))


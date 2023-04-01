# DATABASE PROCESSING

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)



library(Biostrings)
library(tidyverse)

# 1) SSURef_NR99-138.1 ----

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/SSURef_NR99-138.1/"

subpath <- paste0(path,"/SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r_dir")

tax_f <- list.files(path = subpath, pattern = "taxonomy.tsv", full.names = T)
dna_f <- list.files(path = subpath, pattern = "dna-sequences.fasta", full.names = T)

dna <- Biostrings::readDNAStringSet(dna_f) 

length(dna) # From SSURef_NR99-138.1: ~ 55, 348 V9 reads

str(IDS <- names(dna))

IDS <- sapply(strsplit(IDS, " "), `[`, 1)

# Rename DNA names

names(dna) <- IDS

length(unique(IDS)) # Are unique?

head(IDS <- sort(IDS))

# 1.1) V9-SEQS-BASED TAXON FILTERING =====

nrow(tax <- read_tsv(tax_f, col_names = T) %>% filter(`Feature ID` %in% IDS)) # Must be 55,348 groups

# Sanity check

length(dna)==nrow(tax)

# 1.2) FILTER EUK AND CLEAN REDUNDANT TAXONOMY ----

table(sapply(strsplit(tax$Taxon, ";"), `[`, 1))

# d__Archaea: 1540  
# d__Bacteria: 31015 
# d__Eukaryota: 22,793 


into <- c("domain", "phylum", "class", "order", "family", "genus", "specie")

# into <- paste0("rank", 1:7)

tax_df <- tax %>%
  filter(grepl("Eukaryota", Taxon)) %>% # 22,794
  separate(Taxon, sep = ";", into = into) %>% view()
  mutate_at(into, funs(str_replace_all(., c("os__" = "", "[a-z]__" = "", 
    "_sp."=NA_character_, "metagenome"=NA_character_,
    "Incertae_Sedis" = NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_)))) %>%
  mutate_at(vars(all_of(into)), funs(gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(into),  ~na_if(., "")) %>%
  mutate_at(vars(into), funs(str_replace_all(., c("_" = " "))))

# 1.3) ADD SEQ. WIDTHS ----

lengths <- width(dna)

out <- as_tibble(structure(lengths, names = names(dna)), rownames = 'Feature ID') %>% 
  right_join(tax_df) %>% # The euk filter happend here
  dplyr::rename("dna_width"="value")


# 1.4) FILTER EUK SEQUENCES ----


EUK_IDS <- tax_df %>% pull(`Feature ID`) %>% sort()

sum(keep <- names(dna) %in% EUK_IDS) # 22,794

dna <- dna[keep]

# 1.5) OUTPUT -----

taxa_out <- "/SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed.rds"

write_rds(out, file = paste0(subpath, taxa_out))


seqs_out <- "/SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed.fasta"

Biostrings::writeXStringSet(dna, filepath =  paste0(subpath, seqs_out), format = "fasta")


# 2) PR2 V 4.14 SSU  ----
# (pr2_version_4.14.0_SSU-seqs-euk-derep-uniq-1389f-1510r.qza)

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/"

subpath <- paste0(path, "/pr2_version_4.14.0_SSU")

tax_f <- list.files(path = path, pattern = "taxonomy.tsv", full.names = T)
dna_f <- list.files(path = subpath, pattern = "dna-sequences.fasta", full.names = T)

dna <- Biostrings::readDNAStringSet(dna_f) # From pr2_version_4.14.0_SSU: 36,687 V9 reads

str(IDS <- names(dna))

IDS <- sapply(strsplit(IDS, " "), `[`, 1)

length(unique(IDS))

head(IDS <- sort(IDS))


nrow(tax <- read_tsv(tax_f, col_names = T) %>% filter(`Feature ID` %in% IDS)) # Must be 36,687 groups

# 2.2) V9-SEQS-BASED TAXON FILTERING ----

table(sapply(strsplit(tax$Taxon, ";"), `[`, 1))

into <- c("domain", "supergroup", "division", "class", "order", "family", "genus", "specie")

# into <- paste0("rank", 1:8)

tax_df <- tax %>%
  filter(grepl("Eukaryota", Taxon)) %>% # 
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(into, funs(str_replace_all(., c("os__" = "", "[a-z]__" = "", 
    "_sp."=NA_character_, "metagenome"=NA_character_,
    "Incertae_Sedis" = NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_, "[a-z]_X" = NA_character_)))) %>%
  mutate_at(vars(all_of(into)), funs(gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(into),  ~na_if(., "")) %>%
  mutate_at(vars(into), funs(str_replace_all(., c("_" = " "))))

# 2.3) (omit) STACK SILVA RANKS ----

into_silva <- c("domain", "phylum", "class", "order", "family", "genus", "specie")

all_cols <- c("Feature ID", into_silva)

# tax_df %>%  select_at(all_of(all_cols))


# 2.3) ADD SEQ. WIDTHS ----


lengths <- width(dna)

out <- as_tibble(structure(lengths, names = names(dna)), rownames = 'Feature ID') %>% 
  right_join(tax_df) %>% # The euk filter happend here
  dplyr::rename("dna_width"="value")

# 2.4) FILTER EUK SEQUENCES ----


EUK_IDS <- tax_df %>% pull(`Feature ID`) %>% sort()

sum(keep <- names(dna) %in% EUK_IDS) # 36, 677

dna <- dna[keep]

# Sanity check

sum(sort(names(dna)) == EUK_IDS) # 36,677


# 2.5) OUTPUT ----


taxa_out <- "/pr2_version_4.14.0_SSU-1389f-1510r-db_processed.rds"

write_rds(out, file = paste0(subpath, taxa_out))

seqs_out <- "/pr2_version_4.14.0_SSU-1389f-1510r-db_processed.fasta"

Biostrings::writeXStringSet(dna, filepath =  paste0(subpath, seqs_out), format = "fasta")

# 3 GO TO DB_TAXA2WORMS_CONVERTION.R script ----
# 4 GO TO DB_EXPLORATORY_COMPARISON.R script


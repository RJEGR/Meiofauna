# DATABASE PROCESSING

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)



library(Biostrings)
library(tidyverse)

# 1) Recovery V9 SILVA TAXA based on fasta ----

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RESCRIPT"
subpath <- paste0(path, "/SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r")

tax_f <- list.files(path = subpath, pattern = "taxonomy.tsv", full.names = T)
dna_f <- list.files(path = subpath, pattern = "dna-sequences.fasta", full.names = T)

dna <- Biostrings::readDNAStringSet(dna_f) # From SSURef_NR99-138.1: 55, 374 V9 reads

str(IDS <- names(dna))

IDS <- sapply(strsplit(IDS, " "), `[`, 1)

# Rename DNA names

names(dna) <- IDS

length(unique(IDS))

head(IDS <- sort(IDS))

nrow(tax <- read_tsv(tax_f, col_names = T) %>% filter(`Feature ID` %in% IDS)) # Must be 55,374 groups


# 1.2) Split and clean lineage ----

table(sapply(strsplit(tax$Taxon, " "), `[`, 1))

# into <- c("domain", "kingdom", "phylum", "class", "order", "suborder", "family", "genus", "specie")

into <- paste0("rank", 1:9)

tax_df <- tax %>%
  filter(grepl("Eukaryota", Taxon)) %>% # 22,792
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(into, funs(str_replace_all(., c("os__" = "", "[a-z]__" = "", 
    "Incertae_Sedis" = NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_))))

# 1.3) Add widths ----


sum(keep <- names(dna) %in% IDS)

lengths <- width(dna[keep])

out <- as_tibble(structure(lengths, names = names(dna)), rownames = 'Feature ID') %>% 
  right_join(tax_df) %>%
  rename("value"="dna_width")

path_out <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON/"
file_out <- "SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r.rds"

write_rds(out, file = paste0(path_out, file_out))

# 1.4) Output fasta

eukIDS <- tax_df$`Feature ID`

sum(keep <- names(dna) %in% eukIDS)

dna <- dna[keep]

file_out <- "SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r.fasta"

Biostrings::writeXStringSet(dna, filepath =  paste0(path_out, file_out), format = "fasta")

# 2) Recovery V9 - PR2 TAXA based on fasta ----
# (pr2_version_4.14.0_SSU-seqs-euk-derep-uniq-1389f-1510r.qza)

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/"

subpath <- paste0(path, "/pr2_version_4.14.0_SSU")

tax_f <- list.files(path = subpath, pattern = "taxonomy.tsv", full.names = T)
dna_f <- list.files(path = subpath, pattern = "dna-sequences.fasta", full.names = T)

dna <- Biostrings::readDNAStringSet(dna_f) # From pr2_version_4.14.0_SSU: 36,687 V9 reads

str(IDS <- names(dna))

IDS <- sapply(strsplit(IDS, " "), `[`, 1)

length(unique(IDS))

head(IDS <- sort(IDS))

nrow(tax <- read_tsv(tax_f, col_names = T) %>% filter(`Feature ID` %in% IDS)) # Must be 36,687 groups

# 2.2) Split and clean lineage ----

table(sapply(strsplit(tax$Taxon, ";"), `[`, 1))

# into <- c("domain", "supergroup", "division", "class", "order", "family", "genus", "specie")

into <- paste0("rank", 1:8)

tax_df <- tax %>%
  filter(grepl("Eukaryota", Taxon)) %>% # 
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(into, funs(str_replace_all(., c("os__" = "", "[a-z]__" = "", 
    "Incertae_Sedis" = NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_, "[a-z]_X" = NA_character_))))

# 2.3) Add widths ----


sum(keep <- names(dna) %in% IDS) # 36,687

lengths <- width(dna[keep])

out <- as_tibble(structure(lengths, names = names(dna)), rownames = 'Feature ID') %>% 
  right_join(tax_df) %>%
  rename("value"="dna_width")

path_out <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON/"
file_out <- "pr2_version_4.14.0_SSU-1389f-1510r.rds"

write_rds(out, file = paste0(path_out, file_out))

# 1.4) Output fasta

eukIDS <- tax_df$`Feature ID`

sum(keep <- names(dna) %in% eukIDS)

dna <- dna[keep]

file_out <- "pr2_version_4.14.0_SSU-1389f-1510r.fasta"

Biostrings::writeXStringSet(dna, filepath =  paste0(path_out, file_out), format = "fasta")


# 3 GO TO DB_EXPLORATORY_COMPARISON

tax_df %>% view()

find_obs_tax(c(NON_EUK, PHYLA), tax_df) %>% view()


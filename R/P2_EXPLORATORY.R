# EL PROPOSITO DE ESTE SCRIPT ES INSPECCIONAR LA BASE DE DATOS DE RP2, PUES INCLUYE TAXONES DE SEDIMENTOS (V. 4.14) QUE SON DE INTERES DEL PROYECTO DE ALEX

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)
library(tidyverse)


NON_EUK <- c("Amoebozoa",
  "Archaeplastida",
  "Excavata",
  "Fungi",
  "Stramenopiles",
  "Alveolata",
  "Rhizaria")


PHYLA <- c(
  "Annelida",
  "Arthropoda",
  "Brachiapoda",
  "Bryozoa",
  "Chordata",
  "Cnidaria",
  "Echinodermata",
  "Gastrotricha",
  "Gnathostomulida",
  "Hemichordata",
  "Kinorhyncha",
  "Mollusca",
  "Nematoda",
  "Nemertea",
  "Platyhelminthes",
  "Porifera",
  "Rotifera",
  "Xenacoelomorpha")




# Check PR2- V9 - groups ----
# CREATE: pr2_version_4.14.0_SSU-seqs-euk-derep-uniq-1389f-1510r.qza

PR2 <- "~/Documents/MEIOFAUNA_PAPER/pr2_version_4.14.0_SSU/"

tax_f <- list.files(path = PR2, pattern = "taxonomy.tsv", full.names = T)
dna_f <- list.files(path = PR2, pattern = "dna-sequences.fasta", full.names = T)

dna <- Biostrings::readDNAStringSet(dna_f) # 36,687 V9 reads

str(v9_names <- names(dna))

head(v9_names <- sort(v9_names))

tax <- read_tsv(tax_f, col_names = F)

tax <- tax %>% filter(X1 %in% v9_names)

keep <- names(dna) %in% v9_names

lengths <- width(dna[keep])

tax <- as_tibble(structure(lengths, names = names(dna)), rownames = 'X1') %>% right_join(tax)

names(tax) <- c("ID", "dna_width", "Taxon")


into <- paste0("rank", 1:8)

tax_df <- tax %>%
  filter(grepl("Eukaryota", Taxon)) %>% # "Eukaryota$"
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(into, funs(str_replace_all(., c("[a-z]_X" = NA_character_))))

  # filter(grepl(NON_EUK, X2)) %>%

tax_df %>% 
  pivot_longer(cols = all_of(into), values_to = "Taxon") %>%
  # filter(name == "rank1") %>%
  drop_na(Taxon) %>%
  count(name, Taxon) %>% view()




# OTHER



# PR2 <- "~/Documents/MEIOFAUNA_PAPER/pr2_version_4.14.0_SSU/"
# 
# SILVA <- "~/Documents/MEIOFAUNA_PAPER/RESCRIPT/SSURef_NR99-138.1/"

PATH <- "~/Documents/MEIOFAUNA_PAPER/PR2_SILVA/"

tax_f <- list.files(path = PATH, pattern = "SSURef_NR99-138.1.pr2_version_4.14.0_SSU.tax$", full.names = T)
dna_f <- list.files(path = PATH, pattern = "SSURef_NR99-138.1.pr2_version_4.14.0_SSU.fasta", full.names = T)


tax <- read_tsv(tax_f) %>% arrange(`Feature ID`)

tax %>% distinct(`Feature ID`) # todos los ids son distintos (good!)

into <- c("domain", "kingdom", "phylum", "class", "order", "suborder", "family", "genus", "specie")

tax_df <- tax %>% 
  filter(grepl("Eukaryota", Taxon)) %>%
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(into, funs(str_replace_all(., c("os__"=NA_character_, "[a-z]__" = ""))))

# Parse orders
tax %>%  head() %>% view()
tax_df %>% view()

str(V9_ID <- names(Biostrings::readDNAStringSet(dna_f)))
V9_ID <- sapply(strsplit(V9_ID, " "), `[`, 1)
head(V9_ID <- sort(V9_ID))
tax %>% filter(`Feature ID` %in% V9_ID)

widths <- Biostrings::width(Biostrings::readDNAStringSet(dna_f))

# no es claro pq solo se recuperan 22,782 linajes de las 55,374 secuencias



tax_f <- list.files(path = PR2, pattern = "taxonomy.tsv", full.names = T)
dna_f <- list.files(path = PR2, pattern = "dna-sequences.fasta", full.names = T)


tax <- read_tsv(tax_f) %>% arrange(`Feature ID`) # 176,901

str(V9_ID <- names(Biostrings::readDNAStringSet(dna_f)))
V9_ID <- sapply(strsplit(V9_ID, " "), `[`, 1)
head(V9_ID <- sort(V9_ID))
tax %>% filter(`Feature ID` %in% V9_ID)

into <- c("domain", "supergroup", "division", "class", "order", "family", "genus", "specie")


tax_df <- tax %>% 
  filter(grepl("Eukaryota", Taxon)) %>%
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(into, funs(str_replace_all(., c("[a-z]__" = ""))))


tax %>% 
  filter(grepl("Alveolata", Taxon))

tax_df %>% distinct(subg)


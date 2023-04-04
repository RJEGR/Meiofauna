# DB_MERGING.R

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)

library(tidyverse)

# creates a llist than take the vector and fill rank gaps!!!
llist <- function(x) {
  x <- paste(x, sep = ';', collapse = ';')
  x <- list(x)
  x <- unlist(x)
}

# A) Load data  ====


path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"

db2worms_f <- "SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed_to_worms.rds"

db2worms_f <- list.files(path = path, pattern = db2worms_f, full.names = T)

db2worms <- read_rds(db2worms_f) %>% as_tibble()

query_db_f <- "SSURef_NR99-138.1-query2wormsdb.rds"

query_db <- list.files(path = path, pattern = query_db_f, full.names = T)

query_db <- read_rds(query_db)


# 1) FILTER WORMS   ====

updated_query_db <- query_db %>% left_join(db2worms) %>% filter(status == "accepted")

# updated_query_db %>% view()


# NUMBER OF CONVERTED TAXONS:

into <- unique(updated_query_db$WormsRank)

updated_query_db %>% filter(WormsRank != "not found") %>% count(WormsRank) %>% arrange(match(WormsRank, into)) %>% view()


into <- unique(updated_query_db$Rank)

updated_query_db %>% filter(WormsRank != "not found") %>% count(Rank) %>% arrange(match(Rank, into)) # %>% view()
updated_query_db %>% filter(WormsRank != "not found") %>%  unnest(`Feature ID`) %>% distinct(`Feature ID`)

# 2) CREATE TAX FORMAT =====

recode_ranks <- c("phylum (division)" = "phylum","subphylum"= "phylum","infraphylum"= "phylum",
  "subclass" = "class", "superorder"="order")

recode_ranks <- c("phylum (division)" = "pd","subphylum"= "ps","infraphylum"= "pi",
  "subclass" = "class", "superorder"="order")

ordered_ranks <- unique(updated_query_db$Rank) # c("phylum", "class", "order", "family", "genus", "species")

unique(updated_query_db$WormsRank)

updated_query_db %>%
  mutate(WormsRank = recode_factor(WormsRank, !!!recode_ranks)) %>%
  mutate(WormsRank = factor(WormsRank, levels = ordered_ranks)) %>%
  unnest(`Feature ID`) %>%
  mutate_all(list(~ str_replace_all(., c("none" = NA_character_)))) %>%
  drop_na(WormsRank) %>% view()
  group_by(`Feature ID`) %>%
  distinct(WormsRank, valid_name) %>%
  mutate(WormsRank = substr(WormsRank, 1,1), valid_name = paste0(WormsRank, "__", valid_name)) %>%  
  summarise(across(valid_name, .fns = llist)) %>% view()
  dplyr::rename("Taxon" = "valid_name") -> updated_query_db

updated_query_db %>% view()

# SAVE FASTA AND TAX_FILE =====

file_out <- "/SSURef_NR99-138.1_V9_to_worms_taxonomy.tsv"

file_out <- paste0(path, file_out)

write_tsv(updated_query_db, file = file_out)


# B) Load data =====

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"

db2worms_f <- "pr2_version_4.14.0_SSU-1389f-1510r-db_processed_to_worms.rds"

db2worms_f <- list.files(path = path, pattern = db2worms_f, full.names = T)

db2worms <- read_rds(db2worms_f) %>% as_tibble()

query_db_f <- "pr2_version_4.14.0_SSU-query2wormsdb"

query_db <- list.files(path = path, pattern = query_db_f, full.names = T)

query_db <- read_rds(query_db)

# 1) FILTER WORMS   ====

updated_query_db <- query_db %>% left_join(db2worms) %>% filter(status == "accepted") # 6,474 valid worms ids

nrow(updated_query_db)

# NUMBER OF CONVERTED TAXONS:

into <- unique(updated_query_db$Rank)

updated_query_db %>% filter(WormsRank != "not found") %>% count(Rank) %>% arrange(match(Rank, into)) # %>% view()
updated_query_db %>% filter(WormsRank != "not found") %>%  unnest(`Feature ID`) %>% distinct(`Feature ID`)

# 2) CREATE TAX FORMAT =====

unique(updated_query_db$WormsRank)

updated_query_db <- updated_query_db %>%
  unnest(`Feature ID`) %>%
  mutate_all(list(~ str_replace_all(., c("none" = NA_character_)))) %>%
  drop_na(WormsRank) %>%
  group_by(`Feature ID`) %>%
  distinct(WormsRank, valid_name) %>%
  mutate(valid_name = paste0(WormsRank, "__", valid_name)) %>% # WormsRank = substr(WormsRank, 1,2), 
  summarise(across(valid_name, .fns = llist)) %>%
  dplyr::rename("Taxon" = "valid_name")

# C) CONTACT FILES ====

updated_query_db_out <- read_tsv(paste0(path, "/SSURef_NR99-138.1_V9_to_worms_taxonomy.tsv")) %>%
  rbind(.,updated_query_db)

nrow(updated_query_db_out) # 41,283

# CHECK IF distinct FeatureIDS 

updated_query_db_out %>% distinct(`Feature ID`) # 41,283


# 3 ) FILTER SEQUENCES ----
# cat SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed.fasta pr2_version_4.14.0_SSU-1389f-1510r-db_processed.fasta > SILVA_PR2_V9_MULTI.fasta

dna_f <- paste0(path, "/SILVA-PR2-derep-uniq-1389f-1510r.fasta")

dna <- Biostrings::readDNAStringSet(dna_f) 

str(IDS <- updated_query_db_out %>% pull(`Feature ID`) %>% sort())

sum(keep <- names(dna) %in% IDS) # 41,283

dna <- dna[keep]

# Sanity check

sum(sort(names(dna)) == IDS) # 41,283

file_out <- "/CURATED-1389f-1510r_worms_taxonomy.tsv"

file_out <- paste0(path, file_out)

write_tsv(updated_query_db_out, file = file_out)


seqs_out <- "/CURATED-1389f-1510r_worms_sequences.fasta"

Biostrings::writeXStringSet(dna, filepath =  paste0(path, seqs_out), format = "fasta")


# EXIT =====


# updated_query_db %>% unnest(`Feature ID`) %>% filter(`Feature ID` %in% "AB201232.1.1802") %>% view()
# distinct(names_from, values_from ,`Feature ID`)


# 

# summarise(n = length(values_from), 
#   across(values_from, .fns = list)) %>% arrange(desc(n)) %>% head() %>% pull(`Feature ID`)

which_cols <- c("values_from","names_from","Feature ID")

updated_query_db %>% distinct(WormsRank) %>% pull()

updated_query_db %>% unnest(`Feature ID`) %>% filter(Rank == 'order') %>% pull(`Feature ID`) %>% head()

dbt <- updated_query_db %>% unnest(`Feature ID`) %>% filter(`Feature ID` %in% "ACJG01011253.3973.5773")


dbt  %>% view()
# TEST

dbt %>% 
  mutate_all(list(~ str_replace_all(., c("none" = NA_character_)))) %>%
  pivot_longer(any_of(ranks_ws), values_to = "values_from", names_to = "names_from") %>%
  drop_na(values_from) %>%
  distinct(`Feature ID`, names_from, values_from) %>% 
  pivot_wider(names_from = "names_from", values_from = "values_from") %>% view()



updated_query_db %>%  
  unnest(`Feature ID`) %>%
  mutate_all(list(~ str_replace_all(., c("none" = NA_character_)))) %>%
  pivot_longer(any_of(ranks_ws), values_to = "values_from", names_to = "names_from") %>% 
  drop_na(values_from) %>%
  distinct(`Feature ID`, names_from, values_from) %>% 
  # filter(!grepl("incertae sedis",values_from)) %>%
  dplyr::group_by(`Feature ID`, names_from) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)
  pivot_wider(names_from = "names_from", values_from = "values_from") 

updated_query_db %>%
  mutate(values_from = ifelse(valid_name == "none", Taxon, valid_name)) %>% # get worms valid name
  mutate(names_from = ifelse(WormsRank == "not found", Rank, WormsRank)) %>% # try to merge worms and SILVA ranks
  select(values_from, names_from, `Feature ID`, Wormsid) %>%# values_from, names_from
  unnest(`Feature ID`) %>%
  mutate(names_from = ifelse(names_from %in% merge_phy, "phylum", names_from)) %>%
  pivot_wider(names_from = "names_from", values_from = "values_from") %>% view()

# updated_query_db <- 


merge_phy <- c("phylum (division)","subphylum","infraphylum")

updated_query_db %>%
  filter(WormsRank != "not found") %>%
  mutate(values_from = valid_name) %>%
  mutate(names_from = WormsRank) %>% 
  mutate(names_from = ifelse(names_from %in% merge_phy, "phylum", names_from)) %>% view()
  # mutate(names_from = str_replace_all(names_from, "\\([^()]{0,}\\)", "")) %>% #distinct(names_from) %>% pull()
  filter(!grepl("incertae sedis",values_from)) %>%
  # drop_na() %>%
  # select(all_of(which_cols)) %>%
  unnest(`Feature ID`) %>% 
  group_by(`Feature ID`) %>%
  distinct(names_from, values_from) %>%
  # filter(names_from %in% c(ranks_ws, "species")) %>%
  pivot_wider(names_from = "names_from", values_from = "values_from") %>% view()
  select(any_of(c("Feature ID", ranks_ws))) 


# REVISAR PQ HAY DOS CLASIFICACIONES PARA NEMARTEA INSERT SEDIS
# x <- c("AY210452.1.1881","AY238988.1.1820","AY238989.1.1793","KP270787.1.1793","KP270799.1.1968","MK076311.1.1395")

# db1 %>% filter(`Feature ID` %in% x)

# names2wormsdf("Palaeonemertea")


# 2

# which_rank <- "rank9"
# 
# rank_pos <- which(into %in% which_rank)
# 
# rank9 %>% select(any_of(into[1:rank_pos])) %>% 
#   rename("Taxon"= which_rank) %>%
#   distinct(Taxon, .keep_all = T) 
# 
# tax_df %>% distinct(tax) %>% pull()




# (EXPERIMENTAL SO OMIT) HOW TO REVASE A COLUMN NOT PRESENTED IN A ROW

ranks_ws <- c("kingdom","phylum","class","order","family","genus", "species")

# updated_query_db <- 

# query_db %>% left_join(df) %>%

# Generate worms_df ranks

worms_ranks_df <- df %>%
  mutate(values_from = ifelse(valid_name == "none", Taxon, valid_name)) %>%
  select(any_of(c("values_from", ranks_ws))) %>%
  left_join(updated_query_db) %>%
  unnest(`Feature ID`) %>% distinct()
# pivot_longer(any_of(ranks_ws), values_to = "values_from", names_to = "names_from") %>%
# drop_na(values_from) %>% filter(values_from != 'none') %>% distinct() %>%
# right_join(updated_query_db, by = "Taxon")

# updated_query_db %>%
#   select(any_of(c("Rank","Taxon", "Feature ID","values_from", "names_from", ranks_ws))) %>%
#   pivot_longer(cols = any_of(ranks_ws), values_to = "values_from", names_to = "names_from") %>%
#   mutate(values_from = ifelse(names_from == Rank & Taxon != values_from, NA, values_from)) %>%
#   drop_na(values_from) %>% filter(values_from != 'none')


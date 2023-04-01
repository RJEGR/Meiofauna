

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)
library(tidyverse)
library(rstatix)

# 0) Load data and functions ====

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"


db1_f <- list.files(path = path, pattern = "SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed.rds", full.names = T)


db1 <- read_rds(db1_f)

names(db1)[9] <- "species"

# 1) Input Taxa-groups list ----


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



# A) SSURef_NR99-138.1 ----
# 2) Return DB-ranks groups ----

taxa_list <- c(NON_EUK, PHYLA)
taxa_list <- str_to_lower(taxa_list)
taxa_list <- paste(taxa_list, collapse = "|")

# USING ACROSS AND UNNEST FORMAT
# 2.1) PIVOT LONGER FORMAT ----

omit_cols <- c("Feature ID","dna_width")

into <- db1 %>% select_at(vars(-contains(omit_cols))) %>% names()

# OPTION 1 ====
# UNFILTERED RANKS RESULTED IN A SINGLE SET RESULT IN A SET OF 12,592 SEARCHES
# FILTERED RANKS RESULTED IN A SET OF 3,197


query_db <- db1 %>%
  mutate_at(vars(all_of(into)), list(~ str_to_lower(.))) %>%
  filter_at(vars(all_of(into)), any_vars(grepl(taxa_list, .))) %>% 
  mutate_at(vars(all_of(into)), list(~ str_to_sentence(.))) %>%
  pivot_longer(all_of(into), values_to = "Taxon", names_to = "Rank") %>%
  drop_na(Taxon) %>%
  group_by(Rank, Taxon) %>%
  summarise(across(`Feature ID`, .fns = list)) %>%
  ungroup() %>% arrange(match(Rank, into)) 

query_db %>% count(Rank) %>% arrange(match(Rank, into)) 

# OPTION 2) AS A SINGLE RANK  ====

# which_rank <- "phylum"


# query_db <- db1 %>% 
#   mutate_at(vars(all_of(into)), funs(str_to_lower(.))) %>%
#   filter_at(vars(all_of(into)), any_vars(grepl(taxa_list, .))) %>% 
#   mutate_at(vars(all_of(into)), funs(str_to_sentence(.))) %>%
#   rename("Taxon"= which_rank) %>%
#   dplyr::select(`Feature ID`, Taxon) %>%
#   group_by(Taxon) %>%
#   summarise(n = length(`Feature ID`), across(`Feature ID`, .fns = list)) %>%
#   distinct()

# 3) WORMS CONVERTION =====

# query_db <- query_db %>% sample_n(30)

str(query <- query_db %>% pull(Taxon)) # 3197 FOR SILVA

# df <- names2wormsdf(query, accepted = TRUE, marine = TRUE)

file_out <- "/SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed_to_worms.rds"

# write_rds(df, file = paste0(path, file_out))

df <- read_rds(paste0(path, file_out)) %>% as_tibble()

# df %>% view()

# 4) IF TAXON IS NOT IN WORMS, LETS FILL valid_name and WormsRank  ====

# c(,phylum","phylum (division)","subphylum","infraphylum", "class","subclass" ,"order", "superorder", "suborder", "infraorder", "family","genus","species","subspecies","subspecies","variety")

updated_query_db <- query_db %>% left_join(df) %>%
  mutate(values_from = ifelse(valid_name == "none", Taxon, valid_name)) %>% # get worms valid name
  mutate(names_from = ifelse(WormsRank == "not found", Rank, WormsRank)) %>% # try to merge worms and SILVA ranks
  select(values_from, names_from, `Feature ID`) # values_from, names_from

# updated_query_db %>% unnest(`Feature ID`) %>% filter(`Feature ID` %in% "AB201232.1.1802") %>% view()
  # distinct(names_from, values_from ,`Feature ID`)


# (EXPERIMENTAL SO OMIT) HOW TO REVASE A COLUMN NOT PRESENTED IN A ROW

ranks_ws <- c("kingdom","phylum","class","order","family","genus")

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


# 5) UNNEST IDS AND PIVOT WIDER =====

# summarise(n = length(values_from), 
#   across(values_from, .fns = list)) %>% arrange(desc(n)) %>% head() %>% pull(`Feature ID`)

which_cols <- c("values_from","names_from","Feature ID")

# CONS

updated_query_db <- updated_query_db %>%
  filter(!grepl("incertae sedis",values_from)) %>%
  drop_na() %>%
  # select(all_of(which_cols)) %>%
  unnest(`Feature ID`) %>% 
  group_by(`Feature ID`) %>%
  distinct(names_from, values_from) %>%
  # filter(names_from %in% c(ranks_ws, "species")) %>%
  pivot_wider(names_from = "names_from", values_from = "values_from") %>% 
  select(any_of(c("Feature ID", ranks_ws))) 


# REVISAR PQ HAY DOS CLASIFICACIONES PARA NEMARTEA INSERT SEDIS
# x <- c("AY210452.1.1881","AY238988.1.1820","AY238989.1.1793","KP270787.1.1793","KP270799.1.1968","MK076311.1.1395")

# db1 %>% filter(`Feature ID` %in% x)

# names2wormsdf("Palaeonemertea")



# pr2_version_4.14.0_SSU-1389f-1510r-db_processed_to_worms.rds
# SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed.rds

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)
library(tidyverse)

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

file_out <- "/SSURef_NR99-138.1-query2wormsdb.rds"

write_rds(query_db, file = paste0(path, file_out))

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

str(query <- query_db %>% pull(Taxon)) # 3197 FOR SILVA

df <- names2wormsdf(query, accepted = TRUE, marine = TRUE)

# SAVE =====

file_out <- "/SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed_to_worms.rds"

write_rds(df, file = paste0(path, file_out))

quit(save = "no")



# df <- read_rds(paste0(path, file_out)) %>% as_tibble()
# 
# df %>% filter(grepl(taxa_list, str_to_lower(Taxon))) %>% view()
# 
# query_db %>% left_join(df) %>% view()



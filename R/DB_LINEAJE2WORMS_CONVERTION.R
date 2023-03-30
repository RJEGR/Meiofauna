
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)
library(rstatix)


# 0) Load data ====

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"


db2_f <- list.files(path = path, pattern = "pr2_version_4.14.0_SSU-1389f-1510r.rds", full.names = T)
db1_f <- list.files(path = path, pattern = "SSURef_NR99-138.1-dna-seqs-euk-derep-uniq-1389f-1510r.rds", full.names = T)


db1 <- read_rds(db1_f) %>% select_at(vars(contains("rank")))
db2 <- read_rds(db2_f)  %>% select_at(vars(contains("rank")))

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

# 2) Return DB-ranks groups ----

taxa_list <- c(NON_EUK, PHYLA)
taxa_list <- str_to_lower(taxa_list)
taxa_list <- paste(taxa_list, collapse = "|")


into <- db1 %>% select_at(vars(contains("rank"))) %>% names()

db1 <- db1 %>% 
  mutate_at(vars(all_of(into)), funs(str_to_lower(.))) %>%
  filter_at(vars(all_of(into)), any_vars(grepl(taxa_list, .))) %>% 
  mutate_at(vars(all_of(into)), funs(str_to_sentence(.))) %>%
  mutate_at(vars(all_of(into)), funs(gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(into),  ~na_if(., "")) %>%
  mutate_at(vars(into), funs(str_replace_all(., c("_" = " "))))


which_rank <- "rank3"

# rank_pos <- which(into %in% which_rank)

query_db <- db1 %>% rename("Taxon"= which_rank) %>% distinct(Taxon) %>% drop_na()

str(query <- query_db %>% pull())

# Rename tax2worms

tax_rank(head(query), db = "worms")

# (Worms) An authoritative classification and catalogue of marine names


# tax_rank(query[1], db = "worms")

names2wormsdf <- function(mynames, accepted = FALSE, marine = TRUE) {
  
  require(taxize)
  
  worms <- taxize::get_wormsid(mynames, accepted = accepted, ask = FALSE, marine_only = marine)
  
  id <- c(worms)
  
  uri <- attributes(worms)$uri
  
  rankLev <- tax_rank(id, db = "worms")
  
  out <- data.frame(Taxon = mynames, rankLevel = rankLev, wormsid = id, uri = uri )
  
  return(out)
}

# taxize::get_wormsid(query[1], accepted = T, ask = FALSE)

worms <- names2wormsid(head(query), accepted = TRUE)

# source('~/Documents/GitHub/zooplankton_gom/names2worms.R')
# names_ <- names2worms_(query) # take time ...

names_ %>% 
  # drop_na(Phylum_wm) %>% 
  # distinct(ori_name) %>%
  pull() -> marineMetazoa


db1 %>% left_join(names_, by = c("Taxon", "ori_name"))

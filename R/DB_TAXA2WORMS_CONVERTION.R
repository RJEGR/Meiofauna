# Ricardo Gomez-Reyes, 2023

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

query_db <- db1 %>% rename("Taxon"= which_rank) %>% distinct(Taxon) %>% drop_na()

str(query <- query_db %>% pull())


names2wormsdf_ <- function(queries, accepted = FALSE, marine = TRUE) {
  
  # (Worms) An authoritative classification and catalogue of marine names
  
  require(taxize)
  
  worms <- taxize::get_wormsid(queries, accepted = accepted, ask = FALSE, marine_only = marine, messages = FALSE)
  
  not_record <- "not found"
  
  if(attributes(worms)$match == not_record) {
    
    out <- data.frame(Taxon = queries, Rank = not_record, Wormsid = not_record, URI = not_record)
    
  } else {
    
    ID <- c(worms)
    
    URI <- attributes(worms)$uri
    
    rankLev <- taxize::tax_rank(ID, db = "worms")
    
    rankLev <- unlist(rankLev)
    
    names(rankLev) <- NULL
    
    out <- data.frame(Taxon = queries, Rank = rankLev, Wormsid = ID, URI = URI)
    
  }
  
  
  return(out)
}

# Test:
# names2wormsdf_(query[20], accepted = TRUE, marine = TRUE)

names2wormsdf <- function(queries, ...) {
  
  # wrapped version from names2wormsdf_, include ProgressBar
  maxq <- length(queries)
  
  pb <- txtProgressBar(min = 0, max = maxq, style = 3, width = 50, char = "=")
  
  out <- list()
  
  for(i in 1:maxq) {
    
    
    out[[i]] <- names2wormsdf_(queries[i], accepted = FALSE, marine = TRUE)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb) 
  
  out <- do.call(rbind, out)
  
  return(out)
  
}

get_wm_record <- function(wormsid) {
  
  which_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  which_ranks <- str_to_lower(which_ranks)
  
  wormsid <- as.integer(wormsid)
  
  # 1) Choose the fields in WoRMS-history-Records 
  
  which_fields <- c("AphiaID", "scientificname", "status", "rank", 
    "valid_AphiaID", "valid_name", "valid_authority", 
    "isMarine", "isBrackish", "isFreshwater", "isTerrestrial", "modified")

  
  # 2) Retrive records 
  
  records_tibble <- worrms::wm_record(wormsid) 
  
  keep_ranks <- names(records_tibble) %in% which_ranks
  
  # If lower rank lets "stack" the rank levels
  
  which_ranks <- names(records_tibble)[keep_ranks]
  
  all_cols <- c(which_fields, which_ranks)
  
  records_tibble <- records_tibble %>% select_at(all_of(all_cols))
  
  return(records_tibble)
  
}

# Example how to use:

get_wm_record("22565")


  
df <- names2wormsdf(query,accepted = TRUE, marine = TRUE)






db1 %>% left_join(names_, by = c("Taxon", "ori_name"))

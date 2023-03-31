# Ricardo Gomez-Reyes, 2023

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)
library(rstatix)


# 0) Load data ====

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"


db1_f <- list.files(path = path, pattern = "SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed.rds", full.names = T)


db2_f <- list.files(path = path, pattern = "pr2_version_4.14.0_SSU-1389f-1510r-db_processed.rds", full.names = T)


db1 <- read_rds(db1_f)



# db2 <- read_rds(db2_f)  %>% select_at(vars(-contains(omit_cols)))

# STACKS RANKS AS A SILVA =====


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

# USING ACROSS AND UNNEST FORMAT
# 2.1) PIVOT LONGER FORMAT ----

omit_cols <- c("Feature ID","dna_width")

into <- db1 %>% select_at(vars(-contains(omit_cols))) %>% names()

# OPTION 1, 
# UNFILTERED RANKS RESULTED IN A SINGLE SET RESULT IN A SET OF 12,592 SEARCHES
# FILTERED RANKS RESULTED IN A SET OF 3,197

query_db <- db1 %>%
  mutate_at(vars(all_of(into)), funs(str_to_lower(.))) %>%
  filter_at(vars(all_of(into)), any_vars(grepl(taxa_list, .))) %>% 
  mutate_at(vars(all_of(into)), funs(str_to_sentence(.))) %>%
  pivot_longer(all_of(into), values_to = "Taxon", names_to = "Rank") %>%
  drop_na(Taxon) %>%
  group_by(Rank, Taxon) %>%
  summarise(across(`Feature ID`, .fns = list)) %>%
  ungroup() %>% arrange(match(Rank, into)) %>%
  distinct()


# OPTION 2) AS A SINGLE RANK

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

# THEN

query_db %>% 
  unnest(`Feature ID`) %>%
  pivot_wider(names_from = "Rank", values_from = "Name")


db1 %>% filter(`Feature ID` %in% "L26517.1.1880")


# 3) WORMS CONVERTION =====


query_db_test <- query_db %>% sample_n(30)

str(query <- query_db_test %>% pull(Taxon))


names2wormsdf_ <- function(queries, accepted = FALSE, marine = TRUE) {
  
  # (Worms) An authoritative classification and catalogue of marine names
  
  which_fields <- c("AphiaID", "scientificname", "status", "rank", 
    "valid_AphiaID", "valid_name", "valid_authority", 
    "isMarine", "isBrackish", "isFreshwater", "isTerrestrial", "modified")
  
  
  
  require(taxize)
  
  worms <- taxize::get_wormsid(queries, accepted = accepted, ask = FALSE, marine_only = marine, messages = FALSE)
  
  not_record <- "not found"
  
  if(attributes(worms)$match == not_record) {
    
    out <- data.frame(Taxon = queries, WormsRank = not_record, Wormsid = not_record, URI = not_record)
    
    # add emtpy wm_records
    
    df_fields <- data.frame(matrix(ncol = length(which_fields)))
    
    names(df_fields) <- which_fields
    
    df_fields[1,] <- "none"
    
    out <- cbind(out, df_fields)
    
    
  } else {
    
    ID <- c(worms)
    
    URI <- attributes(worms)$uri
    
    rankLev <- taxize::tax_rank(ID, db = "worms")
    
    rankLev <- unlist(rankLev)
    
    names(rankLev) <- NULL
    
    out <- data.frame(Taxon = queries, WormsRank = rankLev, Wormsid = ID, URI = URI)
    
    # add emtpy wm_records
    
    df_fields <- get_wm_record(ID)
    
    out <- cbind(out, df_fields)

    
  }
  
  
  return(out)
}

# Example how to use::
# names2wormsdf_(query[10], accepted = TRUE, marine = TRUE)

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

  # df_fields <- data.frame(matrix(ncol = length(which_fields)))
  # names(df_fields) <- which_fields
  # df_fields[1,] <- "none"
  
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
# get_wm_record("22565")

  
df <- names2wormsdf(query, accepted = TRUE, marine = TRUE)

df %>% as_tibble()

query_db_test %>% left_join(df)

get_wm_record("130130") %>% view()

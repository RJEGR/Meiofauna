# Ricardo Gomez-Reyes, 2023

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
  library(rstatix)
  library(worrms)
  library(taxize)
})


# FUNCTIONS ====

get_wm_record <- function(wormsid) {
  
  which_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  which_ranks <- str_to_lower(which_ranks)
  
  wormsid <- as.integer(wormsid)
  
  # 1) Choose the fields in WoRMS-history-Records 
  
  which_fields <- c("AphiaID", "scientificname", "status", 
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
# get_wm_record("22565")

names2wormsdf_ <- function(queries, accepted = FALSE, marine = TRUE) {
  
  # (Worms) An authoritative classification and catalogue of marine names
  
  which_fields <- c("AphiaID", "scientificname", "status", 
    "valid_AphiaID", "valid_name", "valid_authority", 
    "isMarine", "isBrackish", "isFreshwater", "isTerrestrial", "modified")
  
  which_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  
  which_ranks <- str_to_lower(which_ranks)
  
  all_cols <- c(which_fields, which_ranks)
  
  require(taxize)
  
  worms <- taxize::get_wormsid(queries, accepted = accepted, ask = FALSE, marine_only = marine, messages = FALSE)
  
  not_record <- "not found"
  
  if(attributes(worms)$match == not_record) {
    
    out <- data.frame(Taxon = queries, WormsRank = not_record, Wormsid = not_record, URI = not_record)
    
    # add emtpy wm_records
    
    df_fields <- data.frame(matrix(ncol = length(all_cols)))
    
    names(df_fields) <- all_cols
    
    df_fields[1,] <- "none"
    
    out <- cbind(out, df_fields)
    
    
  } else {
    
    ID <- c(worms)
    
    URI <- attributes(worms)$uri
    
    rankLev <- taxize::tax_rank(ID, db = "worms")
    
    rankLev <- unlist(rankLev)
    
    names(rankLev) <- NULL
    
    out <- data.frame(Taxon = queries, WormsRank = rankLev, Wormsid = ID, URI = URI)
    
    # add wm_records
    
    df_fields <- get_wm_record(ID)
    
    out <- cbind(out, df_fields)
    
  }
  
  system(paste0("mkdir -p ", getwd(), "/NAMES2WORMS_TEMP_DIR"))
  
  FileName <- gsub("[[:space:]]", "", as.character(out[,1:3]))
  
  FileName <- paste(FileName, collapse = "_")
  
  FileName <- digest::digest(FileName, serialize = F) # Create hash name
  
  write_lines(out, file = paste0(getwd(), "/NAMES2WORMS_TEMP_DIR", "/", FileName, ".tmp"))
  
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


# # B) PR2 V 4.14 SSU  ----

# 0) Load data ====

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"


db2_f <- list.files(path = path, pattern = "pr2_version_4.14.0_SSU-1389f-1510r-db_processed.rds", full.names = T)

db2 <- read_rds(db2_f)  # %>% select_at(vars(-contains(omit_cols)))

names(db2)[10] <- "species"

# 2) Return DB-ranks groups ----

taxa_list <- c(NON_EUK, PHYLA)
taxa_list <- str_to_lower(taxa_list)
taxa_list <- paste(taxa_list, collapse = "|")

# 2.1) PIVOT LONGER FORMAT ----

omit_cols <- c("Feature ID","dna_width")

into <- db2 %>% select_at(vars(-contains(omit_cols))) %>% names()

# mutate_at(vars(all_of(into)), funs(str_to_lower(.))) %>%
# mutate_at(vars(all_of(into)), list(~ str_to_lower(.)))


query_db <- db2 %>%
  mutate_at(vars(all_of(into)), list(~ str_to_lower(.))) %>%
  filter_at(vars(all_of(into)), any_vars(grepl(taxa_list, .))) %>% 
  mutate_at(vars(all_of(into)), list(~ str_to_sentence(.))) %>%
  mutate_at(into, list(~ str_replace_all(., c("/" = "_")))) %>%
  pivot_longer(all_of(into), values_to = "Taxon", names_to = "Rank") %>%
  drop_na(Taxon) %>%
  group_by(Rank, Taxon) %>%
  summarise(across(`Feature ID`, .fns = list)) %>%
  ungroup() %>% arrange(match(Rank, into)) 

query_db %>% count(Rank) %>% arrange(match(Rank, into)) 

file_out <- "/pr2_version_4.14.0_SSU-query2wormsdb.rds"

write_rds(query_db, file = paste0(path, file_out))



# IF STOPPED ====

tmp_dir <- paste0(getwd(), "/NAMES2WORMS_TEMP_DIR")
#
tmpFiles <- list.files(path = tmp_dir, pattern = ".tmp$", full.names = T)
#

# str(tmpFiles)

cat("\n")


# if empty file omit (CREATE)

read_lines_ <- function(x) {
  
  which_f <- c("Taxon", "WormsRank", "Wormsid", "URI")
  
  which_fields <- c("AphiaID", "scientificname", "status", 
    "valid_AphiaID", "valid_name", "valid_authority", 
    "isMarine", "isBrackish", "isFreshwater", "isTerrestrial", "modified")
  
  which_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  
  which_ranks <- str_to_lower(which_ranks)
  
  all_cols <- c(which_f, which_fields, which_ranks)
  
  out <- data.frame(matrix(ncol = 21, nrow = length(x)))
  
  maxq <- length(x)
  
  pb <- txtProgressBar(min = 0, max = maxq, style = 3, width = 50, char = "=")
  
  for(i in 1:length(x)) {
    
    
    out[i,] <- read_lines(x[i])
    
    setTxtProgressBar(pb, i)
    
    }
  
  close(pb) 
  
  names(out) <- all_cols
  
  out <- out %>% as_tibble()
  
  return(out)
}


read_lines_(tmpFiles[1:2])

df <- read_lines_(tmpFiles)

# write_rds(df, file = paste0(tmp_dir, "/tmpFiles.rds"))

# df <- read_rds(paste0(tmp_dir, "/tmpFiles.rds"))

#
str(recorded <- df %>% pull(Taxon))
#
query_db <- query_db %>% filter(!Taxon %in% recorded)
#

# 3) WORMS CONVERTION =====

cat("\n")

str(query <- query_db %>% pull(Taxon)) # 18,176 FOR PR2

cat("\n")

df2 <- names2wormsdf(query, accepted = TRUE, marine = TRUE)

out <- rbind(df, df2)

cat("\n")

file_out <- "/pr2_version_4.14.0_SSU-1389f-1510r-db_processed_to_worms.rds"

write_rds(out, file = paste0(path, file_out))

quit(save = "no")






# DB_MERGING.R

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)

library(tidyverse)

# A) Load data  ====

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"

db2worms_f <- "SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed_to_worms.rds"

db2worms_f <- list.files(path = path, pattern = db2worms_f, full.names = T)

db2worms <- read_rds(db2worms_f) %>% as_tibble() 

db2worms <- db2worms %>% mutate_all(list(~ str_replace_all(., c("none" = NA_character_))))

db2worms <- db2worms %>% mutate_all(list(~ str_replace_all(., c("NA" = NA_character_))))


db2worms <- db2worms %>% filter(WormsRank != "not found")

db2worms <- db2worms %>% filter(isMarine == 1)

query_db_f <- "SSURef_NR99-138.1-query2wormsdb.rds"

query_db <- list.files(path = path, pattern = query_db_f, full.names = T)

query_db <- read_rds(query_db)

# 1) FILTER WORMS   ====

updated_query_db <- query_db %>% left_join(db2worms) %>% filter(status == "accepted")

# updated_query_db %>% unnest(`Feature ID`) %>% view()

# IS MARINE? 

db2worms %>% filter(status == "accepted") %>% count(isMarine)

# previous rank is == to wormsRank ?

ranks <- c("kingdom","phylum","class","order","family","genus")

# deal w/ species resolution


findMaxRank <- function(df, colRanks = c("kingdom","phylum","class","order","family","genus")) {
  
  .df <- df %>% select_at(colRanks)
  
  which_rank <-  function(x) { 
    r <- colRanks
    n <- length(r)
    max(seq(1:n)[complete.cases(x)]) 
  }  
  
  df$mr <- apply(.df, 1, which_rank)
  
  df <- df %>%
    mutate(s = ifelse(mr == length(ranks) & valid_name != genus, valid_name, NA)) %>%
    mutate(mr = ifelse(!is.na(s), mr+1, mr))
  
  return(df)
}


.updated_query_db <- updated_query_db %>%
  unnest(`Feature ID`) %>%
  select_at(c("Feature ID", ranks, "valid_name"))


.updated_query_db <- findMaxRank(.updated_query_db)

table(.updated_query_db$mr)

table(findMaxRank(db2worms %>%  select_at(c(ranks, "valid_name")))$mr)



.updated_query_db <- .updated_query_db %>% 
  # mutate(s = ifelse(mr == length(ranks) & valid_name != genus, valid_name, NA)) %>%
  # mutate(mr = ifelse(!is.na(s), mr+1, mr)) %>%
  mutate(s = str_replace_all(s, c(" " = "_"))) %>%
  group_by(`Feature ID`)  %>%
  arrange(desc(mr), .by_group = T) 

.updated_query_db <- .updated_query_db %>% group_by(`Feature ID`) %>% slice_head(n = 1)


table(.updated_query_db$mr)
# 

.updated_query_db <- .updated_query_db  %>% select_at(c("Feature ID", ranks, "s"))

names(.updated_query_db)[which(names(.updated_query_db) %in% ranks)] <- substr(ranks, 1,1) 

.updated_query_db <- .updated_query_db %>%
  mutate(k = ifelse(is.na(k), "k__", paste0("k__", k))) %>%
  mutate(p = ifelse(is.na(p), "p__", paste0("p__", p))) %>%
  mutate(c = ifelse(is.na(c), "c__", paste0("c__", c))) %>%
  mutate(o = ifelse(is.na(o), "o__", paste0("o__", o))) %>%
  mutate(f = ifelse(is.na(f), "f__", paste0("f__", f))) %>%
  mutate(g = ifelse(is.na(g), "g__", paste0("g__", g))) %>%
  mutate(s = ifelse(is.na(s), "s__", paste0("s__", s))) %>%
  unite("Taxon", k:s, sep = ";")
# 

# SAVE FASTA AND TAX_FILE =====

file_out <- "/SSURef_NR99-138.1_V9_to_worms_taxonomy.tsv"

file_out <- paste0(path, file_out)

write_tsv(.updated_query_db, file = file_out)


# B) Load data =====

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"

db2worms_f <- "pr2_version_4.14.0_SSU-1389f-1510r-db_processed_to_worms.rds"

db2worms_f <- list.files(path = path, pattern = db2worms_f, full.names = T)

db2worms <- read_rds(db2worms_f) %>% as_tibble()

db2worms <- db2worms %>% mutate_all(list(~ str_replace_all(., c("none" = NA_character_))))

db2worms <- db2worms %>% mutate_all(list(~ str_replace_all(., c("NA" = NA_character_))))

db2worms <- db2worms %>% filter(WormsRank != "not found")

db2worms <- db2worms %>% filter(isMarine == 1)


query_db_f <- "pr2_version_4.14.0_SSU-query2wormsdb"

query_db <- list.files(path = path, pattern = query_db_f, full.names = T)

query_db <- read_rds(query_db)

# 1) FILTER WORMS   ====

updated_query_db <- query_db %>% left_join(db2worms) %>% filter(status == "accepted") # 6,474 valid worms ids

db2worms %>% filter(status == "accepted") %>% count(WormsRank)

nrow(updated_query_db) # 6330

# 2) CREATE TAX FORMAT =====

.updated_query_db <- updated_query_db %>%
  unnest(`Feature ID`) %>%
  select_at(c("Feature ID", ranks, "valid_name")) 
  

.updated_query_db <- findMaxRank(.updated_query_db)

table(.updated_query_db$mr)
table(findMaxRank(db2worms %>%  select_at(c(ranks, "valid_name")))$mr)


.updated_query_db <- .updated_query_db %>% 
  mutate(s = str_replace_all(s, c(" " = "_"))) %>%
  group_by(`Feature ID`)  %>%
  arrange(desc(mr), .by_group = T) 

.updated_query_db <- .updated_query_db %>% group_by(`Feature ID`) %>% slice_head(n = 1)

table(.updated_query_db$mr)

.updated_query_db <- .updated_query_db  %>% select_at(c("Feature ID", ranks, "s"))

names(.updated_query_db)[which(names(.updated_query_db) %in% ranks)] <- substr(ranks, 1,1) 

.updated_query_db <- .updated_query_db %>%
  # mutate_all(list(~ str_replace_all(., c("NA" = NA_character_)))) %>%
  mutate(k = ifelse(is.na(k), "k__", paste0("k__", k))) %>%
  mutate(p = ifelse(is.na(p), "p__", paste0("p__", p))) %>%
  mutate(c = ifelse(is.na(c), "c__", paste0("c__", c))) %>%
  mutate(o = ifelse(is.na(o), "o__", paste0("o__", o))) %>%
  mutate(f = ifelse(is.na(f), "f__", paste0("f__", f))) %>%
  mutate(g = ifelse(is.na(g), "g__", paste0("g__", g))) %>%
  mutate(s = ifelse(is.na(s), "s__", paste0("s__", s))) %>%
  unite("Taxon", k:s, sep = ";")

# C) CONCAT FILES ====

updated_query_db_out <- read_tsv(paste0(path, "/SSURef_NR99-138.1_V9_to_worms_taxonomy.tsv")) %>%
  rbind(., .updated_query_db)

nrow(updated_query_db_out) # 41,283

# CHECK IF distinct FeatureIDS 

nrow(updated_query_db_out %>% distinct(`Feature ID`)) # 41,283


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



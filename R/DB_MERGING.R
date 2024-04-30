# DB_MERGING.R

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)

library(tidyverse)

# creates a llist than take the vector and fill rank gaps!!!


fill_ranks <- function(x) {
  
  # x: vector of ranks
  # y: backbone of ranks
  
  y <- c("k","p", "c", "o", "f", "g", "s")
  
  # m <- matrix(ncol = length(y))
  
  m <- NULL
  
  for(i in length(y):1) { 
    
    searchrank <- grepl(paste0("^", y[i]), x)
    
    if(sum(searchrank) > 0) {
      
      pos <- which(searchrank)
      
      pos <- sample(pos, size = 1)
      
      rank <- x[pos]
      
      m[i] <- rank
      
    } else {
      
      rank <- paste0(y[i], "__")
      
      m[i] <- rank
      
    }
    
  }
  
  m <- paste0(unique(m), collapse = ";")
  
  return(m)
  
}

# A) Load data  ====

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"

db2worms_f <- "SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed_to_worms.rds"

db2worms_f <- list.files(path = path, pattern = db2worms_f, full.names = T)

db2worms <- read_rds(db2worms_f) %>% as_tibble() 

db2worms <- db2worms %>% filter(WormsRank != "not found")

db2worms <- db2worms %>% filter(isMarine == 1)

query_db_f <- "SSURef_NR99-138.1-query2wormsdb.rds"

query_db <- list.files(path = path, pattern = query_db_f, full.names = T)

query_db <- read_rds(query_db)

query_db %>% left_join(db2worms) %>% view()


# 1) FILTER WORMS   ====

updated_query_db <- query_db %>% left_join(db2worms) %>% filter(status == "accepted")

updated_query_db %>%
  unnest(`Feature ID`) %>% view()
# IS MARINE? 

db2worms %>% filter(status == "accepted") %>% count(isMarine)

# previous rank is == to wormsRank ?

ranks <- c("kingdom","phylum","class","order","family","genus")

# updated_query_db %>% 
#   mutate(q = Rank == WormsRank) %>% 
#   filter(q == 0) %>% # view()
#   count(q)

# updated_query_db %>%
#   unnest(`Feature ID`)  %>%
#   view()

# deal w/ species resolution

updated_query_db %>% filter(Rank == "species") %>% view()

check_df <- updated_query_db %>%
  unnest(`Feature ID`) %>%
  select_at(c("Feature ID", ranks, "scientificname"))



str(cc <- apply(check_df, 1, function(x) { sum(complete.cases(x)) }))

# which(is.na(x))})

check_df <- check_df  %>% mutate(cc = cc-2)

view(check_df)

check_df %>% group_by(`Feature ID`)  %>%
  slice_max(cc) %>% 
  arrange(`Feature ID`) %>% 
  mutate(species = ifelse(cc == 6 & scientificname != genus, scientificname, NA))

check_df <- check_df  %>% select_at(c("Feature ID", ranks, "species", "cc"))

check_df <- check_df %>% ungroup() %>% distinct()

check_df  %>% 
  dplyr::group_by(`Feature ID`) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) %>% left_join(check_df)  %>% view()

updated_query_db %>%
  unnest(`Feature ID`) %>%
  select_at(c("Feature ID", ranks, "scientificname")) %>% 
  arrange(`Feature ID`) %>%
  group_by(`Feature ID`) %>%
  distinct_at(ranks, .keep_all = F)


# updated_query_db %>% view()


# NUMBER OF CONVERTED TAXONS:

into <- unique(updated_query_db$WormsRank)

updated_query_db %>% filter(WormsRank != "not found") %>% 
  count(WormsRank) %>% 
  arrange(match(WormsRank, into)) %>% view()


into <- unique(updated_query_db$Rank)

updated_query_db %>% filter(WormsRank != "not found") %>% count(Rank) %>% arrange(match(Rank, into)) # %>% view()
updated_query_db %>% filter(WormsRank != "not found") %>%  unnest(`Feature ID`) %>% distinct(`Feature ID`)

# 2) CREATE TAX FORMAT =====

updated_query_db %>% count(WormsRank)

ordered_ranks <- c("p", "c", "o", "f", "g", "s") # substr(unique(updated_query_db$Rank),1,1) #

c("kingdom" = "k",
  "phylum" = "p",
  "phylum (division)",
  "subphylum" = "",
  "infraphylum"= "",
  "class" = "c",
  "subclass" = "",
  
  
  )

structure(ordered_ranks, names = c("phylum", "class", "order", "family", "genus", "species"))
  
recode_ranks <- c("phylum (division)" = "phylum","subphylum"= "phylum","infraphylum"= "phylum",
  "subclass" = "class", "superorder"="order")


# 
# Using WormsRank we duplicate records as recode_ranks dont known how to deal w/ sub,infra,super divisions
updated_query_db %>%
  mutate(WormsRank = recode_factor(WormsRank, !!!recode_ranks)) %>%
  unnest(`Feature ID`) %>%
  mutate_all(list(~ str_replace_all(., c("none" = NA_character_)))) %>%
  drop_na(WormsRank) %>%
  group_by(`Feature ID`) %>%
  distinct(valid_name, WormsRank) %>%
  mutate(valid_name = str_replace_all(valid_name, c(" " = "_"))) %>%
  mutate(WormsRank = substr(WormsRank, 1,1), valid_name = paste0(WormsRank, "__", valid_name)) %>%
  mutate(WormsRank = factor(WormsRank, levels = ordered_ranks)) %>%
  filter(WormsRank %in% ordered_ranks) %>%
  summarise(across(valid_name, .fns = fill_ranks)) %>%
  dplyr::rename("Taxon" = "valid_name") -> updated_query_db
# 
# updated_query_db <- updated_query_db %>%
#   unnest(`Feature ID`) %>%
#   mutate_all(list(~ str_replace_all(., c("none" = NA_character_)))) %>%
#   drop_na(Rank) %>% 
#   group_by(`Feature ID`) %>%
#   distinct(valid_name, Rank) %>%
#   mutate(valid_name = str_replace_all(valid_name, c(" " = "_"))) %>%
#   mutate(Rank = substr(Rank, 1,1), valid_name = paste0(Rank, "__", valid_name)) %>%
#   mutate(Rank = factor(Rank, levels = ordered_ranks)) %>%
#   filter(Rank %in% ordered_ranks) %>%
#   pivot_wider(names_from = Rank, values_from = valid_name) %>%
#   mutate(p = ifelse(is.na(p), "p__", p)) %>%
#   mutate(c = ifelse(is.na(c), "c__", c)) %>%
#   mutate(o = ifelse(is.na(o), "o__", o)) %>%
#   mutate(f = ifelse(is.na(f), "f__", f)) %>%
#   mutate(g = ifelse(is.na(g), "g__", g)) %>%
#   mutate(s = ifelse(is.na(s), "s__", s)) %>%
#   unite("Taxon", p:s, sep = ";")


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

db2worms %>% filter(status == "accepted") %>% count(WormsRank)

nrow(updated_query_db)

# NUMBER OF CONVERTED TAXONS:

into <- unique(updated_query_db$WormsRank)

updated_query_db %>% filter(WormsRank != "not found") %>% count(WormsRank) %>% arrange(match(WormsRank, into)) %>% view()

into <- unique(updated_query_db$Rank)

updated_query_db %>% filter(WormsRank != "not found") %>% count(Rank) %>% arrange(match(Rank, into)) # %>% view()
updated_query_db %>% filter(WormsRank != "not found") %>%  unnest(`Feature ID`) %>% distinct(`Feature ID`)

# 2) CREATE TAX FORMAT =====

# unique(updated_query_db$WormsRank)

recode_ranks <- c("infrakingdom"="kingdom", 
  "phylum (division)" = "phylum","subphylum"= "phylum","infraphylum"= "phylum",
  "subclass" = "class", "superorder"="order")

f_ranks <- c("kingdom","phylum", "class", "order", "family", "genus", "species")

ordered_ranks <- c("k","p", "c", "o", "f", "g", "s")

# updated_query_db %>%
#   mutate(WormsRank = recode_factor(WormsRank, !!!recode_ranks)) %>%
#   filter(WormsRank %in% f_ranks) %>%
#   unnest(`Feature ID`) %>%
#   mutate_all(list(~ str_replace_all(., c("none" = NA_character_)))) %>%
#   drop_na(WormsRank) %>% 
#   group_by(`Feature ID`) %>%
#   distinct(WormsRank, valid_name) %>%
#   mutate(valid_name = str_replace_all(valid_name, c(" " = "_"))) %>%
#   mutate(WormsRank = substr(WormsRank, 1,1), valid_name = paste0(WormsRank, "__", valid_name)) %>%
#   # filter(WormsRank %in% ordered_ranks) %>%
#   mutate(WormsRank = factor(WormsRank, levels = ordered_ranks)) %>%
#   summarise(across(valid_name, .fns = fill_ranks)) %>% 
#   dplyr::rename("Taxon" = "valid_name") -> updated_query_db


updated_query_db %>% count(Rank)

.updated_query_db <- updated_query_db %>%
  mutate(WormsRank = recode_factor(WormsRank, !!!recode_ranks)) %>%
  unnest(`Feature ID`) %>%
  mutate_all(list(~ str_replace_all(., c("none" = NA_character_)))) %>%
  drop_na(Rank) %>% 
  group_by(`Feature ID`) %>%
  distinct(valid_name, Rank) %>%
  mutate(valid_name = str_replace_all(valid_name, c(" " = "_"))) %>%
  mutate(Rank = substr(Rank, 1,1), valid_name = paste0(Rank, "__", valid_name)) %>%
  mutate(Rank = factor(Rank, levels = ordered_ranks)) %>%
  filter(Rank %in% ordered_ranks) 

.updated_query_db %>%
  dplyr::group_by(`Feature ID`, Rank) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) %>% count(Rank)

glue <- function(x) paste0(unique(x), collapse = "|")

.updated_query_db %>%
  pivot_wider(names_from = Rank, values_from = valid_name, values_fn = glue) %>%
  filter()
  mutate(p = ifelse(is.na(p), "p__", p)) %>%
  mutate(c = ifelse(is.na(c), "c__", c)) %>%
  mutate(o = ifelse(is.na(o), "o__", o)) %>%
  mutate(f = ifelse(is.na(f), "f__", f)) %>%
  mutate(g = ifelse(is.na(g), "g__", g)) %>%
  mutate(s = ifelse(is.na(s), "s__", s)) %>%
  unite("Taxon", p:s, sep = ";")

# C) CONCAT FILES ====

updated_query_db_out <- read_tsv(paste0(path, "/SSURef_NR99-138.1_V9_to_worms_taxonomy.tsv")) %>%
  rbind(.,updated_query_db)

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


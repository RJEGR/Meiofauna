
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)

library(rstatix)


# 0) Load data ====

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/DB_COMPARISON"


db1_f <- list.files(path = path, pattern = "SSURef_NR99-138.1-dna-seqs-derep-uniq-1389f-1510r-db_processed.rds", full.names = T)


db2_f <- list.files(path = path, pattern = "pr2_version_4.14.0_SSU-1389f-1510r-db_processed.rds", full.names = T)

omit_cols <- c("Feature ID","dna_width")

db1 <- read_rds(db1_f) %>% select_at(vars(-contains(omit_cols)))

db2 <- read_rds(db2_f)  %>% select_at(vars(-contains(omit_cols)))


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


# 2) Return DB groups ----

taxa_list <- c(NON_EUK, PHYLA)

find_obs_tax <- function(filter_list = ..., tax_df, ...) {
  
  # tax_df %>% filter_at(vars(into), any_vars(grepl(w_filter, .)))
  
  filter_list <- str_to_lower(filter_list)
  
  w_filter <- sort(filter_list)
  
  w_filter <- paste(w_filter, collapse = "|")
  
  into <- tax_df %>% select_if(is.character) %>% names()
  
  # 1) Find obs taxon groups
  
  tax_df %>% 
    pivot_longer(cols = all_of(into), values_to = "Taxon") %>%
    mutate(Taxon = str_to_lower(Taxon)) %>%
    filter(grepl(w_filter, Taxon)) %>%
    count(name, Taxon) -> taxon_obs_df
  
  
  taxon_obs_df %>% pull(Taxon) -> taxon_obs
  
  # 2) Which taxon groups are not present in the df (not working yet)
  
  # taxon_obs <- sort(taxon_obs)
  
  # taxon_obs <- paste(taxon_obs, collapse = "|")
  
  # taxon_non_obs <- filter_list[grepl(w_filter, filter_list)]
  
  taxon_non_obs <- taxon_obs[!grepl(w_filter, taxon_obs)]
  
  taxon_non_obs <- ifelse(length(taxon_non_obs) == 0, "none", taxon_non_obs) 
  
  taxon_non_obs_df <- data.frame(name = "none", Taxon = taxon_non_obs, n = 0)
  
  out <- rbind(taxon_obs_df, taxon_non_obs_df) #%>%
    # filter(grepl(w_filter, Taxon))
  
  out$Taxon <- str_to_sentence(out$Taxon)
  
  return(out)
  
}


# Test:
# db <- db1 %>% filter(!grepl("Alveolata",rank2))
# find_obs_tax(c(NON_EUK, PHYLA), db) %>% view()

df1 <- find_obs_tax(c(NON_EUK, PHYLA), db1) %>% mutate(DB = "V9-SSURef_NR99-138.1")
df2 <- find_obs_tax(c(NON_EUK, PHYLA), db2) %>% mutate(DB = "V9-pr2_version_4.14.0_SSU")

rbind(df1, df2) %>% 
  select(-name) %>%
  pivot_wider(names_from = "DB", values_from = "n") %>% 
  view()

# 3) Filtering non assigned ASVS =====


wd <- "~/Documents/MEIOFAUNA_PAPER/INPUTS/"

tax_df <- read_tsv(list.files(path = wd, pattern = 'taxonomy.tsv$', full.names = T))

dna_f <- list.files(path = wd, pattern = "dna-sequences.fasta", full.names = T)

table(sapply(strsplit(tax_df$Taxon, ";"), `[`, 1))

# Unassigned: 1
# d__Eukaryota: 4825 

into <- paste0("rank", 1:9)

tax_df <- tax_df %>%
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(vars(into), funs(str_replace_all(., c("os__" = "", "[a-z]__" = "", 
    "_sp."=NA_character_,
    "Incertae_Sedis" = NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_)))) %>%
  mutate_at(vars(into),  ~na_if(., ''))


# hist(tax_df$Confidence)

# 3.1) Retrive taxonomic interest groups =====


tax_df %>% select_at(vars(contains("rank"))) %>%
  # distinct() %>%
  find_obs_tax(c(NON_EUK, PHYLA), .) %>%
  summarise(sum(n)) 


#  3,146 from a set of 3,273)
# A set of 3272 from 4,826 ASVS were taxon-specific classfied using SILVA 138.1 V9 DB


# 3.2) Remove contaminants source: ----

# In addition, add prevalence and abundance

Freq <- function(x){rowSums(x > 1)}

read_tsv(list.files(path = wd, pattern = 'table_100_80', full.names = T), skip = 1) %>%
  mutate(
    TotalAbundance = rowSums(across(where(is.double))),
    Prevalence = Freq(across(where(is.double)))) %>% 
  select(`Feature ID`, TotalAbundance, Prevalence) %>%
  arrange(desc(Prevalence)) %>% # pull(TotalAbundance) %>% sum()
  filter(!`Feature ID` %in% IDS) -> ab_f

tReads <- sum(ab_f$TotalAbundance)

ab_f <- ab_f %>% mutate(pct_ab = TotalAbundance/tReads)

# tax_df2 %>% left_join(ab_f) %>% view()

# Mammalia: Gorilla_gorilla, Homo_sapiens, Sus_scrofa

contaminant_source <- c("Mammalia", "Aves")

contaminant_source <- paste(contaminant_source, collapse = "|")

contaminant_df <- tax_df %>% 
  filter_at(vars(into), any_vars(grepl(contaminant_source, .))) %>%
  left_join(ab_f)

path_out <- path

file_out <- "/contaminant_sources.csv"

write_excel_csv(contaminant_df, file = paste0(path, file_out))

contaminant_ids <- contaminant_df %>% pull(`Feature ID` )

# 3.2.1) RUN THIS STEP ----

tax_df <- tax_df %>% filter(!`Feature ID` %in% contaminant_ids)

# 3.3) Keep dna and taxa ---- 

taxa_list <- c(NON_EUK, PHYLA)

taxa_list <- paste(taxa_list, collapse = "|")

taxa_list <- str_to_lower(taxa_list)

# the number of unique ASVS is 3,146 (from a set of 3,273)

tax_df %>% 
  mutate_at(vars(into), funs(str_to_lower(.))) %>%
  filter_at(vars(into), any_vars(grepl(taxa_list, .))) %>% 
  pull(`Feature ID`) -> IDS

str(IDS) # 3,146

KEEP_IDS <- tax_df %>% filter(!`Feature ID` %in% IDS) %>% pull(`Feature ID`)


# Lets retrive non assigned asvs and use PR2 as new DB

# 3.2) keep dna ----
# Reverse matching to retrive non classified

dna <- Biostrings::readDNAStringSet(dna_f) # 4,826 ASVS

sum(keep <- names(dna) %in% KEEP_IDS) # 1,628

length(dna <- dna[keep]) # 1,628

# 1680+3146 = 4,826

path_out <- path

file_out <- "/non-assigned-SILVA-dna-sequences.fasta"

Biostrings::writeXStringSet(dna, filepath =  paste0(path_out, file_out), format = "fasta")

# What about other grups -----


tax_df2 <- tax_df %>% 
  filter(`Feature ID` %in% KEEP_IDS)




which_rank <- "rank5"

rank_pos <- which(into %in% which_rank)

lineajedf <- tax_df2 %>% select(any_of(into[1:rank_pos])) %>% 
  rename("Taxon"= which_rank) %>%
  distinct(Taxon, .keep_all = T) 

.tax_subset <- tax_df2 %>% rename("Taxon"= which_rank) 

.tax_subset %>%
  left_join(ab_f) %>%
  group_by(Taxon) %>%
  summarise(TotalAbundance = sum(TotalAbundance), 
    min_freq = min(Prevalence), max_freq = max(Prevalence),
    n_asvs = n(), pct_ab = sum(pct_ab)) %>% arrange(desc(pct_ab)) %>%
  left_join(lineajedf) %>%
  view()


tax_df2 %>% left_join(ab_f) %>% 
  ggplot(aes(x = log10(TotalAbundance), y = Prevalence, color = rank2)) + 
  geom_point()


tax_df2 %>% select(into) -> tax

# 3.x) add Resolution (no posible, debido a los gaps , )

# taxa must be turned as no-rank-propagation

# tax <- tax[complete.cases(tax[2:5]),]

head(res <- !is.na(tax))

tail(res <- rowSums(res)) # max.rank - 

features <- c(nrow(tax) - table(res))

features

# MUST BE CUMMULATIVE DIVISION

pct <- c(features / nrow(tax))

# 5   6   7   8 # <- 1 ASV NO ESTA CLASIFICADO EN EL PRIMER POSICION 
# 648 584 391 324

caption <- "V9-SSURef_NR99-138.1"

df1 <- data.frame(into = into[1:length(pct)], features, pct, g = caption)


df1 %>%
  # mutate(into = factor(into, levels = into)) %>%
  ggplot(aes(x = into, y = features)) +
  geom_path(size = 1.5, alpha=0.6, group = 1) +
  geom_point(size = 3, alpha=0.6) +
  geom_text(aes(label = paste0(round(pct*100, digits = 2), "%")), 
    size = 4, vjust = -1, family = "GillSans") +
  scale_y_continuous(labels = scales::comma) +
  labs(y = "Number of features", x = '', subtitle = 'Features with taxonomic assignation', caption = caption) +
  theme_bw(base_size = 20, base_family = "GillSans") -> ps 

ps + theme(panel.border = element_blank()) -> ps

ps

# 


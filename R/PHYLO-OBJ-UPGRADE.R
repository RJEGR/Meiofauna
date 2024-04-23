# READ TAX AND PROCESS
# READ COUNT AND MERGE TO TAX
# PREPARE PHYLO OBJECT

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/classify-consensus-blast_dir"

# 1

tax <- list.files(path = wd, pattern = "taxonomy.tsv", full.names = T)

dim(tax <- read_tsv(tax)) # 9205

# hist(tax$Consensus)

dim(tax <- filter(tax, Taxon != "Unassigned")) # 3979

# hist(tax$Consensus)

into <- c("k", "p", "c", "o","f", "g", "s")

mutate_ranks <- c("none" = NA_character_, "os__" = "", "[a-z]__" = "", 
                                            "Incertae_Sedis" = NA_character_, 
                                            "uncultured" = NA_character_, 
                                            "Unassigned"=NA_character_,
                                            "Unknown" = NA_character_)
tax <- tax %>% 
  separate(Taxon, sep = ";", into = into) %>%
  mutate_all(list(~ str_replace_all(., mutate_ranks))) %>%
  mutate_all(function(x) {na_if(x,"")}) %>%
  data.frame(row.names = .$`Feature ID`)
  
str(which_asvs_classified <- tax$Feature.ID) # length of 3979

# 2) BIND W/ ABUNDANCE

# wd <- "C:/Users/Israel V/Documents/MEIOFAUNA/raw-seqs/"

# wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/MULTIRUN-ANALYSIS/LIBRARY/"

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/"


ab_f <- list.files(path = wd, pattern = 'multirun_ASVs_count.table', full.names = T)

ab <- read.delim(ab_f, sep = "\t") # %>% as_tibble(rownames = "Feature ID")

# rownames(ab)[ab[colnames(ab) %in% "PE1.F4"] > 0]

sum(keep <- rownames(ab) %in% which_asvs_classified) # must match 3979

dim(ab <- ab[keep,]) # and 163 samples

.all_sam <- names(ab)


# KEEP PIVOTAL DATA

ab %>%
  as_tibble(rownames = "Feature ID") %>%
  # filter(`Feature ID` %in% which_asvs_classified) %>%
  pivot_longer(-`Feature ID`, names_to = "KEYID", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  group_by(KEYID) %>%
  summarise(n_asvs = n(), TotalAbundance = sum(abundance)) %>%
  mutate(LIBRARY_ID = gsub("[.]", "-", KEYID)) -> diversity_db


  
ab %>%
  as_tibble(rownames = "Feature ID") %>%
  # filter(`Feature ID` %in% which_asvs_classified) %>%
  pivot_longer(-`Feature ID`, names_to = "KEYID", values_to = "abundance") %>%
  filter(abundance > 1) %>%
  mutate(Cruise = sapply(strsplit(KEYID, "[.]"), `[`, 1)) %>%
  group_by(Cruise) %>%
  summarise(n_asvs = length(unique(`Feature ID`)),
    n_sam = length(unique(KEYID)),
    TotalAbundance = sum(abundance))

# # BY ASV
# Freq <- function(x){rowSums(x > 1)}
# # 
# ab %>%
#   as_tibble(rownames = "Feature ID") %>% 
#   mutate(
#     TotalAbundance = rowSums(across(where(is.integer))),
#     Prevalence = Freq(across(where(is.integer)))) %>%
#   select(`Feature ID`, TotalAbundance, Prevalence) %>%
#   arrange(desc(Prevalence)) %>% view()
#   filter(`Feature ID` %in% which_asvs_classified) # -> diversity_db

# JOIN TO diversity_db to ambiental parameter db


DB <- read_tsv(list.files(path = wd, pattern = "MANIFEST.tsv", full.names = T))  

MTD <- DB %>% right_join(diversity_db, by = "LIBRARY_ID") %>% 
  data.frame(row.names = .$LIBRARY_ID)

sum(keep <- gsub("[.]", "-", colnames(ab)) %in% MTD$LIBRARY_ID) # match 160 cols

dim(ab <- ab[,keep])

# match sort
sorted_sam <- sort(gsub("[.]", "-", colnames(ab)))

identical( sorted_sam, sort(MTD$LIBRARY_ID))

colnames(ab) <- gsub("[.]", "-", colnames(ab))


MTD %>% count(Region)

# view(MTD)
# 
# fill_sample_gaps <- MTD$LIBRARY_ID[!MTD$LIBRARY_ID %in% DB$LIBRARY_ID]

# MTD %>% count(Region)
# DB$LIBRARY_ID[!DB$LIBRARY_ID %in% MTD$LIBRARY_ID]

identical(sort(names(ab)),sort(rownames(MTD)))

identical(sort(rownames(ab)), sort(rownames(tax)))

library(phyloseq)

phyloseq = phyloseq(otu_table(ab, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')),
                    sample_data(MTD))

saveRDS(phyloseq, paste0(wd, "phyloseq.rds"))







rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/Documents/MEIOFAUNA_PAPER/INPUTS/"
subdir <- paste0(wd, "CURATED_DB_DIR")
# library(Biostrings)
# library(rstatix)

library(tidyverse)

tax_f <- read_tsv(list.files(path = subdir, pattern = 'taxonomy.tsv$', full.names = T))

ab_f <-  read_tsv(list.files(path = wd, pattern = 'table_100_80', full.names = T), skip = 1)


# CLEAN TAXONOMY
# 

into <- paste0("rank", 0:9)

# Res <- function(x) { rowSums(!is.na(x)) }

# "Unassigned","Incertae_Sedis", "uncultured", "Clade_*"

tax_f <- tax_f %>%
  select(-`Feature ID` , -Confidence) %>%
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(into, funs(str_replace_all(., c("os__" = "", "[a-z]__" = "", 
    "Incertae_Sedis" = NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_))))

# "uncultured" = NA_character_
# "Incertae Sedis"=NA_character_
# mutate(Resolution = Res(across(where(is.double)))) %>%

# RESOLUTION PER TAXON ----
# CUANTOS TAXONES DIFERENTES POR RANK

tax_f %>% select(into) -> tax
nt <- function(x) {length(na.omit(unique(x)))}
apply(tax[complete.cases(tax),], 2, nt)

max.rank <-  ncol(tax) 

tail(!is.na(tax))

# , "Phylum", "Class", "Order", "Family"
#"Kingdom","Clade","Subclade","","Genus","subgroup","Specie"

tail(res <- rowSums(!is.na(tax))) # max.rank - 

features <- c(nrow(tax) - table(res))

pct <- c(features / nrow(tax))

caption <- "V9-SSURef_NR99-138.1"

df1 <- data.frame(into, features, pct, g = caption)

df1 %>%
  mutate(into = factor(into, levels = into)) %>%
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


# 2) ----
# Which asvs does not assign to V9-SSURef_NR99-138.1 as reference
# corte a Class


tax %>% count(k)

tax %>% 
  pivot_longer(cols = all_of(into), values_to = "Taxon") %>%
  # filter(name == "rank1") %>%
  drop_na(Taxon) %>%
  count(name, Taxon) %>% view()


# 3) -----
# how many taxa-groups of interest

c(NON_EUK, PHYLA)

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
  
  # 2) Which taxon groups are not present in the df
  
  taxon_obs <- sort(taxon_obs)
  
  taxon_obs <- paste(taxon_obs, collapse = "|")
  
  taxon_non_obs <- filter_list[!grepl(taxon_obs, filter_list)] 
  
  taxon_non_obs_df <- data.frame(name = "none", Taxon = taxon_non_obs, n = 0)
  
  out <- rbind(taxon_obs_df, taxon_non_obs_df) %>%
    filter(grepl(w_filter, Taxon))
  
  out$Taxon <- str_to_sentence(out$Taxon)
  
  return(out)
  
}

find_obs_tax(c(NON_EUK, PHYLA), tax_f) %>% view()

# add comparsing w/ previous

# into <- c("domain", "kingdom", "phylum", "class", "order", "suborder", "family", "genus", "specie")

into <- paste0("rank", 1:9)

taxprev <- read_tsv(list.files(path = wd, pattern = 'taxonomy.tsv.bkp', full.names = T)) %>%
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(into, funs(str_replace_all(., c("D_[0-9]__" = "", "D_1[0-9]__" = "", "Incertae Sedis"=NA_character_, "uncultured" = NA_character_)))) %>%
  select(into)

caption <- "SILVA 132 18S SSU"

tail(resprev <- rowSums(!is.na(taxprev))) # max.rank - 

featuresprev <- c(nrow(taxprev) - table(resprev))

pctprev <- c(featuresprev / nrow(taxprev))

data.frame(into, featuresprev, pctprev)

df2 <- data.frame(into, featuresprev, pctprev, g = caption)
names(df2) <- names(df1)


rbind(df1, df2) %>%
  mutate(into = factor(into, levels = paste0("rank", 1:9))) %>%
  ggplot(aes(x = into, y = features, color = g, group = g)) +
  geom_path(size = 1.5, alpha=0.6) +
  geom_point(size = 3, alpha=0.6) +
  geom_text(aes(label = paste0(round(pct*100, digits = 2), "%")), 
    size = 4, vjust = -1, family = "GillSans") +
  scale_y_continuous(labels = scales::comma) +
  labs(y = "Number of features", x = '', subtitle = 'ASVs with taxonomic assignation') +
  theme_bw(base_size = 20, base_family = "GillSans") -> ps 

ps + theme(panel.border = element_blank(), legend.position = "top") +
  guides(color=guide_legend("",nrow=1)) -> ps

ps
  


# (OMIT ) fill ----
# ranks <- c("rank2","rank3","rank4",'rank5', 'rank6', 'rank7')

tax_f %>%
  select_at(into) %>% 
  mutate(id = 1:nrow(tax_f)) %>%
  pivot_longer(cols = all_of(into)) %>% fill(value) %>% # con esto rellenas espacios vacios o NA !!!
  pivot_wider(names_from = name) %>%
  mutate(`Feature ID` = tax_f$`Feature ID`) %>%
  select(-id) -> tax

# tax %>% left_join(tax_f %>% select(!into), by = "Feature ID") -> tax_f

# PREVALENCE -----
# FILTER BASED ON PREVALENCE

# Using the table above, determine the taxa to filter based on the 0.01 threshold which means threshold*100 of the total ab

tax_f %>% mutate(Level = rank4) -> tax_f

threshold <- 0.001 # el 1% de los reads

sum(tax_f$TotalAbundance)*threshold

table(tax_f$TotalAbundance/sum(tax_f$TotalAbundance) >= threshold)

keepTaxa <- unique(tax_f$Level[tax_f$TotalAbundance/sum(tax_f$TotalAbundance) >= threshold])

prevalence_df <- tax_f[tax_f$Level %in% keepTaxa,]

prevalence_df %>% 
  # filter(TotalAbundance > 1) %>%
  filter(Level %in% keepTaxa) %>%
  group_by(Level) %>%
  rstatix::cor_test(Prevalence, TotalAbundance)

# ggplot(tax_f, aes(TotalAbundance, Prevalence)) + geom_point() + scale_x_log10()


tax_f %>% 
  drop_na(Level) %>% 
  # mutate(col = ifelse(Prevalence == 1, 'Exclusive', 'Intersected')) %>%
  mutate(Level = ifelse(Level %in% keepTaxa, Level, "Other")) %>%
  mutate(Level = factor(Level, levels = c(keepTaxa, "Other"))) %>%
  mutate(Prevalence = Prevalence/nsamples) %>%
  ggplot(aes(TotalAbundance, Prevalence)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_hline(yintercept = threshold, alpha = 0.5, linetype = 2) +
  # geom_vline(xintercept = 0.1, alpha = 0.5, linetype = 2) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  # scale_color_manual(values = c('red', 'blue')) +
  labs(y = "Prevalence (Frac. Samples)", x = "Total Abundance (log10)", color = "Family\nResolution") +
  facet_wrap(~Level, nrow = 1) +
  theme_bw(base_size = 17, base_family = "GillSans") +
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12)) -> p2

p2 + theme(panel.border = element_blank()) 



# tax_f %>% view()
# la esperanza es que todos los asvs/otus lleguen a la clasificacion mas profunda (species), evaluemos cuanto fue asi:

# divergence (omitir esta prueba)

length(unique(tax_f$rank2))

tax_f %>%
  # mutate(Level = rank2) %>%
  drop_na(Level) %>%
  group_by(Level) %>%
  summarise(
    n = n(),
    nR = length(unique(rank7)),
    divergence = 1-(1/length(unique(rank7))),
    nTotal = sum(TotalAbundance)) %>%
  mutate(pct = nTotal/ sum(nTotal) * 100) %>%
  arrange(n) %>%
  mutate(cs = cumsum(n), csR = cumsum(nR)) -> dataV 

dataV %>%
  # filter(divergence > 0.5) %>% # pull(Phylum)
  rstatix::cor_test(n, nTotal)

dataV %>%
  mutate(Level = ifelse(divergence > 0.5, Level, NA)) %>%
  mutate(col = ifelse(Level %in% keepTaxa, "Keep", "Drop")) %>%
  mutate(facet = ifelse(csR > 15, "A", "B")) %>%
  ggplot(aes(divergence, csR)) +
  geom_point(aes(size = pct, color = col)) +
  scale_size(name = "Abundance (%)") +
  # ggsci::scale_color_gsea() +
  ggrepel::geom_text_repel(aes(label = Level), size = 4, family = "GillSans") +
  facet_grid(facet~., scales = "free_y",space = "fixed") +
  # scale_color_manual(values = c('black', 'red')) +
  labs(x = "x", y = 'y') +
  scale_color_manual(values = c('red', 'blue')) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()) -> p1

p1

# Opisthokonta son Fungi y Metazoa en rank4

tax_f %>% group_by(rank4) %>%
  summarise(
    nTotal = sum(TotalAbundance))

# THEN
# TO REMOVE
# rank5 =
# Arachnida (insectos), Tetrapoda (mamiferos),


tax %>% distinct(rank2)

nt <- function(x) {length(na.omit(unique(x)))}

x <- tax[complete.cases(tax),]

apply(tax[complete.cases(tax),], 2, nt)

tax %>% drop_na(rank7) %>% distinct(rank7)

source('~/Documents/GitHub/zooplankton_gom/names2worms.R')

# seleccionamos solo los filos de zooplancton (~14)
# verificamos cuales filos son zoo

query <- pull(unique(taxonomy[, rank.names[4]]))
query <- query[!is.na(query)]
names_ <- names2worms_(query)

names_ %>% filter(Kingdom_wm %in% 'Animalia') %>% 
  drop_na(Phylum_wm) %>% distinct(ori_name) %>%
  pull() -> marineMetazoa

# nsamples <- ab_f %>% select_if(is.double) %>% ncol()


# ab <- as.data.frame(ab_f)
# rownames(ab) <- ab_f$`#OTU ID`
# ab$`#OTU ID` <- NULL
# TotalAbundance <- rowSums(ab)

Freq <- function(x){rowSums(x > 1)}

ab_f %>%
  mutate(
    TotalAbundance = rowSums(across(where(is.double))),
    Prevalence = Freq(across(where(is.double)))) %>% 
  select(`Feature ID`, TotalAbundance, Prevalence) %>%
  arrange(desc(Prevalence)) -> ab


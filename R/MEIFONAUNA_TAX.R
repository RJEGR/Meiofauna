
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/Documents/MEIOFAUNA_PAPER/INPUTS/"

# library(Biostrings)
# library(rstatix)

library(tidyverse)

tax_f <- read_tsv(list.files(path = wd, pattern = 'taxonomy.tsv', full.names = T))
ab_f <-  read_tsv(list.files(path = wd, pattern = 'table_100_80', full.names = T), skip = 1)

nsamples <- ab_f %>% select_if(is.double) %>% ncol()

into <- paste0("rank", 1:7)

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


# CLEAN TAXONOMY
# 

Res <- function(x) { rowSums(!is.na(x)) }

tax_f %>%
  left_join(ab) %>%
  mutate(Resolution = Res(across(where(is.double)))) %>%
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(into, funs(str_replace_all(., c("D_[0-9]__" = "", "D_1[0-9]__" = "", "Incertae Sedis"=NA_character_)))) -> tax_f


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

data.frame(into, features, pct) %>%
  mutate(into = factor(into, levels = into)) %>%
  ggplot(aes(x = into, y = features)) +
  geom_path(size = 1.5, alpha=0.6, group = 1) +
  geom_point(size = 3, alpha=0.6) +
  geom_text(aes(label = paste0(round(pct*100, digits = 2), "%")), 
    size = 4, vjust = -1, family = "GillSans") +
  labs(y = "Number of features", x = '', subtitle = 'Features with taxonomic assignation') +
  theme_bw(base_size = 20, base_family = "GillSans") -> ps 

ps + theme(panel.border = element_blank()) -> ps

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

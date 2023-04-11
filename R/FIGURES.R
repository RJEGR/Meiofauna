# DATAVIZ

library(tidyverse)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# 0) LOAD DATA =====

wd <- "~/Documents/MEIOFAUNA_PAPER/INPUTS/"

# subdir <- paste0(wd, "CURATED_DB_DIR")


subdir <- paste0(wd, "classify-consensus-blast_dir")

subdir2 <- paste0(wd, "BKP/classify-consensus-blast_dir/")

caption_1 <- "A) Enriched"
caption_2 <- "B) SSURef_NR99_138"# "V9-SSURef_NR99-138.1"
caption_3 <- "C) SSURef_NR99_132"# "SILVA 132 18S SSU"



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


taxa_list <- c(NON_EUK, PHYLA)
taxa_list <- paste(taxa_list, collapse = "|")
taxa_list <- str_to_lower(taxa_list)



tax_f <- read_tsv(list.files(path = subdir, pattern = 'taxonomy.tsv$', full.names = T)) %>% mutate(DB = caption_1)
tax_f_1 <- read_tsv(list.files(path = subdir2, pattern = 'taxonomy.tsv$', full.names = T)) %>% mutate(DB = caption_2)
tax_f_2 <- read_tsv(list.files(path = wd, pattern = 'taxonomy.tsv.bkp$', full.names = T)) %>% mutate(DB = caption_3)


contaminant_source_f <- "~/Documents/MEIOFAUNA_PAPER/DB_COMPARISON/contaminant_sources.csv"
contaminant_source <- read_csv(contaminant_source_f, col_names = T) %>% distinct(`Feature ID`) %>% pull()

# 1) SPLIT TAX ======

# tax_f %>% view()

into <- c("k","p", "c", "o", "f", "g", "s")

tax1 <- tax_f %>%
  # select(-`Feature ID` , -Consensus) %>%
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(vars(all_of(into)), list(~ str_replace_all(., c("[a-z]__" = "", "Unassigned"=NA_character_)))) %>%
  mutate_at(vars(all_of(into)), list(~ gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(all_of(into)),  ~na_if(., ""))

into2 <- paste0("rank", 1:9)

tax2 <- tax_f_1 %>%
  separate(Taxon, sep = ";", into = into2) %>% 
  mutate_at(vars(all_of(into2)), list(~str_replace_all(., c("os__" = "", "[a-z]__" = "", "Incertae Sedis"=NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_, "_sp."=NA_character_, "metagenome"=NA_character_, 
    "environmental"=NA_character_, "Group|group"=NA_character_)))) %>%
  mutate_at(vars(all_of(into2)), list(~ gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(all_of(into2)),  ~na_if(., ""))


  
tax3 <- tax_f_2 %>%
  separate(Taxon, sep = ";", into = into2) %>% 
  mutate_at(vars(all_of(into2)), list(~str_replace_all(., c("D_[0-9]__" = "", "D_1[0-9]__" = "", "Incertae Sedis"=NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_, "_sp.|sp."=NA_character_, "metagenome"=NA_character_,
    "environmental"=NA_character_, "Group|group"=NA_character_)))) %>%
  # mutate_at(vars(all_of(into2)), list(~ gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(all_of(into2)),  ~ na_if(., ""))

# FILTER BY =====

# TEORICAMENTE, LOS GRUPOS DE taxa_list son los unicos presentes en el objeto tax1
# tax1 <- tax1 %>% 
#   mutate_at(vars(into), list(~ str_to_lower(.))) %>%
#   filter_at(vars(into), any_vars(grepl(taxa_list, .))) %>%
#   mutate_at(vars(all_of(into)), list(~ str_to_sentence(.))) # %>%
  # filter(!`Feature ID` %in% contaminant_source)

tax2 <- tax2 %>% 
  mutate_at(vars(into2), list(~ str_to_lower(.))) %>%
  filter_at(vars(into2), any_vars(grepl(taxa_list, .))) %>%
  mutate_at(vars(all_of(into2)), list(~ str_to_sentence(.))) %>% 
  filter(!`Feature ID` %in% contaminant_source) 
  

tax3 <- tax3 %>% 
  mutate_at(vars(into2), list(~ str_to_lower(.))) %>%
  filter_at(vars(into2), any_vars(grepl(taxa_list, .))) %>% 
  mutate_at(vars(all_of(into2)), list(~ str_to_sentence(.))) %>%
  filter(!`Feature ID` %in% contaminant_source)


# tax2 %>% view()
# tax3 %>% view()

# VIZ =====

rbind(select(tax1, Consensus, DB), 
      select(tax2, Consensus, DB),
      select(tax3, Consensus, DB))  %>% 
  ggplot(aes(Consensus, color = DB, fill = DB)) + 
  geom_histogram(alpha=0.5, position="dodge",bins = 50) +
  # geom_rug() +
  # scale_color_brewer(palette="Dark2")  +
  # scale_fill_brewer(palette="Dark2") +
  theme_bw(base_size = 20, base_family = "GillSans") +
  theme(panel.border = element_blank(), legend.position = "top") +
  labs(y = "# ASVs")

# 2)

smmy1 <- tax1 %>% 
  pivot_longer(cols = all_of(into), values_to = "Taxon") %>%
  drop_na(Taxon) %>%
  count(name) %>% arrange(match(name, into)) %>% mutate(pct = n/nrow(tax1), DB = caption_1) %>%
  mutate(name = factor(name, levels = into))

smmy2 <- tax2 %>% 
  pivot_longer(cols = all_of(into2), values_to = "Taxon") %>% 
  drop_na(Taxon) %>% 
  count(name) %>% arrange(match(name, into2)) %>% mutate(pct = n/nrow(tax2), DB = caption_2)


smmy3 <- tax3 %>% 
  pivot_longer(cols = all_of(into2), values_to = "Taxon") %>%
  drop_na(Taxon) %>% 
  count(name) %>% arrange(match(name, into2)) %>% mutate(pct = n/nrow(tax3), DB = caption_3)


# ESTA IMAGEN NO ES TAN COMPARABLE UNA DB VS LA OTRA
rbind(smmy1) %>% # , smmy2, smmy3
  ggplot(aes(x = name, y = n, color = DB, group = DB)) +
  facet_wrap( DB ~., scales = 'free_x') +
  geom_path(size = 1.5, alpha=0.6) +
  geom_point(size = 3, alpha=0.6) +
  geom_text(aes(label = paste0(round(pct*100, digits = 2), "%")), 
    size = 4, vjust = -1, family = "GillSans") +
  scale_y_continuous(labels = scales::comma) +
  theme_bw(base_size = 20, base_family = "GillSans") +
  theme(panel.border = element_blank(), legend.position = "top")


# SPLIT TAX BY CLADES ====

METAZOAN <- str_to_lower(paste(PHYLA, collapse = "|"))

UNI_EUK <- str_to_lower(paste(NON_EUK, collapse = "|"))

METAZOAN_DF <- tax1 %>%
  mutate_at(vars(into), list(~ str_to_lower(.))) %>%
  filter_at(vars(into), any_vars(grepl(METAZOAN, .))) %>%
  mutate_at(vars(into), list(~ str_to_sentence(.))) %>%
  mutate(Clade = "Metazoans")


UNI_EUK_DF <- tax1 %>%
  mutate_at(vars(into), list(~ str_to_lower(.))) %>%
  filter_at(vars(into), any_vars(grepl(UNI_EUK, .))) %>%
  mutate_at(vars(into), list(~ str_to_sentence(.))) %>%
  mutate(Clade = "Unicellular eukaryotes")


tax1 <- rbind(METAZOAN_DF,UNI_EUK_DF)

write_tsv(tax1, file = paste0(subdir, "/CURATED-classify-consensus-blast-asvs.tsv"))

# 


METAZOAN_DF <- tax2 %>%
  mutate_at(vars(into2), list(~ str_to_lower(.))) %>%
  filter_at(vars(into2), any_vars(grepl(METAZOAN, .))) %>%
  mutate_at(vars(into2), list(~ str_to_sentence(.))) %>%
  mutate(Clade = "Metazoans")


UNI_EUK_DF <- tax2 %>%
  mutate_at(vars(into2), list(~ str_to_lower(.))) %>%
  filter_at(vars(into2), any_vars(grepl(UNI_EUK, .))) %>%
  mutate_at(vars(into2), list(~ str_to_sentence(.))) %>%
  mutate(Clade = "Unicellular eukaryotes")


tax2 <- rbind(METAZOAN_DF,UNI_EUK_DF)

# tax 3


METAZOAN_DF <- tax3 %>%
  mutate_at(vars(into2), list(~ str_to_lower(.))) %>%
  filter_at(vars(into2), any_vars(grepl(METAZOAN, .))) %>%
  mutate_at(vars(into2), list(~ str_to_sentence(.))) %>%
  mutate(Clade = "Metazoans")


UNI_EUK_DF <- tax3 %>%
  mutate_at(vars(into2), list(~ str_to_lower(.))) %>%
  filter_at(vars(into2), any_vars(grepl(UNI_EUK, .))) %>%
  mutate_at(vars(into2), list(~ str_to_sentence(.))) %>%
  mutate(Clade = "Unicellular eukaryotes")


tax3 <- rbind(METAZOAN_DF,UNI_EUK_DF)


clade_1 <- tax1 %>% count(Clade) %>% mutate(DB = caption_1, pct = n/sum(n))
clade_2 <- tax2 %>% count(Clade) %>% mutate(DB = caption_2, pct = n/sum(n))
clade_3 <- tax3 %>% count(Clade) %>% mutate(DB = caption_3, pct = n/sum(n))

rbind(clade_1, clade_2, clade_3) %>% view()
  mutate(DB = factor(DB, levels = unique(DB))) %>%
  ggplot(aes(x = pct, y = DB, fill = Clade)) +
  geom_col() +
  # geom_col(position = "dodge") +
  scale_x_continuous(labels = scales::percent) +
  theme_bw(base_size = 20, base_family = "GillSans") +
  theme(panel.border = element_blank(), legend.position = "top") +
  labs(y = "", x = "% ASVs") +
  ggsci::scale_fill_uchicago()

# PLOTS =======

ab_f <-  read_tsv(list.files(path = wd, pattern = 'table_100_80', full.names = T), skip = 1)

MTD <-  read_tsv(list.files(path = wd, pattern = "mapping-file-corregido.tsv", full.names = T)) %>% select(`#SampleID`, Depth, Region) %>%
  mutate(Region = factor(Region, levels = c("Yucatan", "NW Shelf", "NW Slope", "Deep-sea")))


# PREVALENCE ====
Freq <- function(x){rowSums(x > 1)}

ab_f %>%
  mutate(
    TotalAbundance = rowSums(across(where(is.double))),
    Prevalence = Freq(across(where(is.double)))) %>% 
  select(`Feature ID`, TotalAbundance, Prevalence) %>%
  arrange(desc(Prevalence)) -> ab

out <- tax1 %>% left_join(ab)

# out %>% view()

n_sam <- ab_f %>% select_if(is.double) %>% ncol()

out %>% 
  # drop_na(k) %>%
  ggplot(aes(TotalAbundance, Prevalence/n_sam)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  # scale_color_manual(values = c('red', 'blue')) +
  labs(y = "Prevalence (Frac. Samples)", x = "Total Abundance (log10)", color = "Consensus", caption = caption_1) +
  facet_grid(Clade ~ k) +
  theme_bw(base_size = 17, base_family = "GillSans") +
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12)) -> p2

p2 + theme(panel.border = element_blank(), legend.position = "top")


out %>%
  group_by(Clade, p) %>%
  summarise(TotalAbundance = sum(TotalAbundance)) %>% view()


# BARPLOT =======

which_sam <- ab_f %>% select_if(is.double) %>% names()

# agglom_lev <- "k"
# low_ab <- 0.01

df <- ab_f %>% 
  pivot_longer(cols = all_of(which_sam), 
    values_to = "ab", 
    names_to = "#SampleID") %>%
  filter(ab > 0) %>%
  left_join(MTD) %>%
  right_join(tax1) 


df %>% distinct(p)

barTax <- function(df, agglom_lev = "k", low_ab = 0.01) {

  library(ggsci)
  
  df <- df %>%
    rename( "Level" = agglom_lev) %>%
    drop_na(Level) %>%
    group_by(`#SampleID`) %>%
    mutate(RA = ab/sum(ab)) %>%
    mutate(Level = ifelse(RA < low_ab, "ZLow", Level))
  
  
  labels <- df %>% pull(Level) %>% unique() %>% sort()
  colourCount = length(labels)
  
  if(colourCount > 7) {
    library(ggsci)
    getPalette <- colorRampPalette(pal_uchicago()(7))(colourCount)
  } else {
    getPalette <- pal_uchicago("default")(colourCount)
  }
  
  # pal_locuszoom
  # 
  
  getPalette[length(getPalette)] <- "Black"
  labels[length(labels)] <- "Low abundance"
  
  
  
  df %>%
    ggplot(aes(x = RA, y = `#SampleID`, fill = Level)) +
    facet_grid( Region ~., scales = "free", space = "free", switch = "y") +
    geom_col() +
    scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_fill_manual(agglom_lev, labels = labels, values = getPalette) +
    labs(y = "", x = "Relative abundance (%)") +
    theme_classic(base_size = 16, base_family = "GillSans") +
    theme(axis.ticks.y = element_blank(), 
      axis.text.y = element_blank(), axis.line.y = element_blank(),
      strip.background = element_rect(fill = 'grey', color = 'white')) +
    guides(fill = guide_legend(""))
  
}

ggsavepath <- paste0(wd, '/Figures/')


ps <- barTax(df, agglom_lev = "k", low_ab = 0.01)

ggsave(ps, path = ggsavepath, filename = 'barplot_k.png', device = png, 
  width = 7, height = 8)

ps <-  barTax(df, agglom_lev = "p", low_ab = 0.1)


ggsave(ps, path = ggsavepath, filename = 'barplot_p.png', device = png, 
  width = 7, height = 8)

# 2) Heatmap of p =====

agglom_lev <- "p"

heatmap_df <- ab_f %>% 
  right_join(tax1) %>%
  rename( "Level" = agglom_lev) %>%
  group_by(Level, Clade) %>% 
  summarise_at(vars(which_sam), sum) %>%
  ungroup() %>%
  drop_na(Level)

hclust_df <- data.frame(heatmap_df %>% select_at(which_sam), 
  row.names = heatmap_df$Level)

hclust <- hclust(dist(hclust_df), "complete")

heatmap_df <- heatmap_df %>%
  pivot_longer(cols = which_sam, names_to = '#SampleID', values_to = 'ab') %>%
  left_join(MTD) %>%
  filter(ab > 0) %>%
  group_by(`#SampleID`) %>% 
  mutate(RA = ab/sum(ab)) 


# heatmap_df %>% summarise(sum(RA))

ps <- heatmap_df %>%
  ggplot(aes(y = Level, x = `#SampleID`, fill = RA)) +
  geom_raster() +
  facet_grid(Clade~Region, scales = "free", space = "free") +
  scale_fill_viridis_c(name = "Relative Abundance", na.value = "white", limits = c(0, 1)) +
  # ggh4x::scale_y_dendrogram(hclust = hclust, guide = ggh4x::guide_dendro(position = "left")) +
  labs(x = '', y = '') +
  guides(fill = guide_colorbar(barwidth = unit(3, "in"),
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(size = 12))) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    panel.border = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.text.x = element_blank(), axis.line.x = element_blank(),
    strip.background = element_rect(fill = 'grey', color = 'white')) 

ggsave(ps, path = ggsavepath, filename = 'heatmap_p.png', device = png, 
  width = 8, height = 6.5)

# 

color_vector <- as.character(unique(MTD$Region))

getPalette <- c("#4DAF4A", "#313695", "lightblue", "#E41A1C")


axis_col <- structure(getPalette, names = color_vector)

heatmap_df %>%
  group_by(Region, Level, Clade) %>%
  summarise(TotalAbundance = sum(ab)) %>%
  group_by(Level) %>%
  mutate(RA = TotalAbundance/sum(TotalAbundance)) %>%
  ggplot(aes(x = RA, y = Level, fill = Region)) +
  facet_grid(Clade~., scales = "free", space = "free") +
  geom_col() +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Relative abundance", y = "p") +
  scale_fill_manual(values = axis_col ) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey', color = 'white')) +
  guides(fill = guide_legend("")) -> ps

ggsave(ps, path = ggsavepath, filename = 'barplot_p_regions.png', device = png, 
  width = 8, height = 6.5)

# ABUNDANCE =====

w_cols <- ab_f %>% select_if(is.double) %>% names()

ab_f %>%
  right_join(tax1, by = "Feature ID") %>% 
  pivot_longer(cols = all_of(w_cols), names_to = "#SampleID") %>%
  filter(value > 0) %>%
  left_join(MTD ) %>%
  group_by(Region, Clade, p) %>% 
  summarise(value = sum(value)) %>% 
  # pivot_wider(names_from = p, values_from = value, values_fill = 1) %>%
  pivot_wider(names_from = Region, values_from = value, values_fill = 0) %>%
  ungroup() -> dat


ab_f %>%
  left_join(tax_f, by = "Feature ID") %>% 
  pivot_longer(cols = all_of(w_cols), names_to = "#SampleID") %>%
  filter(value > 0) %>%
  group_by(Clade) %>%
  summarise(n = length(unique(`Feature ID`)), Total = sum(value)) %>% 

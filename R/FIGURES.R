# DATAVIZ

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




tax_f <- read_tsv(list.files(path = subdir, pattern = 'taxonomy.tsv$', full.names = T)) %>% mutate(DB = caption_1)
tax_f_1 <- read_tsv(list.files(path = subdir2, pattern = 'taxonomy.tsv$', full.names = T)) %>% mutate(DB = caption_2)
tax_f_2 <- read_tsv(list.files(path = wd, pattern = 'taxonomy.tsv.bkp$', full.names = T)) %>% mutate(DB = caption_3)


rbind(tax_f, tax_f_1, tax_f_2) %>% ggplot(aes(Consensus, color = DB, fill = DB)) + 
  geom_histogram(alpha=0.5, position="dodge",bins = 50) +
  # geom_rug() +
  # scale_color_brewer(palette="Dark2")  +
  # scale_fill_brewer(palette="Dark2") +
  theme_bw(base_size = 20, base_family = "GillSans") +
  theme(panel.border = element_blank(), legend.position = "top")

# 1) SPLIT TAX ======

# tax_f %>% view()

into <- c("k","p", "c", "o", "f", "g", "s")

tax1 <- tax_f %>%
  select(-`Feature ID` , -Consensus) %>%
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(vars(all_of(into)), list(~ str_replace_all(., c("[a-z]__" = "", "Unassigned"=NA_character_)))) %>%
  mutate_at(vars(all_of(into)), list(~ gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(all_of(into)),  ~na_if(., ""))

into2 <- paste0("rank", 1:9)

tax2 <- tax_f_1 %>%
  select(-`Feature ID` , -Consensus) %>%
  separate(Taxon, sep = ";", into = into2) %>% 
  mutate_at(vars(all_of(into2)), list(~str_replace_all(., c("D_[0-9]__" = "", "D_1[0-9]__" = "", "Incertae Sedis"=NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_)))) %>%
  mutate_at(vars(all_of(into2)), list(~ gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(all_of(into2)),  ~na_if(., ""))


  
tax3 <- tax_f_2 %>%
  select(-`Feature ID` , -Consensus) %>%
  separate(Taxon, sep = ";", into = into2) %>% 
  mutate_at(vars(all_of(into2)), list(~str_replace_all(., c("D_[0-9]__" = "", "D_1[0-9]__" = "", "Incertae Sedis"=NA_character_, "uncultured" = NA_character_, "Unassigned"=NA_character_)))) %>%
  mutate_at(vars(all_of(into2)), list(~ gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(all_of(into2)),  ~na_if(., ""))

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

rbind(smmy1, smmy2, smmy3) %>% 
  ggplot(aes(x = name, y = pct, color = DB, group = DB)) +
  facet_wrap( DB ~., scales = 'free_x') +
  geom_path(size = 1.5, alpha=0.6) +
  geom_point(size = 3, alpha=0.6) +
  geom_text(aes(label = paste0(round(pct*100, digits = 2), "%")), 
    size = 4, vjust = -1, family = "GillSans") +
  scale_y_continuous(labels = scales::comma) +
  theme_bw(base_size = 20, base_family = "GillSans") +
  theme(panel.border = element_blank(), legend.position = "top")


# BARPLOT =======


ab_f <-  read_tsv(list.files(path = wd, pattern = 'table_100_80', full.names = T), skip = 1)

MTD <-  read_tsv(list.files(path = wd, pattern = "mapping-file-corregido.tsv", full.names = T))


Freq <- function(x){rowSums(x > 1)}

ab_f %>%
  mutate(
    TotalAbundance = rowSums(across(where(is.double))),
    Prevalence = Freq(across(where(is.double)))) %>% 
  select(`Feature ID`, TotalAbundance, Prevalence) %>%
  arrange(desc(Prevalence)) -> ab

out <- tax_f %>%
  # select(-`Feature ID` , -Consensus) %>%
  separate(Taxon, sep = ";", into = into) %>% 
  mutate_at(vars(all_of(into)), list(~ str_replace_all(., c("[a-z]__" = "", "Unassigned"=NA_character_)))) %>%
  mutate_at(vars(all_of(into)), list(~ gsub("[[:space:]]", "", .))) %>%
  mutate_at(vars(all_of(into)),  ~na_if(., "")) %>% left_join(ab)

out %>% view()

out %>% ggplot(aes())

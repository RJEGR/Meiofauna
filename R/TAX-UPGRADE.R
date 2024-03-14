# datviz tax
# contrast old version from taxonomy.csv (122 samples) and this taxonmy file constructed from 189 sam.

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/Downloads/"

library(phyloseq)
library(microbiome)
library(tidyverse)

phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))

tax_f <- phyloseq %>% tax_table() %>% as(., "matrix") %>% as_tibble()

into <- c("k", "p", "c", "f", "g", "s")

# find contaminants?

contaminant_source <- c("Mammalia", "Aves")

contaminant_source <- paste(contaminant_source, collapse = "|")

contaminant_df <- tax_f %>% 
  filter_at(vars(all_of(into)), any_vars(grepl(contaminant_source, .)))
  # left_join(ab_f)

tax_f %>% select(all_of(into)) -> tax

# tax <- tax %>% mutate_all(function(x) {na_if(x,"")})

features <- colSums(!is.na(tax))


pct <- c(features / nrow(tax))

caption <- "feaures-from-189-samples"

df1 <- data.frame(into, features, pct, g = caption)

tax %>% drop_na(k)

df1 %>%
  mutate(into = factor(into, levels = into)) %>%
  ggplot(aes(x = into, y = features)) +
  geom_path(size = 1.5, alpha=0.6, group = 1) +
  geom_point(size = 3, alpha=0.6) +
  geom_text(aes(label = paste0(round(pct*100, digits = 2), "%")), 
    size = 4, vjust = -1, family = "GillSans") +
  scale_y_continuous(labels = scales::comma) +
  labs(y = "Number of sequence features", x = '', caption = caption) +
  theme_bw(base_size = 14, base_family = "GillSans") -> ps 

ps + theme(panel.border = element_blank()) -> ps
ps

tax %>%  drop_na(g) %>% view()

# 2)

nt <- function(x) {length(na.omit(unique(x)))}

tax %>% distinct(c) %>% drop_na()

# k  p  c  f  g  s 
# 5 47 48 52 70  3 

ntax <- apply(tax, 2, nt)

data.frame(x = names(ntax), ntax) %>%
  mutate(x = factor(x, levels = into)) %>%
  ggplot(aes(x = x, y = ntax)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=ntax), size = 3) +
  geom_text(aes(label = ntax),vjust = -0.5, hjust = 0.5 ) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank()) +
  labs(x = "", subtitle = 'Features with taxonomic assignation', y = "Number of taxons")

# # 3) Remove contaminants source: ----
# as in script DB_EXPLORATORY COMPARISON

# In addition, add prevalence and abundance

Freq <- function(x){rowSums(x > 1)}

phyloseq %>% otu_table() %>% as(., "matrix") %>% as_tibble(rownames = "Feature.ID") %>%
  mutate(
    TotalAbundance = rowSums(across(where(is.integer))),
    Prevalence = Freq(across(where(is.integer)))) %>% 
  select(`Feature.ID`, TotalAbundance, Prevalence) %>%
  arrange(desc(Prevalence)) -> ab_f

tReads <- sum(ab_f$TotalAbundance)

ab_f <- ab_f %>% mutate(pct_ab = TotalAbundance/tReads)

ab_f %>% 
  left_join(tax_f) %>%
  ggplot(aes(x = log10(TotalAbundance), y = Prevalence)) + 
  geom_point() + facet_wrap(~ k, scales = "free")
  

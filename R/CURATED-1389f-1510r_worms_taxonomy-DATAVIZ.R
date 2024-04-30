#
# TEST DATABASE COMPOSITION
# 

library(tidyverse)


wd <- "~/Documents/MEIOFAUNA_PAPER/DB_COMPARISON/"


tax_f <- read_tsv(list.files(path = wd, pattern = 'CURATED-1389f-1510r_worms_taxonomy.tsv$', full.names = T))


into <- c("k", "p", "c", "f", "g", "s")

# Res <- function(x) { rowSums(!is.na(x)) }
# 
# tax_f <- tax_f %>%
#   separate(Taxon, sep = ";", into = into)
# 
# 
# # RESOLUTION PER TAXON ----
# # CUANTOS TAXONES DIFERENTES POR RANK
# 
# tax_f %>% select(all_of(into)) -> tax
# 

into <- c("k", "p", "c", "o","f", "g", "s")

mutate_ranks <- c("none" = NA_character_, "os__" = "", "[a-z]__" = "", 
  "Incertae_Sedis" = NA_character_, 
  "uncultured" = NA_character_, 
  "Unassigned"=NA_character_,
  "Unknown" = NA_character_)

taxdb <- tax_f %>% 
  separate_wider_delim(cols = Taxon, delim = ";", names = into, too_few = "align_start") %>%
  # separate(Taxon, sep = ";", into = into) %>%
  select(all_of(into)) %>%
  mutate_all(list(~ str_replace_all(., mutate_ranks))) %>%
  mutate_all(function(x) {na_if(x,"")}) 



taxdb <- taxdb %>% mutate_all(function(x) {na_if(x,"")})

# max.rank <-  ncol(tax) 

# res <- colSums(!is.na(tax))

# tail(res <- rowSums(!is.na(tax))) # max.rank - 
# features <- nrow(tax) - colSums(!is.na(tax))

features <- colSums(!is.na(taxdb))


pct <- c(features / nrow(taxdb))

caption <- "CURATED-1389f-1510r_worms_taxonomy"

df1 <- data.frame(into, features, pct, g = caption)

# tax %>% drop_na(p)

df1 %>%
  mutate(into = factor(into, levels = into)) %>%
  ggplot(aes(x = into, y = features)) +
  geom_path(size = 1.5, alpha=0.6, group = 1) +
  geom_point(size = 3, alpha=0.6) +
  geom_text(aes(label = paste0(round(pct*100, digits = 2), "%")), 
            size = 4, vjust = -1, family = "GillSans") +
  scale_y_continuous(labels = scales::comma) +
  labs(y = "Number of sequences", x = '', caption = caption) +
  theme_bw(base_size = 14, base_family = "GillSans") -> ps 

ps + theme(panel.border = element_blank()) -> ps

# ps

# how many taxa

nt <- function(x) {length(na.omit(unique(x)))}

taxdb %>% distinct(c) %>% drop_na()

# k    p    c    f    g    s 
# 5  110  431  782 1306 1751 

ntax <- apply(taxdb, 2, nt)

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
  labs(x = "", subtitle = 'Features with taxonomic assignation', y = "Number of taxons") -> p2

library(patchwork)

p <- p2 / plot_spacer() / ps + plot_layout(heights = c(0.05, -0.01, 0.05))

ggsave(p, filename = 'CURATED-1389f-1510r_worms_taxonomy.png', path = path, width = 5, height = 10, device = png)


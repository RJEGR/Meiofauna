

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/Downloads/"

library(phyloseq)
library(microbiome)
library(tidyverse)

phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))

# set the contrast (look at the contrasts)
sample_data(phyloseq) %>% with(., table(Region))

ancombc_data1 <-  phyloseq %>%
  # subset_samples(Tissue=="Foregut") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) # %>%
# aggregate_taxa(., "p")


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

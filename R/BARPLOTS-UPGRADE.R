

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/"

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/"


library(phyloseq)
library(microbiome)
library(tidyverse)

# phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))
phyloseq <- read_rds(paste0(wd, '/ps.rds'))


# BARPLOT =======
colnames(tax_table(phyloseq))
colnames(sample_data(phyloseq))

# phyloseq %>%
#   aggregate_taxa(., "k") %>%
#   prune_samples(sample_sums(.) > 0, .) %>%
#   prune_taxa(taxa_sums(.) > 0, .) %>%
#   transform_sample_counts(function(x) x / sum(x)) %>%
#   plot_bar(x="LIBRARY_ID", y="Abundance", fill= "k")

# ps %>% microbiome::psmelt2(ps)

# ps_to_df <- function(ps) {
# as_tibble(as(otu_table(ps), "matrix"), rownames = "rowid")
# as_tibble(as(tax_table(ps), "matrix"), rownames = "rowid")
# as_tibble(as(tax_table(ps), "matrix"), rownames = "rowid")
#   
# }
  

reg_levels <- c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")

.df <- phyloseq %>%
  # aggregate_taxa(., "k") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  microbiome::psmelt2(sample.column = "SampleID", feature.column = "Feature") %>%
  filter(value > 0) %>%
  mutate(Region = factor(Region, levels = reg_levels))

.df <- .df %>%
  mutate(SampleID = KEYID.x) %>%
  filter(!grepl("CHAPO",SampleID ))

.df %>% filter(k %in% "Plantae") %>% group_by(p) %>% summarise(n = n(), N = sum(value))
  
.df <- .df %>% filter(!k %in% "Plantae")

filter_df <- function(df, agglom_lev = "k", low_ab = 0.01) {
  
  recode_to <- c(`ZZLow` = "Low abundance", `ZZNo hit` = "No hit")
    
  df %>%
    dplyr::rename( "Level" = all_of(agglom_lev)) %>%
    # drop_na(Level) %>%
    group_by(SampleID) %>%
    mutate(RA = value/sum(value)) %>%
    ungroup() %>%
    mutate(Level = ifelse(RA < low_ab, "ZZLow", Level)) %>%
    mutate(Level = ifelse(is.na(Level), "ZZNo hit", Level)) %>%
    mutate(Level = dplyr::recode(Level, !!!recode_to))
}

filter_df(.df, agglom_lev = "p") %>% group_by(SampleID) %>% summarise(sum(RA))

barTax <- function(df, agglom_lev = "k", low_ab = 0.01) {
  
  # df <- ps %>%
  #   aggregate_taxa(., agglom_lev) %>%
  #   prune_samples(sample_sums(.) > 0, .) %>%
  #   prune_taxa(taxa_sums(.) > 0, .) %>%
  #   microbiome::psmelt2()
  
  library(ggsci)
  
  df <- df %>% filter_df(agglom_lev, low_ab)
  
  
  labels <- df %>% filter(!Level %in% c("No hit", "Low abundance")) %>% pull(Level) %>% unique() %>% sort()
  
  # labels[which(labels %in% "Low abundance")] <- NA
  # labels[which(labels %in% "No hit")] <- NA
  # labels <- labels[!is.na(labels)]
  
  colourCount <- length(labels)
  
  if(colourCount > 7) {
    library(ggsci)
    getPalette <- colorRampPalette(pal_igv()(7))(colourCount)
  } else {
    getPalette <- pal_igv("default")(colourCount)
  }
  
  # pal_locuszoom
  # pal_ucscgb
  
  getPalette <- c(getPalette, "gray90", "black")
  
  labels <- c(labels, c("No hit", "Low abundance"))
  
  scale_fill_val <- structure(getPalette, names = labels)
  # 
  # getPalette[length(getPalette)] <- "Black"
  # labels[length(labels)] <- "Low abundance"
  # 
  # getPalette[length(getPalette)-1] <- "grey90"
  # labels[length(labels)-1] <- "No hit"
  
  df %>%
    mutate(Level = factor(Level, levels = labels)) %>%
    ggplot(aes(y = RA, x = SampleID, fill = Level)) +
    facet_grid( ~ Region, scales = "free", space = "free") +
    geom_col() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_fill_manual(agglom_lev, values = scale_fill_val) +
    labs(y = "Relative abundance (%)", x = "Samples") +
    theme_classic(base_size = 16, base_family = "GillSans") +
    theme(
      # axis.ticks.x = element_blank(),
      # axis.text.x = element_blank(),
      axis.ticks.length.x = unit(-1.5, "mm"),
      # width of tick marks in mm
      axis.ticks = element_line(size = .3),
      axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = .5, size = 4,
        margin = unit(c(t = 0.5, r = 0, b = 0, l = 0), "mm")),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
          strip.background = element_rect(fill = 'grey', color = 'white')) +
    guides(fill = guide_legend("")) 
  
}

# ggsavepath <- paste0(wd, '/Figures/')

# filtering low confidence assignments

# conf_val <- quantile(as.numeric(unique(.df$Confidence)), probs = seq(0, 1, 0.05))
#   
# conf_val <- conf_val[names(conf_val) %in% "75%"]
# 
# .df %>% mutate(Confidence = as.numeric(Confidence)) %>%
#   dplyr::rename( "Level" = all_of("k"))
#   mutate(Level = ifelse(Confidence >= conf_val, Level, "No hit"))
# 
  
ps <- barTax(.df, agglom_lev = "k", low_ab = 0) + 
  theme(legend.position = "top") +
  guides(fill = guide_legend(title = "", nrow = 1 ))

ggsave(ps, path = wd, filename = 'barplot_k_2.png', 
  device = png, 
       width = 12, height = 6)

barTax(.df, agglom_lev = "p", low_ab = 0.1)

# dvz <- .df %>% mutate(p = ifelse(as.numeric(Confidence) >= 0.5, p, "No hit"))


# TOP PHYLA 

out <- filter_df(.df) %>% filter(Level != "Low abundance")

MTD <- sample_data(phyloseq)

# color_vector <- as.character(unique(MTD$Region))

"#082F6B"

getPalette <- c("#000056", "#2E71A7","#60A4CF", "#9ECAE1", "#C8E0EF", "#E6F0F9")

axis_col <- structure(getPalette, names = reg_levels)

agglom_lev <- "p"

out <- out %>%
  dplyr::rename( ".Level" = "Level") %>%
  dplyr::rename( "Level" = agglom_lev) %>%
  mutate(Level = ifelse(as.numeric(Confidence) >= 0.5, Level, "No hit")) %>%
  drop_na(Level) %>%
  mutate(Level = ifelse(.Level == "Low abundance","Low abundance", Level))
  # mutate(Level = ifelse(RA < low_ab, "ZLow", Level))


ps <- .df %>%
  dplyr::rename( "Level" = agglom_lev) %>%
  # mutate(Level = ifelse(as.numeric(Confidence) >= 0.5, Level, NA)) %>%
  drop_na(Level) %>%
  # group_by(SampleID) %>%
  group_by(Region, Level) %>%
  summarise(TotalAbundance = sum(value)) %>%
  group_by(Level) %>%
  mutate(RA = TotalAbundance/sum(TotalAbundance)) %>%
  ggplot(aes(x = RA, y = Level, fill = Region)) +
  # facet_grid(Clade~., scales = "free", space = "free") +
  geom_col(position = position_stack(reverse = T)) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(limits = rev) +
  labs(x = "Relative abundance (%)", y = "Phylum") +
  scale_fill_manual(values = axis_col ) +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey', color = 'white')) +
  guides(fill = guide_legend(""))

ggsave(ps, path = wd, filename = 'barplot_p.png', 
  device = png, 
  width = 6, height = 7.5)


# NESTED BAR
# devtools::install_github("gmteunisse/ggnested")
# EX. https://github.com/gmteunisse/fantaxtic?tab=readme-ov-file
library(ggnested)


ps <- filter_df(.df, agglom_lev = "p", low_ab = 0) %>%
  ggnested(
  aes(x = SampleID,
    main_group = k,
    sub_group = Level)) +
  geom_bar(aes(weight = RA)) +
  facet_grid(~ Region, scales = "free", space = "free") +
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme( axis.ticks.x = element_blank(),
    axis.text.x = element_blank())

ggsave(ps, path = wd, filename = 'ggnested.png', 
  device = png, 
  width = 12, height = 6)

pal <- nested_palette(.df %>% drop_na(p), group = "k", subgroup = "p")

pal %>% distinct(k, group_colour) %>% pull(k, name = group_colour)
pal %>% distinct(p, subgroup_colour) %>% pull(p, name = subgroup_colour)


devtools::install_github("gmteunisse/fantaxtic")
require("fantaxtic")

data(GlobalPatterns)

top_asv <- top_taxa(GlobalPatterns, n_taxa = 10)
plot_nested_bar(ps_obj = top_asv$ps_obj,
  top_level = "k",
  nested_level = "p")


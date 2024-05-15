
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/"

library(tidyverse)

phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))

ab_f <- phyloseq %>%
  # transform_sample_counts(., function(x) x / sum(x)) %>%
  otu_table() %>% as("matrix") %>% as_tibble(rownames = "Feature.ID")

vars_to_numeric <- c("Depth","Latitude",	"Longitude", 
  "Clay",	"Silt",	"Sand",	"IC",	"TOC",	"CN",	"Oxygen", 
  "Finas",	"Medias",	"Muy_finas",	"Gruesas")

datTraits <- phyloseq %>% sample_data() %>% as(., "matrix")

as.zero <- function(x) { ifelse(is.na(x), 0, x)}

datTraits <- datTraits %>% 
  as_tibble() %>% 
  mutate_at(c(vars_to_numeric), as.numeric) %>%
  mutate_at(c(vars_to_numeric), as.zero) %>%
  # data.frame(row.names = "Feature.ID") %>%
  select(all_of(c("LIBRARY_ID", vars_to_numeric, "Region"))) %>%
  mutate(Region = factor(Region, levels = c("Yucatan", "NW Shelf", "NW Slope", "Deep-sea")))



tax_f <- phyloseq %>% tax_table() %>% as(., "matrix") %>% as_tibble()

# 1) Input Taxa-groups list ----


UNI_EUK <- c("Amoebozoa",
  "Archaeplastida",
  "Excavata",
  "Fungi",
  "Stramenopiles",
  "Alveolata",
  "Rhizaria")


METAZOAN <- c(
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



METAZOAN <- str_to_lower(paste(METAZOAN, collapse = "|"))
UNI_EUK <- str_to_lower(paste(UNI_EUK, collapse = "|"))


into <- c("k", "p", "c", "o", "f", "g", "s")

METAZOAN_DF <- tax_f %>%
  mutate_at(vars(all_of(into)), list(~ str_to_lower(.))) %>%
  filter_at(vars(all_of(into)), any_vars(grepl(METAZOAN, .))) %>%
  mutate_at(vars(all_of(into)), list(~ str_to_sentence(.))) %>%
  mutate(Clade = "Metazoans")


UNI_EUK_DF <- tax_f %>%
  mutate_at(vars(all_of(into)), list(~ str_to_lower(.))) %>%
  filter_at(vars(all_of(into)), any_vars(grepl(UNI_EUK, .))) %>%
  mutate_at(vars(all_of(into)), list(~ str_to_sentence(.))) %>%
  mutate(Clade = "Unicellular eukaryotes")

.tax_f <- select(rbind(METAZOAN_DF,UNI_EUK_DF), Feature.ID, Clade)

tax_f <- tax_f %>% left_join(.tax_f, by = "Feature.ID")

# Color scale ----

reg_levels <- c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")
getPalette <- c("#000056", "#2E71A7","#60A4CF", "#9ECAE1")
axis_col <- structure(getPalette, names = reg_levels)

# axis_col[names(axis_col) %in% "Deep-sea"] <- "#313695"


as_tibble(axis_col, rownames = "Region") %>%
  rename("axis_col" = "value") %>%
  right_join(datTraits, by = "Region", multiple = "all") %>% 
  pull(axis_col, name = LIBRARY_ID) -> true_species_cols


# mtd %>% view()

# NMDS =====

# https://www.davidzeleny.net/anadat-r/doku.php/en:pca_r

agglom_lev <- "p"

which_sam <- colnames(ab_f)

# 1) Agglomerate data (best solution)

data <- ab_f %>% 
  right_join(tax_f, by = "Feature.ID") %>%
  # filter(Clade %in% "Metazoans") %>%
  rename( "Level" = agglom_lev) %>%
  group_by(Level) %>% 
  # summarise_at(vars(all_of(which_sam)), sum) %>%
  summarise_if(is.integer, sum) %>%
  mutate(Level = ifelse(is.na(Level), "No hit", Level)) %>%
  # drop_na(Level) %>% 
  data.frame(., row.names = .$Level) %>%
  select(-Level)

data <- apply(data, 2, function(x) x / sum(x))

# data <- vegan::rarefy(data, min(rowSums(data)))

colSums(data)

# data <- vegan::decostand(data, method = "hellinger") # chord or hellinger

# data <- vegan::rda(data)

set.seed(202405)

mMDS <- vegan::metaMDS(t(data), distance = "bray", k = 2, trymax = 1000)

# cca_df <- vegan::cca(data, )


# vegan::stressplot(mMDS)
  
w_mtd <- c("LIBRARY_ID", "Region") # "LIBRARY_ID","Zone", "Description", "Profundidad", "Depth", 

vegan::scores(mMDS, tidy = T) %>% 
  as_tibble(rownames = "LIBRARY_ID") %>% 
  mutate(`LIBRARY_ID` = gsub("[.]","-", LIBRARY_ID)) %>% 
  filter(score %in% 'sites') %>%
  left_join(datTraits) -> MDSdf

vegan::scores(mMDS, tidy = T) %>% 
  as_tibble() %>% 
  filter(!score %in% 'sites') -> MDSdf_sp

ggplot() +
  geom_point(data = MDSdf, aes(x = NMDS1, y = NMDS2, color = Region), size = 5, alpha = 0.7) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  guides(color=guide_legend("",nrow=1)) +
  scale_color_manual(values = axis_col) +
  geom_text(data = MDSdf_sp, 
    aes(x = NMDS1, y = NMDS2, label = label), 
    family = "GillSans")

# PCA =====
# METAZOAN Euk (Meiofauna)

tax_f %>% count(Clade)

# #SampleID to LIBRARY_ID
# Feature ID to Feature.ID

PCA_out <- function(ab_f, w_clade = ...) {
  
  data <- ab_f %>%
    right_join(tax_f, by = "Feature.ID") %>%
    filter(Clade %in% w_clade) %>%
    # rename( "Level" = agglom_lev) %>%
    # drop_na(Level) %>%
    mutate_at(vars(all_of(which_sam)), function(x) {1E6 * x/sum(x)}) %>%
    select_at(vars(all_of(which_sam), `Feature.ID`)) %>%
    data.frame(., row.names = .$`Feature.ID`) %>%
    select(-`Feature.ID`)
  

  PCA = prcomp(t(log2(data+1)), center = FALSE, scale. = FALSE)
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])
  
  PCAdf %>%
    mutate(`LIBRARY_ID` = rownames(.)) %>%
    mutate(`LIBRARY_ID` = gsub("[.]","-", `LIBRARY_ID`)) %>%
    left_join(datTraits, by = "LIBRARY_ID") %>%
    mutate(Clade = w_clade) -> PCAdf
  
  
  PCAdf %>%
    ggplot(., aes(PC1, PC2)) +
    facet_grid(~ Clade) +
    geom_point(aes(color = Region), size = 5, alpha = 0.7) +
    labs(caption = '') +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme_bw(base_family = "GillSans", base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
    coord_fixed(ratio = sd_ratio) +
    guides(color=guide_legend("",nrow=1)) +
    scale_color_manual(values = axis_col)
  
  
  return(PCAdf)
  
}

p1 <- PCA_out(ab_f, w_clade = "Unicellular eukaryotes")
p2 <- PCA_out(ab_f, w_clade = "Metazoans")


library(patchwork)
 
psave <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = 'top')

psave

ggsavepath <- paste0(wd, '/Figures/')

ggsave(psave, path = ggsavepath, filename = 'PCA.png', device = png, width = 10, height = 10)


# as distance tree -----


dist.method <- 'euclidean'
linkage.method <- 'complete'

dat <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

hc <- hclust(dist(dat, method = dist.method), 
  method = linkage.method)

dend <- hc %>% as.dendrogram()


# Find special clusters:
library(dynamicTreeCut)

clusters <- cutreeDynamic(hc, distM = as.matrix(dist(dat)), method = "tree")
# we need to sort them to the order of the dendrogram:
clusters <- clusters[order.dendrogram(dend)]

clusters_numbers <- unique(clusters) - (0 %in% clusters)
n_clusters <- length(clusters_numbers)


# devtools::install_github('talgalili/dendextend')

library(dendextend)

dend2 <- dend %>%
  hang.dendrogram(hang = 0.01) %>%
  ladderize() %>%
  set("labels_cex", 0.7) %>%
  # set("nodes_pch", c(NA,NA,NA,19)) %>%
  # set("nodes_col", 3) %>%
  # set("leaves_col", 2) %>%
  # branches_attr_by_clusters(clusters, values = axis_col) %>%
  color_labels(col =  true_species_cols)

dend2 %>% plot()

legend("topleft", legend = names(axis_col), fill = axis_col)

# colored_bars(clusters, dend, sort_by_labels_order = FALSE)

# BIND AB, MTD and tax
# Using Region as discriminatory


# candisc canonical discriminant analysis (CDA) -----
library(candisc)

w_cols <- ab_f %>% select_if(is.integer) %>% names()

agglom_lev <- "p" 
  
ab_f %>%
  left_join(tax_f, by = "Feature.ID") %>% 
  # filter(Clade %in% "Metazoans") %>%
  pivot_longer(cols = all_of(w_cols), names_to = "LIBRARY_ID") %>%
  filter(value > 0) %>%
  left_join(datTraits) %>%
  rename( "Level" = agglom_lev) %>%
  group_by(Region, Level) %>% 
  summarise(value = sum(value)) %>%
  mutate(Level = ifelse(is.na(Level), "No hit", Level)) %>%
  # mutate(value = log2(value + 1)) %>%
  mutate(value = as.numeric(value)) %>%
  # pivot_wider(names_from = p, values_from = value, values_fill = 0) %>%
  pivot_wider(names_from = Region, values_from = value, values_fill = 0) %>%
  ungroup() -> dat

dat <- dat %>%
  data.frame(., row.names = .$Level) %>%
  select(-Level)

dat <- apply(dat, 2, function(x) x / sum(x))

# dat <- round(dat)

colSums(dat)

dat <- dat %>% data.frame(p = rownames(.), .)

mod <- lm(cbind(`Yucatan`,`NW.Shelf`,`NW.Slope`, `Deep.sea`) ~  p, data=dat)

manova <- Anova(mod, type = as.character(2))

can <- candisc(mod, data = dat)


manova <- Anova(mod, type = as.character(2))


is.na(dat)

dat %>% view()

library(candisc)

# Ex

iris.mod <- lm(cbind(Petal.Length, Sepal.Length, Petal.Width, Sepal.Width) ~ Species, data=iris)
iris.can <- candisc(iris.mod, data=iris)
# plot(iris.can, col=col, pch=pch)
heplot(iris.can)


depVars <- dat %>% select_if(is.double) %>% names()

indepVars <- "Region"

form <- formula(paste('cbind(',
  paste(depVars, collapse = ','),
  ') ~ ',
  paste(indepVars, collapse = '+')))

mod <- lm(form, data=dat, na.action=na.exclude)

mod <- lm("Loukozoa ~ Region", data=dat, na.action=na.exclude)

# mod <- lm(cbind(Yucatan, `NW Shelf`, `NW Slope`,`Deep-sea`) ~  p, data=dat)

car::Anova(mod, test="Wilks")

# ERROR HERE ======

manova <- Anova(mod, type = as.character(2))

can <- candisc(mod, data = dat, term = 'Region', ndim = "2")


can <- candisc(mod, data = dat, term = 'p', ndim = "2")

mod <- lm(cbind(`Annelida`,`Arthropoda`,`Bryozoa`, `Cnidaria`) ~  Region, data=dat)

can <- candisc(mod, data = dat)



# heplot(candisc(mod, data = data, term = 'rank3', ndim = "1"))

# plot(can, col = axis_col)

ellipse_df <- can$means %>%
  data.frame() %>%
  rownames_to_column(var= "Region")

data.frame(can$structure+can$coeffs.std ) %>%
  rownames_to_column(var= 'Metric') -> segment_df

can$scores -> scores_df

x_label <- paste0("Can1 (", round(can$pct[1], 1), "%)")
y_label <- paste0("Can2 (", round(can$pct[2], 1), "%)")

ggplot() + 
  geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
  geom_point(data = scores_df, 
    aes(x = Can1,y = Can2, color = Region),  alpha = 0.5) +
  scale_color_manual("", values = axis_col) +
  scale_fill_manual("", values = axis_col) +
  ggforce::geom_mark_ellipse(
    data = ellipse_df, aes(Can1,y = Can2, color = Region, fill = Region)) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  geom_segment(data = segment_df, 
    aes(x = 0, y = 0, xend = Can1, yend = Can2),
    arrow = arrow(length = unit(0.07, "cm"))) +
  geom_text(data = segment_df, 
    aes(x = Can1+sign(Can1)*0.2, 
      y = Can2+sign(Can1)*0.2, label = Metric), 
    family = "GillSans") +
  coord_fixed() +
  labs(x = x_label, y = y_label) +
  theme(panel.border = element_blank(), legend.position = 'top')  -> p4

# p4 + facet_grid(~ Region) + theme(
#   strip.background = element_rect(fill = 'white', color = 'white')) +
#   xlim(c(-3, 3)) +
#   ylim(c(-2.5, 2.5)) -> p4


ggsave(p4, path = ggsavepath, 
  filename = 'candisc_108hpf.png',width = 4,height = 4, dpi = 300)



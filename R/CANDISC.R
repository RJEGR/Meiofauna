
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/Documents/MEIOFAUNA_PAPER/INPUTS/"

library(tidyverse)

ab_f <-  read_tsv(list.files(path = wd, pattern = 'table_100_80', full.names = T), skip = 1)

mtd <- read_tsv(list.files(path = wd, pattern = 'mapping-file-corregido.tsv', full.names = T))

w_mtd <- c("#SampleID","Zone", "Description", "Profundidad", "Depth", "Region")

# Color scale ----


color_vector <- as.character(unique(mtd$Region))
n <- length(color_vector)
# getPalette <- RColorBrewer::brewer.pal(n, 'Set1')
getPalette <- c("#4DAF4A", "#313695", "lightblue", "#E41A1C")

axis_col <- structure(getPalette, names = color_vector)

# axis_col[names(axis_col) %in% "Deep-sea"] <- "#313695"


as_tibble(axis_col, rownames = "Region") %>%
  rename("axis_col" = "value") %>%
  right_join(mtd, by = "Region") %>% 
  pull(axis_col, name = `#SampleID`) -> true_species_cols


# mtd %>% view()

# PCA -----
# https://www.davidzeleny.net/anadat-r/doku.php/en:pca_r

data <- ab_f %>% select_if(is.double) 

# data <- vegan::decostand(data, method = "hellinger") # chord or hellinger

PCA = prcomp(t(log10(data+1)), center = FALSE, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

PCAdf %>%
  mutate(`#SampleID` = rownames(.)) %>%
  left_join(mtd %>% select(all_of(w_mtd))) -> PCAdf


PCAdf %>%
  ggplot(., aes(PC1, PC2)) +
  geom_point(aes(color = Region), size = 5, alpha = 0.7) +
  labs(caption = '') +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'top') +
  # coord_fixed(ratio = sd_ratio) +
  guides(color=guide_legend("",nrow=1)) +
  scale_color_manual(values = axis_col) 
  # see::scale_color_colorhex() 

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


w_cols <- ab_f %>% select_if(is.double) %>% names()

ab_f %>%
  left_join(tax_f %>% select(`Feature ID`, rank3), by = "Feature ID") %>% 
  pivot_longer(cols = all_of(w_cols), names_to = "#SampleID") %>%
  filter(value > 0) %>%
  left_join(mtd %>% select(all_of(w_mtd))) %>%
  group_by(Region, Depth, rank3) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = rank3, values_from = value, values_fill = 0) %>%
  # pivot_wider(names_from = Region, values_from = value, values_fill = 0) %>%
  # pivot_wider(names_from = Depth, values_from = value, values_fill = 0) %>%
  ungroup() -> dat

# top abundance ranks

ab_f %>%
    left_join(tax_f %>% select(`Feature ID`, rank3), by = "Feature ID") %>% 
    pivot_longer(cols = all_of(w_cols), names_to = "#SampleID") %>%
    filter(value > 0) %>%
    left_join(mtd %>% select(all_of(w_mtd))) %>%
    group_by(rank3) %>% summarise(value = sum(value))  %>% arrange(desc(value)) %>%
    pull(rank3) %>% head(8)

library(candisc)

# data %>% filter(hpf == 108) -> dat

wnames <- dat %>% select_if(is.double) %>% names()

# mod <- lm(cbind(`0-100`,`100-200`,`200-500`,`500-1000`,`1000-1500`,`1500-2000`,`2000-2500`, `2500-3000`, `3000-3500`, `3500-4000`) ~ Region , data=dat)

# mod <- lm(cbind(`Deep-sea`,`NW Shelf`,`NW Slope`,Yucatan) ~  Depth, data=dat)

mod <- lm(cbind(`Nucletmycea`,`Holozoa`,`Chloroplastida`,`Alveolata`, `Rhizaria`, `Discoba`, `Stramenopiles`) ~  Region, data=dat)

can <- candisc(mod, data = dat, term = 'Region', ndim = "2")

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



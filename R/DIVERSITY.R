# GENERATE ALPHA DIVERSITY


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/"


reg_levels <- c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")
getPalette <- c("#000056", "#2E71A7","#60A4CF", "#9ECAE1")
axis_col <- structure(getPalette, names = reg_levels)

library(phyloseq)
library(tidyverse)

phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))

# rarefy firts

x<- as(otu_table(phyloseq), "matrix")

rs <- rowSums(x)

quantile(rs)

library(vegan)
S <- specnumber(x)
Srar <- rarefy(x, min(rs))
plot(S, Srar, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

r <- rarecurve(t(x), step=50, cex=0.5, tidy = T)

ps <- r %>% 
  full_join(sample_data(phyloseq), by = c("Site" = "LIBRARY_ID")) %>%
  as_tibble() %>%
  mutate(Region = factor(Region, levels = reg_levels)) %>%
  ggplot(aes(x = Sample, y = Species, group = Site)) + 
  geom_path(aes(color = Region)) + # 
  labs(x = "Reads", y = "Number of ASVs") +
  scale_x_continuous(labels = scales::comma) +
  # geom_vline(xintercept =   median(sample_sums(phyloseq)), linetype = "dashed") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  scale_color_manual(values = axis_col) +
  facet_wrap(~ Region, nrow = 1, scales = "free") +
  theme(legend.position='bottom', legend.justification = "right")

ggsave(ps, path = wd, filename = 'rarefied-asvs.png', width = 10, height = 2.5, device = png, dpi = 300)

# rarefy without replacement
ps.rarefied = rarefy_even_depth(phyloseq, 
  rngseed=1, sample.size=min(sample_sums(phyloseq)), replace=F, trimOTUs = F)

# Alpha & Beta ----

measures <- c("Observed", "Chao1", "Shannon", "InvSimpson")

R <- estimate_richness(ps.rarefied, measures = measures)

# calculate evenness index using vegan package

Evenness <- R$Shannon / log(R$Observed)

# Evenness <- diversity(x, index = "shannon") / log(specnumber(x)) 

head(R <- cbind(R, Evenness))

rownames(R) <- gsub("[.]", "-", rownames(R))

MTD <- sample_data(phyloseq)

MTD <- MTD[match(rownames(MTD), rownames(R)),]

R <- data.frame(MTD,R)

sample_data(phyloseq) <- R 


# saveRDS(phyloseq, paste0(wd, "/phyloseq.rds"))


R %>% 
  pivot_longer(cols = c(measures, "Evenness")) %>%
  mutate(Region = factor(Region, levels = reg_levels)) %>%
  mutate(name = factor(name, levels = c(measures, "Evenness")))%>% 
  mutate(x = as.numeric(`X..TOC`)) %>%
  ggplot(aes(x = x, y = value)) +
  facet_wrap(~ name, scales = "free_y", nrow = 1) +
  geom_point()

# plots

reg_levels <- c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")


p <- R %>% 
  pivot_longer(cols = c(measures, "Evenness")) %>%
  mutate(Region = factor(Region, levels = reg_levels)) %>%
  mutate(name = factor(name, levels = c(measures, "Evenness"))) %>%
  ggplot(aes(x = Region, y = value)) +
  facet_wrap(~ name, scales = "free_y", nrow = 1) + 
  stat_boxplot(geom ='errorbar', width = 0.15,
    position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6), outlier.shape=NA) +
  stat_summary(fun = mean, geom="point", shape=20, 
    size = 3, color="red", fill="red") +
  coord_cartesian(ylim=c(0,NA)) +
  labs(y = "Index") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position='bottom', legend.justification = "right",
    panel.grid = element_blank(), legend.text = element_text(size = 7),
    strip.background = element_rect(fill = 'grey88'))

p

ggsave(p, path = wd, filename = 'diversity-indexes.png', width = 12, height = 2.5, device = png, dpi = 300)

# but also
library(microbiome)
tab <- alpha(phyloseq, index = "all")


# rarefacion =====
# TOMA MUCHO TIEMPO

# iNEXT computes the following two types of rarefaction (interpolation) and extrapolation (prediction) and the associated 95% confidence intervals:
  
# 1. Sample-size-based rarefaction and extrapolation: diversity estimates for rarefied and extrapolated samples up to a maximum size (double reference sample size by default or a user-specified endpoint); see below.
# 2. Coverage-based rarefaction and extrapolation: diversity estimates for rarefied and extrapolated samples for sample coverage up to a maximum coverage that is obtained from the double reference sample size by default or a user-specified endpoint; see below).

library(iNEXT)

# test <- t(ab)[1,]
# i.zero <- which(test == 0)
# test.no.zero <- test[-i.zero]

ab <- as(otu_table(phyloseq), "matrix")

x <- apply(ab, 2, function(x) x[-which(x == 0)])

out <- iNEXT(x, q=c(0), datatype="abundance")


# Sample-size-based R/E curves, separating by "site""
ggiNEXT(out, type=1, facet.var="none", grey = T)
ggiNEXT(out, type=2, facet.var="site", grey = T)
ggiNEXT(out, type=2, facet.var="none", color.var="site")



# test specpool
# 


x <- as(t(otu_table(phyloseq)), "matrix")

pool <- with(MTD, specpool(x, Region))
pool

estimateR(x[1:5])

op <- par(mfrow=c(1,2))

pool <- with(MTD, specpool(x, Region))



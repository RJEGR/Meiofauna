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

# Rarefaction curves ----

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

# rarefy without replacement (not recommended, McMurdie and Holmes (2014) Waste Not, Want Not ...)

# ps.rarefied = rarefy_even_depth(phyloseq, rngseed=1, sample.size=min(sample_sums(phyloseq)), replace=F, trimOTUs = F)

# # barplot(sample_sums(ps.rarefied))
# barplot(sample_sums(phyloseq))


# Deal w/ spurious singletons -----

tax <- as(tax_table(phyloseq), "matrix") %>% as_tibble(rownames = "Feature ID")

Freq <- function(x){rowSums(x > 1)}

prevelancedf <- as(otu_table(phyloseq), "matrix") %>%
  as_tibble(rownames = "Feature ID") %>%
  mutate(
    TotalAbundance = rowSums(across(where(is.integer))),
    Prevalence = Freq(across(where(is.integer)))) %>%
  select(`Feature ID`, TotalAbundance, Prevalence) %>%
  arrange(desc(Prevalence)) %>%
  left_join(tax)

mutate_ranks <- c("NA_character_" = "No hit")

prevelancedf %>%
  # mutate_all(list(~ str_replace_all(., mutate_ranks)))
  mutate(k = ifelse(is.na(k), "No hit", k)) %>%
  group_by(k) %>% 
  cor_test(vars = "TotalAbundance", vars2 = c("Prevalence"))

n_sam <- ncol(as(otu_table(phyloseq), "matrix"))

prevelancedf %>% 
  mutate(Level = k) %>%
  # drop_na(Level) %>%
  mutate(Level = ifelse(is.na(Level), "No hit", Level)) %>%
  # mutate(x = as.double(Confidence), color = ifelse(x < 0.25, "Low conf." ," Conf.")) %>%
  # mutate(color = ifelse(TotalAbundance < 2, "Low ab." ," ab")) %>%
  ggplot(aes(x = TotalAbundance, y = Prevalence/n_sam)) + # color = color
  geom_point(size = 2, alpha = 0.7) + 
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2, color = "blue") +
  # geom_vline(xintercept = 1, alpha = 0.5, linetype = 2) +
  scale_x_log10(
                breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "ASV Size", y = "Prevalence (Frac. Samples)") +
  # scale_color_manual(values = c('blue', 'red')) +
  facet_wrap( ~ Level, nrow = 1) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  theme_bw(base_size = 17, base_family = "GillSans") +
  theme(panel.grid.minor = element_blank(),
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 0, hjust = 1, vjust = 1, size = 12)) -> p2

p2 + theme(panel.border = element_blank(), legend.position = "top")

sum(keep <- taxa_sums(phyloseq) > 1)

# prevelancedf %>% filter(TotalAbundance <= 2) %>%
#   ggplot(aes(as.double(Confidence))) + geom_histogram()


ps <- prune_taxa(keep, phyloseq)


# Alpha & Beta ----

measures <- c("Observed", "Chao1", "Shannon", "InvSimpson")

R <- estimate_richness(ps, measures = measures)

# calculate evenness index using vegan package

Evenness <- R$Shannon / log(R$Observed)

# Evenness <- diversity(x, index = "shannon") / log(specnumber(x)) 

head(R <- cbind(R, Evenness))

rownames(R) <- gsub("[.]", "-", rownames(R))

MTD <- sample_data(ps)

MTD <- MTD[match(rownames(MTD), rownames(R)),]

dim(R <- data.frame(MTD,R))
dim(MTD)

sample_data(ps) <- R 


saveRDS(ps, paste0(wd, "/ps.rds"))


# plots

R <- sample_data(read_rds(paste0(wd, "/ps.rds"))) %>% as(., "matrix")



vars_to_numeric <- c("Depth","Latitude",	"Longitude", 
  "Clay",	"Silt",	"Sand",	"IC",	"TOC",	"CN",	"Oxygen", 
  "Finas",	"Medias",	"Muy_finas",	"Gruesas")

measures <- c("Observed", "Chao1", "Shannon", "InvSimpson", "Evenness")

R <- R %>% as_tibble() %>% mutate_at(c(vars_to_numeric, measures), as.double) 


reg_levels <- c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")

# Statistical test----

# 1) test if gaussianity (parcial FALSE)

df <- R %>% 
  pivot_longer(cols = measures, names_to = "rich_name", values_to = "Index") %>%
  mutate(Region = factor(Region, levels = reg_levels)) %>%
  select(Region, rich_name, Index)

df %>% 
  group_by(Region, rich_name) %>% 
  rstatix::shapiro_test(Index) %>%
  mutate(gauss = ifelse(p > 0.05, TRUE, FALSE)) %>% 
  arrange(rich_name)
  
# visualize 


qqfun <- function(x) {
  x <- x[x > 0]
  qq <- qqnorm(x, plot.it = F) 
  qq %>% as_tibble()
}

df %>% 
  # group_by(Region, rich_name) %>%
  reframe(qqfun(Index)) %>% # instead of summarise, use reframe
  # group_by(Region, rich_name) %>%
  mutate(outlier = is_outlier(y)) %>%
  ggplot(aes(x, y)) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
    se = TRUE, na.rm = TRUE) +
  geom_point(aes(color = outlier), size = 2.5, alpha = 0.8) + # , label.y = 2.5
  # facet_wrap(~Region) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", p.accuracy = 0.001) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  labs(x = "Expected (qqplot)", y = "Observed Richness") +
  scale_color_manual(name = expression("Outlier-"~sigma),  values = c('black', 'red')) +
  theme(legend.position = "top")

# 2) Homocelasticidad (partial FALSE) ----
# Before doing parametric or not test, lets to analyze homogeneity of variance across experimental using Levene’s test:

df %>% 
  group_by(rich_name) %>% 
  levene_test(Index ~ as.factor(Region)) %>%
  mutate(hom_var = ifelse(p > 0.05, TRUE, FALSE))

# If the p-value or significance is above  0.05, we can assume homogeneous variance based on the available data. Ex. The p value of 0.144 is greater than 0.05 so the null hypothesis is maintained and there is no difference between the variances


# Statistical priori test----
# Using welch when not assuming equal variance

df %>% 
  group_by(rich_name) %>%
  rstatix::welch_anova_test(Index ~ Region) %>%
  # rstatix::anova_test(r_ind_adj ~ pH) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() -> prior.welch.stats

# view(prior.welch.stats)

# Kruskal-Wallis test by rank is a non-parametric alternative to one-way ANOVA test, which extends the two-samples Wilcoxon test in the situation where there are more than two groups. It’s recommended when the assumptions of one-way ANOVA test are not met. 
# kruskal-wallis p/ mas de dos muestras independientes

df %>% 
  group_by(rich_name) %>%
  rstatix::kruskal_test(Index ~ Region) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p") -> prior.kruskal.stats

# view(prior.kruskal.stats) 


# POSTERIORI TEST

df %>%
  group_by(rich_name) %>%
  pairwise_wilcox_test(Index ~ Region,  conf.level = 0.95) %>% # if not parametric
  # pairwise_t_test(Index ~ Region) %>% # if parametric
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() -> post.test

post.test <- post.test %>% filter(p.adj.signif != "ns")

post.test %>% add_xy_position(x = "Region", scales = "free_y", dodge = 0.6, step.increase = 0.5) -> stats

# post.test %>% filter(p.adj.signif != "ns") %>%
#   pivot_longer(cols = c("group1","group2"), values_to = "Region") %>%
#   select(Region, rich_name, p.adj.signif) %>%
#   right_join(df, by = c("Region", ""))


p <- df %>%
  filter(rich_name %in% unique(stats$rich_name)) %>%
  mutate(rich_name = factor(rich_name, levels = measures)) %>%
  ggplot(aes(x = Region, y = Index)) +
  facet_wrap(~ rich_name, scales = "free_y", nrow = 1) +
  # stat_boxplot(geom ='errorbar', width = 0.15,  position = position_dodge(0.6)) +
  # geom_boxplot(width = 0.6, position = position_dodge(0.6), outlier.shape=NA) +
  geom_jitter(width=0.1,alpha=0.2, height = 0.1, size = 0.7, shape = 1) +
  # stat_summary(fun = mean, geom="point", shape=20, 
  #   size = 3, color="red", fill="red") +
  coord_cartesian(ylim=c(0,NA)) +
  labs(y = "Alpha diversity", x = "") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position='bottom', legend.justification = "right",
    panel.grid = element_blank(), legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
      margin = unit(c(t = 0.5, r = 0, b = 0, l = 0), "mm")))

# stats$y.position <- stats$y.position+2

p <- p + 
  ggpubr::stat_pvalue_manual(stats, label = "p.adj.signif", family = "GillSans",
    remove.bracket = F, tip.length = 0.01,  hide.ns = T, size = 2.5) 

p <- p + theme(strip.background = element_rect(fill = 'grey88', color = 'white'),
  panel.border = element_blank(), legend.position = 'top') 


ggsave(p, path = wd, filename = 'diversity-indexes2.png', width = 5.5, height = 2, device = png, dpi = 300)

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



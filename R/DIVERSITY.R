# GENERATE ALPHA DIVERSITY


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/"


library(phyloseq)
library(tidyverse)

phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))


# Alpha & Beta ----

measures <- c("Observed", "Chao1", "Shannon", "InvSimpson")

R <- estimate_richness(phyloseq, measures = measures)

rownames(R) <- gsub("[.]", "-", rownames(R))

MTD <- sample_data(phyloseq)

MTD <- MTD[match(rownames(MTD), rownames(R)),]

R <- data.frame(MTD,R)

sample_data(phyloseq) <- R 


saveRDS(phyloseq, paste0(wd, "/phyloseq.rds"))


# plots

reg_levels <- c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")


p <- R %>% 
  pivot_longer(cols = measures) %>%
  mutate(Region = factor(Region, levels = reg_levels)) %>%
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


ggsave(p, path = wd, filename = 'diversity-indexes.png', width = 12, height = 2.5, device = png, dpi = 300)

# test specpool
# 


x <- as(t(otu_table(phyloseq)), "matrix")

pool <- with(MTD, specpool(x, Region))
pool

estimateR(x[1:5])

op <- par(mfrow=c(1,2))

pool <- with(MTD, specpool(x, Region))



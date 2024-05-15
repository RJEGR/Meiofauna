# GENERATE NETWORK ANALYSIS USING WGCNA
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/"

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/"

library(phyloseq)
library(microbiome)
library(WGCNA)
library(flashClust)
library(tidyverse)


ps <- read_rds(paste0(wd, '/ps.rds'))

rank_names(ps)

ab_f <- ps %>% otu_table() %>% 
  as("matrix") %>% as_tibble(rownames = "Feature.ID")
  
tax_f <- ps %>% tax_table() %>% as(., "matrix") %>% as_tibble()

which_sam <- colnames(ab_f)

# 1) Agglomerate data (best solution)

agglom_lev <- "c"


datExpr <- ab_f %>% 
  right_join(tax_f, by = "Feature.ID") %>%
  rename( "Level" = all_of(agglom_lev)) %>%
  mutate(Level = ifelse(as.double(Confidence) > 0.8, Level, NA)) %>%
  mutate(Level = ifelse(is.na(Level), "No hit", Level)) %>%
  group_by(Level) %>% 
  # summarise_at(vars(all_of(which_sam)), sum) %>%
  summarise_if(is.integer, sum) %>%
  # drop_na(Level) %>% 
  data.frame(., row.names = .$Level) %>%
  select(-Level)

dim(datExpr)

datExpr <- t(datExpr) # log2(count+1) # 

str(datExpr)

cat("\n:::::\n")

gsg = goodSamplesGenes(datExpr, verbose = 3)

gsg$allOK

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

# When you have a lot of rows in the matrix (ex. genes or taxons) use the following code

max_power <- 30

cat(powers <- c(c(1:10), seq(from = 10, to = max_power, by=1)))

powers = unique(powers)

allowWGCNAThreads()

cor_method =  "cor" # by default WGCNA::cor(method =  'pearson') is used, "bicor"

corOptionsList = list(use ='p') # maxPOutliers = 0.05, blocksize = 20000

sft = pickSoftThreshold(datExpr, 
  powerVector = powers, 
  corFnc = cor_method,
  corOptions = corOptionsList,
  verbose = 5, 
  networkType = "unsigned")



soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 2)

power_pct <- soft_values[which.max(soft_values)]


softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]

meanK <- sft$fitIndices[softPower, 5]

# hist(sft$fitIndices[,5])

softPower <- min(softPower)

cat("\nsoftPower value", softPower, '\n')


title1 = 'Scale Free Topology Model Fit, R^2'
title2 = 'Mean Connectivity'

caption = paste0("Lowest power for which the scale free topology index reaches the ", power_pct*100, " %")

sft$fitIndices$mean.k.

sft$fitIndices %>% 
  mutate(scale = -sign(slope)*SFT.R.sq) %>%
  select(Power, mean.k., scale) %>% pivot_longer(-Power) %>%
  mutate(name = ifelse(name %in% 'scale', title1, title2)) %>%
  ggplot(aes(y = Power, x = value)) +
  geom_text(aes(label = Power), size = 2) +
  geom_abline(slope = 0, intercept = softPower, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = min(meanK), linetype="dashed", alpha=0.5) +
  labs(y = 'Soft Threshold (power)', x = '', 
    caption = caption) +
  facet_grid(~name, scales = 'free_x', switch = "x") +
  # scale_x_continuous(position = "top") +
  theme_bw(base_family = "GillSans",base_size = 16) -> psave

psave


# 2) 
# Build a adjacency "correlation" matrix

enableWGCNAThreads()

# specify network type
# softPower <- 30 # if asvs level

adjacency <- adjacency(datExpr, 
  power = softPower, 
  type = "unsigned")


TOM <- TOMsimilarity(adjacency, TOMType = "unsigned") # specify network type


dissTOM = 1 - TOM

# heatmap(dissTOM) # if more than 100 rows, omit

# Generate Modules ----
# Generate a clustered gene tree

geneTree = flashClust(as.dist(dissTOM), method="average")

dev.off()

plot(geneTree, xlab="", sub="", 
  main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

# This sets the minimum number of rows to cluster into a module

minClusterSize <- 1

# Using the median or mean k
# minClusterSize <- abs(sft$fitIndices[softPower, 5])

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
  deepSplit = 2, pamRespectsDendro = FALSE,
  minClusterSize = minClusterSize)

length(table(dynamicMods))


dynamicColors = labels2colors(dynamicMods)
names(dynamicColors) <- colnames(datExpr)
MEList = moduleEigengenes(datExpr, 
  colors= dynamicColors,
  softPower = softPower)

MEs = MEList$eigengenes

MEDiss= 1 - cor(MEs)

METree = flashClust(as.dist(MEDiss), method= "average")

plot(METree, main = "Clustering of module eigengenes",
  xlab = "", sub = "")


# Now we will see if any of the modules should be merged. I chose a height cut of 0.30, corresponding to a similarity of 0.70 to merge:

MEDissThres = 0.3

abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors

mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, dynamicColors, c("Modules"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

caption <- c("Dynamic Tree Cut", "Merged dynamic", "\n(cutHeight: ", MEDissThres,")")

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
  caption, dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors <- mergedColors

colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1


MEs <- mergedMEs


# how modules where obtained:
nm <- table(moduleColors)


cat('Number of mudules obtained\n :', length(nm))

# w/ p level result in 22 modules

# write_rds(as_tibble(moduleColors, rownames = "Family") %>% 
#             dplyr::rename("WGCNA" = "value"), 
#           file = paste0(path, "WGCNA_MIRS.rds"))

# Plot TOM
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
TOMplot(plotTOM, dendro = geneTree, Colors = moduleColors)

#

# datTraits <- phyloseq %>% sample_data() %>% with(., table(LIBRARY_ID, Region))

vars_to_numeric <- c("Depth","Latitude",	"Longitude", 
  "Clay",	"Silt",	"Sand",	"IC",	"TOC",	"CN",	"Oxygen", 
  "Finas",	"Medias",	"Muy_finas",	"Gruesas")


# measures <- c("Observed", "Chao1", "Shannon", "InvSimpson", "Evenness")


datTraits <- ps %>% sample_data() %>% as(., "matrix")

as.zero <- function(x) { ifelse(is.na(x), 0, x)}

datTraits <- datTraits %>% as_tibble(rownames = "rows") %>% 
  mutate_at(c(vars_to_numeric), as.numeric) %>%
  mutate_at(c(vars_to_numeric), as.zero) %>%
  data.frame(row.names = "rows") %>%
  select(all_of(c(vars_to_numeric))) 

rownames(datTraits) <- gsub("-", ".",rownames(datTraits))

identical(rownames(datExpr), rownames(datTraits))

# Recalculate MEs with color labels

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

MEs = orderMEs(MEs0)

names(MEs) <- str_replace_all(names(MEs), '^ME', '')

# heatmap(as(MEs, "matrix"))

moduleTraitCor = WGCNA::cor(MEs, datTraits, use= "p")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datTraits))

moduleTraitCor %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'moduleTraitCor') -> df1

moduleTraitPvalue %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'corPvalueStudent') %>%
  right_join(df1) -> df1

hclust <- hclust(dist(moduleTraitCor), "complete")

hc_order <- hclust$labels[hclust$order]

df1 %>%
  mutate(star = ifelse(corPvalueStudent <.001, "***", 
    ifelse(corPvalueStudent <.01, "**",
      ifelse(corPvalueStudent <.05, "*", "")))) -> df1



df1 <- df1 %>% mutate(name = factor(name, levels = vars_to_numeric))

lo = floor(min(df1$moduleTraitCor))
up = ceiling(max(df1$moduleTraitCor))
mid = (lo + up)/2

# add module size

reads <- colSums(datExpr)
Total <- sum(reads)

# Sanity check:

identical(names(colSums(datExpr)), names(moduleColors))

data.frame(reads, moduleColors) %>% 
  as_tibble(rownames = "Name") %>% 
  group_by(moduleColors) %>% 
  summarise(n = n(), reads = sum(reads)) %>%
  dplyr::rename('module' = 'moduleColors') -> stats

module_size <- structure(stats$n, names = stats$module)

module_size <- module_size[match(hc_order, names(module_size))]

identical(hc_order, names(module_size))

hc_order <- paste0(hc_order, " (", module_size, ")")


library(ggh4x)

df1 %>%
  mutate(moduleTraitCor = round(moduleTraitCor, 2)) %>%
  mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), moduleTraitCor)) %>%
  # mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), '')) %>%
  ggplot(aes(y = module, x = name, fill = moduleTraitCor)) +
  geom_tile(color = 'white', size = 0.7, width = 1) +
  # geom_raster() +
  geom_text(aes(label = star),  vjust = 0.5, hjust = 0.5, size= 2.5, family =  "GillSans") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  ggh4x::scale_y_dendrogram(hclust = hclust, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(title = "Modules", labels = hc_order, label_size = 10, label_family = "GillSans")) +
  labs(x = 'Parameters', y = "") +
  guides(fill = guide_colorbar(barwidth = unit(5, "in"),
    barheight = unit(0.1, "in"), label.position = "top",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(size = 10))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey89', color = 'white')) -> p1

p1 <- p1 + theme(
  axis.line.x = element_blank(),
  axis.line.y = element_blank(),
  axis.text.y = element_text(hjust = 1.2),
  axis.ticks.length = unit(5, "pt"),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
    margin = unit(c(t = 0.5, r = 0, b = 0, l = 0), "mm")))

p1 <- p1 + theme(panel.spacing.x = unit(-0.5, "mm"))

# p1

ggsave(p1, filename = 'moduleTraitCor_p.png', path = wd, width = 8, height = 8, device = png, dpi = 300)

# BARPLOT

stats %>% pull(n, name = module)


p2 <- stats %>% 
  mutate(module = factor(module, levels = hclust$labels[hclust$order])) %>%
  ggplot(aes(y = module)) + #  fill = DE, color = DE
  scale_x_continuous("Número de taxons") +
  geom_col(aes(x = n), width = 0.95, position = position_stack(reverse = TRUE)) +
  # geom_col(aes(x = reads_frac), width = 0.95, fill = "grey")
  # scale_fill_manual(name = '', values = c("#303960", "#647687", "#E7DFD5")) + # grey90
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.title.y = element_blank(), 
    axis.text.y= element_blank(),
    axis.ticks.y =element_blank(), 
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
      margin = unit(c(t = 0.5, r = 0, b = 0, l = 0), "mm")))


library(patchwork)

# p1 + plot_spacer() + p2 + plot_layout(widths = c(5,-0.5, 10))  #& theme(plot.margin = 0)

psave <- p1 +  plot_spacer() + p2 + plot_layout(widths = c(7, -0.25, 1.5)) + labs(caption = '* corPvalueStudent < 0.05 ') 

psave






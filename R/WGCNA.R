# GENERATE NETWORK ANALYSIS USING WGCNA
# EX. https://www.polarmicrobes.org/weighted-gene-correlation-network-analysis-wgcna-applied-to-microbial-communities/
# use also Sparse Estimation of Correlations among Microbiomes (SECOM) for correlation analysis. 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "~/MEIOFAUNA/INPUTS/"

library(phyloseq)
library(microbiome)
library(WGCNA)
library(flashClust)
library(tidyverse)


phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))

datExpr <- phyloseq %>%
  # subset_samples(Tissue=="Foregut") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  aggregate_taxa(., "p") %>%
  # transform_sample_counts(., function(x) x / sum(x)) %>%
  otu_table() %>%
  as("matrix")

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
k=softConnectivity(datE=datExpr,power=3)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

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

# sft <- read_rds(paste0(path, 'SoftThreshold_',cor_method, '.rds'))

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 2)

hist(soft_values)

power_pct <- quantile(soft_values, probs = 0.95)

softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]

meanK <- sft$fitIndices[softPower,5]

hist(sft$fitIndices[,5])

softPower <- min(softPower)

cat("\nsoftPower value", softPower, '\n')


title1 = 'Scale Free Topology Model Fit,signed R^2'
title2 = 'Mean Connectivity'

caption = paste0("Lowest power for which the scale free topology index reaches the ", power_pct*100, " %")

sft$fitIndices %>% 
  mutate(scale = -sign(slope)*SFT.R.sq) %>%
  select(Power, mean.k., scale) %>% pivot_longer(-Power) %>%
  mutate(name = ifelse(name %in% 'scale', title1, title2)) %>%
  ggplot(aes(y = Power, x = value)) +
  geom_text(aes(label = Power), size = 2) +
  geom_abline(slope = 0, intercept = softPower, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = min(meanK), linetype="dashed", alpha=0.5) +
  labs(y = 'Soft Threshold (power)', x = '', 
       caption = caption) +
  facet_grid(~name, scales = 'free_x', switch = "x") +
  # scale_x_continuous(position = "top") +
  theme_light(base_family = "GillSans",base_size = 16) -> psave

psave


# 2) 
# Build a adjacency "correlation" matrix

enableWGCNAThreads()

# specify network type
# softPower <- 2

adjacency <- adjacency(datExpr, 
                       power = softPower, 
                       type = "unsigned")


TOM <- TOMsimilarity(adjacency) # specify network type


dissTOM = 1 - TOM

# heatmap(dissTOM) # if more than 100 rows, omit

# Generate Modules ----
# Generate a clustered gene tree

geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", 
     main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

# This sets the minimum number of genes to cluster into a module

minClusterSize <- 2

cutHeight <- 0.6

dynamicMods <- cutreeDynamic(dendro= geneTree, 
                             distM = dissTOM,
                             method = "hybrid",
                             deepSplit = 2, 
                             cutHeight = cutHeight,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minClusterSize)

dynamicColors = labels2colors(dynamicMods)
names(dynamicColors) <- colnames(datExpr)
MEList = moduleEigengenes(datExpr, 
                          colors= dynamicColors,
                          softPower = softPower)

MEs = MEList$eigengenes

MEDiss= 1 - cor(MEs)

METree = flashClust(as.dist(MEDiss), method= "average")

# Set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0

MEDissThres = 0.8

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors = merge$colors

mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, dynamicColors, c("Modules"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

caption <- c("Dynamic Tree Cut", "Merged dynamic", "\n(cutHeight: ", cutHeight,")")

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), caption, dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors <- mergedColors

colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1


MEs <- mergedMEs


# how modules where obtained:
nm <- table(moduleColors)


cat('Number of mudules obtained\n :', length(nm))


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

datTraits <- phyloseq %>% sample_data() %>% with(., table(LIBRARY_ID, Region))

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

df1 %>%
  mutate(star = ifelse(corPvalueStudent <.001, "***", 
                       ifelse(corPvalueStudent <.01, "**",
                              ifelse(corPvalueStudent <.05, "*", "")))) -> df1


df1 %>%
  mutate(moduleTraitCor = round(moduleTraitCor, 2)) %>%
  mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), moduleTraitCor)) %>%
  # mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), '')) %>%
  ggplot(aes(y = module, x = name, fill = moduleTraitCor)) +
  geom_tile(color = 'white', size = 0.7, width = 1) +
  # geom_raster() +
  geom_text(aes(label = star),  vjust = 0.5, hjust = 0.5, size= 4, family =  "GillSans") +
  ggsci::scale_fill_gsea(name = "", reverse = T, na.value = "white") +
  # scale_fill_viridis_c(name = "Membership", na.value = "white") +
  ggh4x::scale_y_dendrogram(hclust = hclust) +
  labs(x = '', y = 'Module') +
  guides(fill = guide_colorbar(barwidth = unit(3.5, "in"),
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
  axis.ticks.length = unit(5, "pt"))

p1 <- p1 + theme(panel.spacing.x = unit(-0.5, "mm"))

p1

# BARPLOT

reads <- colSums(datExpr)

Total <- sum(reads)

# Sanity check:

identical(names(colSums(datExpr)), names(moduleColors))

data.frame(reads, moduleColors) %>% 
  as_tibble(rownames = "Name") %>% 
  group_by(moduleColors) %>% 
  summarise(n = n(), reads = sum(reads)) %>%
  dplyr::rename('module' = 'moduleColors') -> stats

p2 <- stats %>% 
  mutate(module = factor(module, levels = hclust$labels[hclust$order])) %>%
  ggplot(aes(y = module)) + #  fill = DE, color = DE
  scale_x_continuous("Número de miRs") +
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
        axis.ticks.length = unit(5, "pt"))


library(patchwork)

# p1 + plot_spacer() + p2 + plot_layout(widths = c(5,-0.5, 10))  #& theme(plot.margin = 0)

psave <- p1 +  plot_spacer() + p2 + plot_layout(widths = c(7, -0.25, 1.5)) + labs(caption = '* corPvalueStudent < 0.05 ') 

psave

save.image()


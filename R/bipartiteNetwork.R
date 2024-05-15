
# Network analysis represents a powerful way to visual- ize and analyze the association between species and their habitat, and it has become popular in ecology due to the prevalence of networks in nature (Pavlopoulos et al. 2018). Pollinator–plant associations, parasite–host relationships, plant–animal mutualisms, and food webs have been success- fully modeled with network structures, resulting in an abun- dance of quantitative approaches to describe networks (Bas- compte et al. 2003; Casey et al. 2019; Olesen et al. 2007; Pozas-Schacre et al. 2021; Runghen et al. 2021). Two of the most commonly estimated network metrics are nestedness, the loss of interactions in smaller subsets of ecological net- works, and modularity, the division of networks into smaller interacting clusters (Newman 2006; Payrató-Borràs et al. 2020; Rodríguez-Gironés and Santamaría, 2006). These net- work metrics are ideally suited to measure specialization in a goby-sponge network because they can identify the occur- rence of highly asymmetrical and non-random interactions between two groups (Bascompte et al. 2003; Newman 2006).

# We computed nestedness using the “nested” function from the bipartite package for both the null models and our empirical networks. Nestedness assesses the relative occur- rence of generalist nodes, which are linked to a wide range of other nodes, and specialist nodes, which are linked to a subset of those nodes (Pavlopoulos et al. 2018) ... We computed modularity using the “computeModules” function from the bipartite package. Modularity measures how strongly a network clusters into groups of interacting nodes.

library(bipartite)

data(Safariland)


class(Safariland)

plotweb(Safariland, col.high = "salmon", col.low = "violet") 

# own data

path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/"


file_out <- "WGCNA.p.taxa"

g <- read_rds(paste0(path, "/",file_out, ".g.rds"))

ps <- read_rds(paste0(path, '/ps.rds'))

library(phyloseq)

infoNodes <- as_tibble(as(tax_table(ps), "matrix"), rownames = "nodeName")

g <- g %>% activate("nodes") %>% left_join(infoNodes, by = c('nodeName'))

agglom_lev <- "p"

web <- g %>% activate("nodes") %>% as_tibble() %>%
  filter(k != "Plantae") %>%
  dplyr::rename( "Level" = all_of(agglom_lev)) %>%
  dplyr::rename( "Module" = all_of("nodeAttr[nodesPresent, ]")) %>%
  mutate(Level = ifelse(as.numeric(Confidence) >= 0.5, Level, "No hit")) %>%
  mutate(Level = ifelse(!is.na(Level), Level, "No hit")) %>%
  group_by(Module, Level) %>%
  count(sort = T) %>%
  filter(Level != "No hit") %>%
  pivot_wider(names_from = Level, values_from = n, values_fill = 0) %>%
  data.frame(row.names = .$Module) %>%
  select(-Module) %>% as("matrix")

plotweb(web, col.high = "salmon", col.low = "violet", method = "cca", text.rot=90) 

visweb(web, type = "nested") 

networklevel(web) # index = "nestedness"

# indexes
# Conectancia: Demuestra los posibles enlaces que se puede dar en la red indicando la proporción de interacciones realizadas entre especies.
# 
# Robustez: Evalua la estabilidad de la red, su tolerancia ante perturbaciones.
# 
# Anidamiento: Describe cualitativamente y cuantitativamente el cómo se mantienen las relaciones entre especies en una comunidad.
# 
# Especialización H2´: Indica que tan diversas son esas interacciones considerando toda la matriz y cada especie.
# 
# Asimetría: Fuerza e interaccion entre los 2 niveles tróficos. Osea el nivel de dependencia del nivel superior con el inferior. En este caso del polinizador con la planta.
# 
# Diversidad de Shannon: Cuantifica biodiversidad específica reflejando heterogeneidad en las comunidades analizadas considerando número de especies y abundancia relativa.

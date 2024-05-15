# network ====

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(igraph)
library(tidygraph)
library(ggraph)
library(tidyverse)

# path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/"
  
path <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/raw-seqs-bkp/filtN/cutadapt/Illumina/filterAndTrim/"

x <- readr::read_rds(file.path(path, "WGCNA.rds"))

exportNet <- function(TOMobj, moduleColors, threshold = 0.9) {
  
  nodeNames <- names(moduleColors)
  nodeAttr <- moduleColors
  
  
  # file_out <- gsub('.RData', '', basename(TOMFile))
  
  file_out <- "WGCNA.taxa"
  
  edgeFile <- paste0(path, "/",file_out, ".edges.txt")
  nodeFile <- paste0(path, "/",file_out, ".nodes.txt")
  
  
  Net = WGCNA::exportNetworkToCytoscape(TOMobj,
                                 edgeFile = edgeFile,
                                 nodeFile = nodeFile,
                                 weighted = TRUE,
                                 threshold = threshold,
                                 nodeNames = nodeNames, nodeAttr = nodeAttr)
  
  cat('\nEdges: ',nrow(Net$edgeData), 'and Nodes: ',nrow(Net$nodeData),'\n')
  
  graph <- tidygraph::tbl_graph(nodes = Net$nodeData, edges = Net$edgeData, directed = FALSE)
  
  return(graph)
  
  # if(is.null(query)) {
  #   graph
  # } else {
  #   graph %>% activate(nodes) %>% filter(nodeName %in% query)
  # }
  
  
  # return(Net)
  
}

# run once 

TOM <- x$TOM

moduleColors <- x$moduleColors

# g <- exportNet(TOM, moduleColors, threshold = 0.5)

file_out <- "WGCNA.p.taxa"

edgeFile <- paste0(path, "/",file_out, ".edges.txt")
nodeFile <- paste0(path, "/",file_out, ".nodes.txt")


g <- read_rds(paste0(path, "/",file_out, ".g.rds"))


g <- g %>% activate("nodes") %>%  
  mutate(betweenness = betweenness(.), degree = centrality_degree(),
    membership = igraph::cluster_louvain(., igraph::E(g)$weight)$membership,
    pageRank = page_rank(.)$vector)


ps <- read_rds(paste0(path, '/ps.rds'))

library(phyloseq)

infoNodes <- as_tibble(as(tax_table(ps), "matrix"), rownames = "nodeName")

g <- g %>% activate("nodes") %>% left_join(infoNodes, by = c('nodeName'))

agglom_lev <- "p"

# plot n asvs (by tax group) per module

g %>% activate("nodes") %>% as_tibble() %>%
  filter(k != "Plantae") %>%
  dplyr::rename( "Level" = all_of(agglom_lev)) %>%
  dplyr::rename( "Module" = all_of("nodeAttr[nodesPresent, ]")) %>%
  mutate(Level = ifelse(as.numeric(Confidence) >= 0.5, Level, "No hit")) %>%
  mutate(Level = ifelse(!is.na(Level), Level, "No hit")) %>%
  group_by(Module, Level) %>%
  count(sort = T) %>%
  # pivot_wider(names_from = Level, values_from = n, values_fill = 0) %>% view()
  filter(Level != "No hit") %>%
  mutate(Module = factor(Module, levels = hc_order)) %>%
  ggplot(aes(x = Module, y = Level, fill = n)) +
  scale_fill_viridis_c() +
  geom_tile(color = 'white', size = 0.7, width = 1) +
  labs(x = 'Modules', y = "Phylum") +
  guides(fill = guide_colorbar(title = "Number of ASVs", barwidth = unit(5, "in"),
    barheight = unit(0.1, "in"), label.position = "top",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(size = 10))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(hjust = 1.2),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
      margin = unit(c(t = 0.5, r = 0, b = 0, l = 0), "mm")),
    panel.spacing.x = unit(-0.5, "mm")) -> psave


# plot n edges per tax group
get_summary_stats_df <- g %>% activate("nodes") %>% as_tibble() %>%
  filter(k != "Plantae") %>%
  dplyr::rename( "Level" = all_of(agglom_lev)) %>%
  dplyr::rename( "Module" = all_of("nodeAttr[nodesPresent, ]")) %>%
  mutate(Level = ifelse(as.numeric(Confidence) >= 0.5, Level, "No hit")) %>%
  mutate(Level = ifelse(!is.na(Level), Level, "No hit")) %>%
  filter(degree > 0) %>%
  group_by(Level) %>%
  select(degree, betweenness) %>%
  # summarise(n = sum(degree)) 
  rstatix::get_summary_stats(type = c("common"))


get_summary_stats_df %>%
  filter(Level != "No hit") %>%
  ggplot(aes(x = Level, y = mean)) +
  facet_grid(variable ~ ., scales = "free_y") +
  geom_col() +
  labs(y = "Mean connectivity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
    margin = unit(c(t = 0.5, r = 0, b = 0, l = 0), "mm")))

# p1

ggsave(psave, filename = 'moduleTraitCor_tax.png', path = wd, width = 8, height = 8, device = png, dpi = 300)



# g <- tidygraph::tbl_graph(nodes = read.delim(nodeFile), edges = read.delim(edgeFile), directed = FALSE)

# g %>% activate("edges")  %>% as_tibble() %>% ggplot(aes(weight)) + geom_histogram()

# g %>% activate("edges")  %>% as_tibble() %>% ggplot(aes(weight)) + stat_ecdf()

g <- g %>% activate("edges") %>% 
  mutate(weight = ifelse(weight < 0.9, NA, weight)) %>% filter(!is.na(weight))

# g <- g %>% activate("edges") %>% 
#   mutate(weight = ifelse(weight < 0.9, NA, weight)) %>% filter(!is.na(weight))
# 

# write_rds(g, file = paste0(path, "/",file_out, ".g.rds"))

# g %>% activate("edges")  %>% as_tibble() %>% ggplot(aes(weight)) + geom_histogram()



str(scale_col <- g %>% activate("nodes") %>% as_tibble() %>% distinct(`nodeAttr[nodesPresent, ]`) %>% pull())

g %>% activate("nodes")  %>% as_tibble() %>% ggplot(aes(betweenness, degree)) + geom_point()

g %>% activate("nodes")  %>% as_tibble() %>% # count(betweenness) %>%
  ggplot(aes(betweenness)) + 
  geom_histogram()

  
# g %>% activate("nodes") %>% as_tibble() %>%
#   # ggplot(aes(degree)) + geom_histogram()
#   group_by(`nodeAttr[nodesPresent, ]`) %>%
#   count(degree) %>%
#   ggplot(aes(y = n, x = degree, color = `nodeAttr[nodesPresent, ]`)) +
#   geom_point() +
#   geom_line(orientation = "x") +
#   facet_wrap(~ `nodeAttr[nodesPresent, ]`) +
#   theme_bw(base_family = "GillSans", base_size = 14) +
#   scale_color_manual(values = structure(scale_col, names = scale_col)) +
#   theme(legend.position = "none") +
#   labs(x = "Centrality Degree", y = "Number of ASVs")
  


# .g <- g

g <- g %>% activate("nodes") %>% filter(betweenness > 1)


Levels <- scale_col

g <- g %>% activate("nodes") %>%
  mutate(`nodeAttr[nodesPresent, ]` = factor(`nodeAttr[nodesPresent, ]`, levels = Levels))

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')



ggraph(layout) +
  facet_nodes(~ `nodeAttr[nodesPresent, ]`) +
  geom_node_point(aes(color = k)) +
  geom_edge_arc(aes(edge_alpha = weight), strength = 0.1) + # edge_width
  # scale_edge_width(range = c(0.3, 1)) +
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top") 
  # coord_fixed() 

ggraph(layout) +
  scale_color_manual('', values = structure(scale_col, names = scale_col) ) +
  geom_node_point(aes(color = `nodeAttr[nodesPresent, ]`, size = degree)) +
  geom_edge_arc(aes(edge_alpha = weight), strength = 0.1) + # edge_width
  # ggrepel::geom_text_repel(data = layout, aes(x, y, label = nodeName), max.overlaps = 50, family = "GillSans") +
  scale_edge_width(range = c(0.3, 1)) +
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top") +
  coord_fixed() -> psleft

ggsave(psleft, filename = 'WGCNA_TAX2NETWORK.png', path = path, width = 10, height = 10, device = png, dpi = 300)


# quit()

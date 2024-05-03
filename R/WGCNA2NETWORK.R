# network ====
# 
# module_members_1 <- c("grey", "black", "yellow", "pink", "turquoise")
# 
# module_members_2 <- c("brown", "red", "green", "blue")

library(igraph)
library(tidygraph)
library(ggraph)

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

# g <- exportNet(TOM, moduleColors, threshold = 0)

file_out <- "WGCNA.p.taxa"

edgeFile <- paste0(wd, "/",file_out, ".edges.txt")
nodeFile <- paste0(wd, "/",file_out, ".nodes.txt")


g <- tidygraph::tbl_graph(nodes = read.delim(nodeFile), edges = read.delim(edgeFile), directed = FALSE)

g %>% activate("edges")  %>% as_tibble() %>% ggplot(aes(weight)) + geom_histogram()

g <- g %>% activate("edges") %>% mutate(weight = ifelse(weight < 0.3, NA, weight)) %>% filter(!is.na(weight))

g %>% activate("edges")  %>% as_tibble() %>% ggplot(aes(weight)) + geom_histogram()

scale_col <- g %>% activate("nodes") %>% as_tibble() %>% distinct(`nodeAttr.nodesPresent...`) %>% pull()

g <- g %>% activate("nodes") %>%  
  mutate(betweenness = betweenness(.), degree = centrality_degree(),
         membership = igraph::cluster_louvain(., igraph::E(g)$weight)$membership,
         pageRank = page_rank(.)$vector)


g %>% activate("nodes") %>% as_tibble() %>% 
  # ggplot(aes(degree)) + geom_histogram()
  group_by(`nodeAttr.nodesPresent...`) %>%
  count(degree) %>%
  ggplot(aes(y = n, x = degree, color = `nodeAttr.nodesPresent...`)) + 
  geom_point() +
  geom_line(orientation = "x") +
  facet_wrap(~ `nodeAttr.nodesPresent...`) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  scale_color_manual(values = structure(scale_col, names = scale_col)) +
  theme(legend.position = "none") +
  labs(x = "Centrality Degree", y = "Number of ASVs")
  

Levels <- scale_col

g <- g %>% activate("nodes") %>%
  mutate(`nodeAttr.nodesPresent...` = factor(`nodeAttr.nodesPresent...`, levels = Levels))

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')


ggraph(layout) +
  scale_color_manual('', values = structure(scale_col, names = scale_col) ) +
  geom_edge_arc(aes(edge_alpha = weight), strength = 0.1) + # edge_width
  geom_node_point(aes(color = `nodeAttr.nodesPresent...`, linewidth = degree)) +
  # ggrepel::geom_text_repel(data = layout, aes(x, y, label = nodeName), max.overlaps = 50, family = "GillSans") +
  scale_edge_width(range = c(0.3, 1)) +
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "top") +
  coord_fixed() -> psleft

ggsave(psleft, filename = 'WGCNA_TAX2NETWORK.png', path = path, width = 10, height = 10, device = png, dpi = 300)


quit

node_labs <- paste(LETTERS[1:length(Levels)], ") ", Levels, sep = "")
names(node_labs) <- Levels


ggraph(layout) +
  facet_nodes(~ `nodeAttr.nodesPresent...`, labeller = labeller(`nodeAttr.nodesPresent...` = node_labs), 
              scales = "free") +
  scale_color_manual('', values = structure(scale_col, names = scale_col) ) +
  geom_edge_arc(aes(edge_alpha = weight), strength = 0.1) + # edge_width
  geom_node_point(aes(color = `nodeAttr.nodesPresent...`, size = degree)) +
  ggrepel::geom_text_repel(data = layout, aes(x, y, label = nodeName), max.overlaps = 50, family = "GillSans", size = 2) +
  scale_edge_width(range = c(0.3, 1)) +
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  # coord_fixed() +
  guides(color = "none") +
  ggforce::geom_mark_hull(
    aes(x, y, group = as.factor(`nodeAttr.nodesPresent...`)),
    color = NA, fill = "grey76",
    concavity = 4,
    con.size = 0.3,
    con.linetype = 2,
    expand = unit(2, "mm"),
    alpha = 0.25)  +
  guides(fill = FALSE) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = 'grey89', color = 'white')) -> ps


#
library(patchwork)

psave <- psleft + ps + patchwork::plot_layout(widths = c(0.7, 1.2))

ggsave(psave, filename = 'WGCNA_MIRS2NETWORK_2.png', path = path, width = 12, height = 8, device = png, dpi = 300)

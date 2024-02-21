# Generate a correlation sample-sample heatmap

wd <- "~/MEIOFAUNA/INPUTS/"

MTD <-  read_tsv(list.files(path = wd, pattern = "mapping-file-corregido.tsv", full.names = T)) %>% 
  select(`#SampleID`, Depth, Region) %>%
  dplyr::rename("LIBRARY_ID" = "#SampleID") %>%
  mutate(Region = factor(Region, levels = c("Yucatan", "NW Shelf", "NW Slope", "Deep-sea")))

ab_f <- list.files(path = wd, pattern = 'table_100_80', full.names = T)

DATA <-  read_tsv(ab_f, skip = 1) 

DATA <- DATA %>% select_if(is.double) 

DATA <- as(DATA, "matrix")

sample_cor = cor(DATA, method='pearson', use='pairwise.complete.obs')

h <- heatmap(sample_cor, col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Rowv)

hc_order <- hc_samples$labels[h$rowInd]

#

sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  left_join(MTD) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order))) -> sample_cor_long

library(ggh4x)

P <- sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) +  
  geom_tile(linewidth = 0.2) +
  ggsci::scale_fill_material(name = "", "indigo") +
  scale_x_discrete(position = 'bottom') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = hc_order, label_size = 2.5, label_family = "GillSans")) +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) +
  theme_bw(base_size = 7, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "bottom",
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 2.5),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) 

P <- P + guides(
  fill = guide_colorbar(barwidth = unit(1.5, "in"),
                        barheight = unit(0.05, "in"), label.position = "bottom",
                        alignd = 0.5,
                        ticks.colour = "black", ticks.linewidth = 0.5,
                        frame.colour = "black", frame.linewidth = 0.5,
                        label.theme = element_text(family = "GillSans", size = 7)))


# TOP

TOPDF <- sample_cor_long %>%
  distinct(LIBRARY_ID, Region, Depth  ) %>%
  # dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  # mutate(label = ifelse(pH %in% "Low", "*", "")) %>%
  mutate(y = 1)


color_vector <- as.character(unique(MTD$Region))

getPalette <- c("#4DAF4A", "#313695", "lightblue", "#E41A1C")

axis_col <- structure(getPalette, names = color_vector)


topplot <- TOPDF %>%
  ggplot(aes(y = y, x = LIBRARY_ID, color = as.factor(Depth))) +
  geom_point(shape = 15, size = 2) +
  #geom_text(aes(label = label),  vjust = -0.7, hjust = 0.5, size = 1.5, family =  "GillSans", color = "#d73027") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  # ggh4x::guide_dendro()
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)) +
  guides(color = guide_legend(title = "", ncol = 4 )) +
  theme_bw(base_family = "GillSans", base_size = 7) +
  scale_color_manual(values = axis_col ) +
  theme(legend.position = 'top',
        panel.border = element_blank(), 
        plot.background = element_rect(fill='transparent', color = 'transparent'),
        plot.margin = unit(c(0,0,0,0), "pt"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank()) 

library(patchwork)

psave <- topplot/ plot_spacer() /P + plot_layout(heights = c(0.6, -0.5, 5))

ggsave(psave, filename = 'SAMPLE_HEATMAP.png', path = wd, width = 4, height = 4, device = png, dpi = 300)


color_vector <- as.character(unique(MTD$Region))

getPalette <- c("#4DAF4A", "#313695", "lightblue", "#E41A1C")

axis_col <- structure(getPalette, names = color_vector)

MTD %>% count(Depth) %>% pull(Depth)

depthL <- c("0-100","100-200", "200-500",
            "500-1000", "1000-1500", "1500-2000",
            "2000-2500","2500-3000","3000-3500","3500-4000") 

MTD %>% count(Depth , Region) %>% 
  mutate(Depth = factor(Depth, levels = depthL)) %>%
  ggplot(aes(x = Depth, y = n,color = Region)) +
  geom_point(size = 4, shape = 1) +
  # ggrepel::geom_label_repel(aes(label = n, size = 1)) +
  guides(color = guide_legend(title = "", ncol = 4 )) +
  theme_classic(base_family = "GillSans", base_size = 7) +
  scale_color_manual(values = axis_col ) +
  theme(legend.position = 'top') +
  labs(y = "N samples") -> p

ggsave(p, filename = 'SAMPLES-DEPTH.png', path = wd, width = 5, height = 4, device = png, dpi = 300)

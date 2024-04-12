# GENERATE MAPS



rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "/Users/cigom/Documents/MEIOFAUNA_PAPER/RDADA2-OUTPUT/"



reg_levels <- c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")
getPalette <- c("#000056", "#2E71A7","#60A4CF", "#9ECAE1")
axis_col <- structure(getPalette, names = reg_levels)

library(phyloseq)
library(tidyverse)

phyloseq <- read_rds(paste0(wd, '/phyloseq.rds'))

DB <- sample_data(phyloseq) %>% as(., "matrix") %>% as_tibble()

station <- c('LIBRARY_ID','Sampling.station',	'Cruise',	'Region')
coords <- c('Latitude','Longitude', 'Depth')
vars <- c("Observed", "Chao1", "Shannon", "InvSimpson")


cols <- c(station, coords, vars)

COORDS <- DB %>% 
  select_at(vars(all_of(contains(cols)))) %>% 
  mutate_at(vars(contains(c(coords, vars))), as.double) %>%
  drop_na(Latitude,Longitude)


# COORDS %>% summarise(minLon = min(Longitude), maxLon = max(Longitude),
#   )
# 


limits <- c(min(COORDS$Longitude), max(COORDS$Longitude), 
          min(COORDS$Latitude), max(COORDS$Latitude))

# https://mikkovihtakari.github.io/ggOceanMaps/articles/ggOceanMaps.html

library(ggOceanMaps)

#limits are given longitude min/max, latitude min/max

# #EEEAC4

p <- basemap(
  limits = limits,
  bathymetry = TRUE, 
  # bathy.style = "rcb",
  bathy.style = "contour_blues",
  land.col = "grey88") + # "#EEEAC4"
  ggspatial::annotation_scale(location = "br", text_family = "GillSans") + 
  # ggspatial::annotation_north_arrow(location = "tr", which_north = "true") +
  labs(x = "", y = "") +
  # scale_y_continuous() + scale_x_continuous() +
  guides(fill = guide_legend(title="Bathymetry (m)",  nrow = 1),
         color = guide_legend(title="Bathymetry (m)",  nrow = 1)) +
  theme(
    text = element_text(family = "GillSans"))

 
c("Deep-sea", "NW Slope", "NW Shelf", "Yucatan")

p1 <- p +  
  ggspatial::geom_spatial_point(data = subset(COORDS, Region == "Deep-sea"), 
  aes(x = Longitude, y = Latitude, shape = Cruise), 
    shape = 21, size = 1.5, stroke = 1.2, fill = "white", color = "#2426FA") +
  facet_grid(~ Region)

p2 <- p +  
  ggspatial::geom_spatial_point(data = subset(COORDS, Region == "NW Slope"), 
    aes(x = Longitude, y = Latitude, shape = Cruise), 
    shape = 21, size = 1.5, stroke = 1.2, fill = "white", color = "#2426FA") +
  facet_grid(~ Region)


p3 <- p +  
  ggspatial::geom_spatial_point(data = subset(COORDS, Region == "NW Shelf"), 
    aes(x = Longitude, y = Latitude, shape = Cruise), 
    shape = 21, size = 1.5, stroke = 1.2, fill = "white", color = "#2426FA") +
  facet_grid(~ Region)

p4 <- p +  
  ggspatial::geom_spatial_point(data = subset(COORDS, Region == "Yucatan"), 
    aes(x = Longitude, y = Latitude, shape = Cruise), 
    shape = 21, size = 1.5, stroke = 1.2, fill = "white", color = "#2426FA") +
  facet_grid(~ Region) 



library(patchwork)

pmap <- p1 +p2 + p3+ p4 + 
  patchwork::plot_layout(guides = "collect", nrow = 1) & 
  theme_bw(base_family = "GillSans") +
  theme(legend.position='bottom', legend.justification = "right",
    panel.grid = element_blank(), legend.text = element_text(size = 7),
    strip.background = element_rect(fill = 'grey88'))

# scale_color_manual(values = axis_col ) +
# theme(legend.position = "top") 

ggsave(pmap, path = wd, filename = 'maps-region-panel.png', width = 10, height = 3.5, device = png, dpi = 300)


basemap(
  limits = limits,
  bathymetry = TRUE, 
  bathy.style = "poly_blues",
  land.col = "#eeeac4", gla.col = "cadetblue") +
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") +
  labs(x = "", y = "") +
  # scale_y_continuous() + scale_x_continuous() +
  guides(fill=guide_legend(title="Bathymetry (m)")) +
  theme(
    text = element_text(family = "GillSans"))

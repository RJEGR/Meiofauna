# GENERATE MAPS


# https://mikkovihtakari.github.io/ggOceanMaps/articles/ggOceanMaps.html

library(ggOceanMaps)
#limits are given longitude min/max, latitude min/max
p <- basemap(limits = c(-118, -100, 20, 43),
  bathymetry = TRUE, 
  bathy.style = "rcb",
  land.col = "grey76") +
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") +
  labs(x = "Longitud", y = "Latitud") +
  guides(fill=guide_legend(title="Profundidad (m)")) +
  theme(
    text = element_text(family = "GillSans"))

P +  ggspatial::geom_spatial_point(data = dt, aes(x = lon, y = lat), color = "red")


ggsave(p, path = path, filename = 'FIGURE_4_TESIS_1.png', width = 10, height = 3.5, device = png, dpi = 300)

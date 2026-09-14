source("1.load_data.R")

#barplot

barplot_metabarcoding <- abundance_bar_plot(
  table = table_fung,
  metadata = metadata_fungi,
  taxonomy_db = "silva",
  level = "phylum", 
  x_col = "Type_of_soil",
  label = "Phylum", 
  x_axis_title = "Source", 
  add_remained = TRUE
)

barplot_metagenomic <- abundance_bar_plot(
  table = table_meta,
  metadata = metadata_meta,
  taxonomy_db = "Kraken2",
  level = "phylum",
  x_col = "Polygon",
  x_axis_title = "Polygons",
  label = "Phylum",
  add_remained = TRUE)


both <- cowplot::plot_grid(barplot_metabarcoding, barplot_metagenomic, labels = "AUTO")

ggsave("plots/barplots.png", width = 10, height = 6, units = "in", dpi = 300)


#heatmaps
metadata_fungi1 <- metadata_fungi[
  order(metadata_fungi$Type_of_soil),
]

sample_order <- metadata_fungi1$SAMPLEID

table_fung1 <- table_fung[, c(
  sample_order[sample_order %in% colnames(table_fung)],
  "taxonomy"
)]



heat_metabarcoding <- abundance_heatmap_plot(
  table = table_fung1,
  metadata = metadata_fungi, condition1 = "Type_of_soil",  
  show_column_names = FALSE,

  top_n = 20)

# Ordenar metadata por Polygon
metadata_meta1 <- metadata_meta[
  order(metadata_meta$Polygon),
]

# Obtener el orden de las muestras
sample_order <- metadata_meta1$SAMPLEID

# Ordenar table_meta por Polygon y dejar taxonomy al final
table_meta1 <- table_meta[, c(
  sample_order[sample_order %in% colnames(table_meta)],
  "taxonomy"
)]
heat_metagenomic <- abundance_heatmap_plot(
  table = table_meta1,
  metadata = metadata_meta1,
  condition1 = "Polygon",
  cluster = FALSE,
  show_column_names = FALSE,
  top_n = 20
)

both2 <- cowplot::plot_grid(heat_metabarcoding, heat_metagenomic, labels = "AUTO", rel_widths = c(1.5,1))

ggsave("plots/heats.png", width = 16, height = 6, units = "in", dpi = 300)


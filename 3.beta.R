source("1.load_data.R")

# "compositional" usa una instancia Monte Carlo de ALDEx2 al azar:
# set.seed() antes de cada llamada para obtener siempre el mismo resultado.
set.seed(123)
beta_fungi <- beta_ord_plot(table = table_fung,
              metadata = metadata_fungi,
              group_col = "Source",
              shape_col = "Treatment",
              distance = "compositional",
              ordination = "PCA")

set.seed(123)
betaperm_fungi <- beta_test_table(table = table_fung,
                metadata= metadata_fungi,
                formula_str = "Source",
                distance = "compositional",
                test = "permanova",
                permutations = 999)



beta_meta <- beta_ord_plot(table = table_meta,
              metadata = metadata_meta,
              distance = "jaccard",
              ordination = "PCoA",
              group_col = "Polygon")

beta_perm_meta <- beta_test_table(table = table_meta,
                metadata = metadata_meta,
                formula_str = "Polygon",
                distance = "jaccard",
                test = "permanova",
                permutations = 999)



beta<- cowplot::plot_grid(beta_fungi, betaperm_fungi,
                          beta_meta,  beta_perm_meta,
                          labels = c("A", "B", "C", "D"),
                          label_y = 1, ncol = 2, nrow = 2, align = "hv")

ggsave("plots/beta.png", width = 12, height = 10, units = "in", dpi = 300)


beta_decay_plot(
  table    = table_meta,
  metadata = metadata_meta,
  lat_col  = "Latitude",
  lon_col  = "Longitude",
  distance = "jaccard"
)

betapart <- beta_partition_ord_plot(
  table = table_fung,
  metadata = metadata_fungi,
  family = "jaccard",
  group_col = "Source",
  shape_col = "Treatment", 
  save_table = FALSE)


# --- Figura 1: con que compartimento comparte mas la raiz --------------------
# Pares ENTRE compartimentos (Rhizosphere-Roots y Bulk soil-Roots), con pares
# de muestras del mismo tratamiento (una faceta por tratamiento).
# beta_dissimilarity_plot, partition = "shared" (presencia/ausencia): cuantas
# features comparte cada par. Se espera que la raiz comparta mas con la
# rizosfera, la zona de suelo que modifica, que con el suelo lejano.
# Resultado: si (mediana Rhizosphere-Roots 3, 3, 6 contra Bulk soil-Roots
# 0, 1, 1 en TC, TD, TED), y en TED comparte el doble con la rizosfera.

betadis <- beta_dissimilarity_plot(
  table                 = table_fung,
  metadata              = metadata_fungi,
  comparison_condition1 = c("Rhizosphere_vs_Roots", "Bulk soil_vs_Roots"),
  condition1_col        = "Source",
  condition2_col        = "Treatment",
  group_colors          = c("Rhizosphere_vs_Roots" = "#E69F00",
                            "Bulk soil_vs_Roots"   = "#56B4E9"),
  partition             = "shared",
  family                = "jaccard",
  x_axis_title          = "Section",
  show_x_labels         = FALSE,
  stat                  = "wilcox.test")
betadis


# --- Figura 2: la raiz tambien comparte los dominantes con la rizosfera? -----
# Las mismas comparaciones que la Figura 1, ahora con beta_turnover_plot
# (numeros de Hill): q = 0 cuenta todos los taxones igual; q = 1 y q = 2 dan
# mas peso a los comunes y dominantes. Menos recambio = mas parecidas.
# Resultado: Rhizosphere-Roots siempre por debajo de Bulk soil-Roots (la raiz
# se parece mas a la rizosfera). En TC, con q = 1 y q = 2 las dos quedan en
# ~1: los dominantes de la raiz son propios. En TD y TED, Rhizosphere-Roots
# baja a ~0.6-0.75: la raiz comparte dominantes con la rizosfera.

betaturn <- beta_turnover_plot(
  table                 = table_fungi,
  metadata              = metadata_fungi,
  comparison_condition1 = c("Rhizosphere_vs_Roots", "Bulk soil_vs_Roots"),
  condition1_col        = "Source",
 # condition2_col        = "Treatment",
  group_colors          = c("Rhizosphere_vs_Roots" = "#E69F00",
                            "Bulk soil_vs_Roots"   = "#56B4E9"),
  x_axis_title          = "Section",
  show_x_labels         = FALSE,
  stat                  = "wilcox.test")
betaturn

# Figura compuesta: A = shared, B = recambio con Hill; una sola leyenda abajo
# (es la misma en los dos paneles)
leyenda <- cowplot::get_legend(betadis + ggplot2::theme(legend.position = "bottom"))
paneles <- cowplot::plot_grid(
  betadis  + ggplot2::theme(legend.position = "none"),
  betaturn + ggplot2::theme(legend.position = "none"),
  labels = c("A", "B"), rel_widths = c(1.5, 1))
betafig <- cowplot::plot_grid(paneles, leyenda, ncol = 1, rel_heights = c(1, 0.08))
betafig

ggsave("plots/betaturn.png", betafig, width = 14, height = 7, units = "in", dpi = 300)

# Los p-valores de las dos figuras comparan pares de muestras, que no son
# independientes (cada muestra entra en muchos pares): se reportan como
# descriptivos. La prueba formal entre grupos es la PERMANOVA de beta_test_table().


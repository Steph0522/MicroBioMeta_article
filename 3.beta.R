source("1.load_data.R")

# structure and beta diversity

set.seed(123)
beta_fungi <- beta_ord_plot(table = table_fung,
              metadata = metadata_fungi,
              group_col = "Source",
              shape_col = "Treatment",
              distance = "compositional",
              mc_samples = 128,
              ordination = "PCA")

set.seed(123)
betaperm_fungi <- beta_test_table(table = table_fung,
                metadata= metadata_fungi,
                formula_str = "Source*Treatment",
                distance = "compositional",
                mc_samples = 128,
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



beta<- cowplot::plot_grid(beta_fungi, ggplot2::autoplot(betaperm_fungi),
                          beta_meta,  ggplot2::autoplot(beta_perm_meta),
                          labels = c("A", "B", "C", "D"),
                          label_y = 1, ncol = 2, nrow = 2, align = "hv")

#ggsave("plots/beta.png", width = 12, height = 11, units = "in", dpi = 300)

# shared and turnover



betapart <- beta_partition_ord_plot(
  table = table_fung,
  metadata = metadata_fungi,
  family = "jaccard",
  group_col = "Source",
  shape_col = "Treatment", 
  save_table = FALSE)




betadis <- beta_dissimilarity_plot(
  table                 = table_fung,
  metadata              = metadata_fungi,
  comparison_condition1 = c("Rhizosphere_vs_Roots", "Bulk soil_vs_Roots"),
  condition1_col        = "Source",
  condition2_col        = "Treatment",
  group_colors          = c("Rhizosphere_vs_Roots" = "#56B4E9",
                            "Bulk soil_vs_Roots"   = "#E69F00"),
  partition             = "shared",
  family                = "jaccard",
  x_axis_title          = "Section",
  show_x_labels         = FALSE,
  stat                  = "wilcox.test")
betadis


betaturn <- beta_turnover_plot(
  table                 = table_fungi,
  metadata              = metadata_fungi,
  comparison_condition1 = c("Rhizosphere_vs_Roots", "Bulk soil_vs_Roots"),
  condition1_col        = "Source",
  group_colors          = c("Rhizosphere_vs_Roots" = "#56B4E9",
                            "Bulk soil_vs_Roots"   = "#E69F00"),
  x_axis_title          = "Section",
  show_x_labels         = FALSE,
  stat                  = "wilcox.test")
betaturn

leyenda <- cowplot::get_legend(betadis + ggplot2::theme(legend.position = "bottom"))
paneles <- cowplot::plot_grid(
  betadis  + ggplot2::theme(legend.position = "none"),
  betaturn + ggplot2::theme(legend.position = "none"),
  labels = c("A", "B"), rel_widths = c(1.5, 1))
betafig <- cowplot::plot_grid(paneles, leyenda, ncol = 1, rel_heights = c(1, 0.08))
betafig

#ggsave("plots/betaturn.png", betafig, width = 8, height = 6, units = "in", dpi = 300)


#beta decay
decay<-beta_decay_plot(
  table    = table_meta,
  metadata = metadata_meta,
  lat_col  = "Latitude",
  lon_col  = "Longitude",
  distance = "jaccard"
)

#ggsave("plots/betadecay.png",  width = 6, height = 4, units = "in", dpi = 300)




#plot
ord_row <- cowplot::plot_grid(beta_fungi, beta_meta, labels = c("A", "B"),
                              nrow = 1, align = "h")
tab_row <- cowplot::plot_grid(ggplot2::autoplot(betaperm_fungi), ggplot2::autoplot(beta_perm_meta), labels = c("C", "D"),
                              nrow = 1)
cmp_legend <- cowplot::get_legend(betadis + ggplot2::theme(legend.position = "bottom"))
cmp_row <- cowplot::plot_grid(
  betadis  + ggplot2::labs(x = NULL) + ggplot2::theme(legend.position = "none"),
  betaturn + ggplot2::labs(x = NULL) + ggplot2::theme(legend.position = "none"),
  decay,
  labels = c("E", "F", "G"), nrow = 1, rel_widths = c(1.1, 0.8, 1.5))
beta_fig <- cowplot::plot_grid(ord_row, tab_row, cmp_row, cmp_legend, ncol = 1,
                               rel_heights = c(1, 0.35, 0.8, 0.06))

ggsave("plots/beta_fig.png", beta_fig, width = 16, height = 14.5,
       units = "in", dpi = 300, bg = "white")

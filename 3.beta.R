source("1.load_data.R")


beta_fungi <- beta_ord_plot(table = table_fung,
              metadata = metadata_fungi,
              group_col = "Source",
              shape_col = "Treatment",
              distance = "compositional",
              ordination = "PCA")

betaperm_fungi <- beta_test_table(table = table_fung,
                metadata= metadata_fungi,
                formula_str = "Source",
                method = "compositional",
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
                method = "jaccard",
                test = "permanova",
                permutations = 999)



beta<- cowplot::plot_grid(beta_fungi, betaperm_fungi,
                          beta_meta,  beta_perm_meta,
                          labels = c("A", "B", "C", "D"),
                          label_y = 1, ncol = 2, nrow = 2, align = "hv")

ggsave("plots/beta.png", width = 10, height = 7, units = "in", dpi = 300)


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
  index = "jaccard",
  group_col = "Source",
  shape_col = "Treatment", 
  save_table = FALSE)


groups_fung <- setdiff(unique(as.character(metadata_fungi$Source)), "Roots")
vs_roots_fung <- paste0(pmin(groups_fung, "Roots"), "_vs_", pmax(groups_fung, "Roots"))

betadis <- beta_dissimilarity_plot(table = table_fung,
                        metadata = metadata_fungi,
                        comparison_condition1 = c("Rhizosphere_vs_Roots", "Rhizosphere_vs_Rhizosphere", "Bulk soil_vs_Roots"),
                        group_colors = c("Rhizosphere_vs_Roots" = "#E69F00",
                                         "Rhizosphere_vs_Rhizosphere" = "#56B4E9",
                                         "Bulk soil_vs_Roots" = "#009E73"),
                        condition1_col = "Source",
                        condition2_col = "Treatment",
                        x_axis_title = "Samples",
                        show_x_labels = FALSE,
                        x_label_angle = 45,
                       # stat = "kruskal.test",
                        partition = "turnover",
                        family = "jaccard",
                        save_table = FALSE)

betaturn <- beta_turnover_plot(
  table                 = table_fungi,
  metadata              = metadata_fungi,
  show_x_labels = FALSE,
  comparison_condition1 = c("Rhizosphere_vs_Roots",
                            "Rhizosphere_vs_Rhizosphere",
                            "Bulk soil_vs_Roots"),
  comparison_condition2 = c("TC_vs_TC", "TD_vs_TD", "TED_vs_TED"),
  condition1.x          = "Source.x",
  condition1.y          = "Source.y",
  condition2.x          = "Treatment.x",
  condition2.y          = "Treatment.y",
  facet_by              = "Treatment.x",
  #color_facets_x        = "#5D478B",
  color_axis_x          = c("Rhizosphere_vs_Roots" = "#E69F00",
                            "Rhizosphere_vs_Rhizosphere" = "#56B4E9",
                            "Bulk soil_vs_Roots" = "#009E73"))
cowplot::plot_grid(betadis, betaturn)

ggsave("plots/betaturn.png", width = 15, height = 8, units = "in", dpi = 300)


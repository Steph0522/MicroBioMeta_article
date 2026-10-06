source("1.load_data.R")

metadata_fungi_compar <-  metadata_fungi %>%
  filter(Type_of_soil == "Rhizosphere" | Type_of_soil =="Roots")

table_fungi_compar <- table_fung[match(metadata_fungi_compar$SAMPLEID, colnames(table_fung))]

table_fung_compar <- merge_feature_taxonomy(table_fungi_compar, taxonomy_fungi)

rhizo_roots_colors <- c(Rhizosphere = "#56B4E9", Roots = "#009E73")

# aldex heatmap
set.seed(123)
a1 <- aldex_heatmap_plot(table = table_fung_compar,
                         metadata = metadata_fungi_compar,
                         p_adjust_method = "none", pval_threshold = 0.05,
                         # no effect-size filter: same taxa (p < 0.05) as the volcano
                         effect_threshold = 0,
                         group_col = "Source",
                         group_colors = rhizo_roots_colors)

# aldex volcano
set.seed(123)
a2 <- aldex_volcano_plot(table = table_fung_compar,
                         metadata = metadata_fungi_compar,
                         group_col = "Source",
                         p_adjust_method = "none",
                         type ="volcano",
                         label_size = 5,
                         filter_uncultured = TRUE,
                         col_sup = rhizo_roots_colors[["Rhizosphere"]],
                         col_inf = rhizo_roots_colors[["Roots"]])

# ancom_plot
a3 <- ancombc_plot(
  table          = table_fung_compar,
  metadata       = metadata_fungi_compar,
  group_col      = "Source",
  level          = "Genus",
  min_prevalence = 0.3,
  p_adjust_method = "none",
  bar_colors     = unname(rhizo_roots_colors[c("Roots", "Rhizosphere")])
)

# random forest
set.seed(123)
a4 <- random_forest_lollipop_plot(
  table = table_fung,
  metadata = metadata_fungi,
  top_n = 10,
  variable_to_predict = "Source")

# bubble plot
a5 <- ratios_bubble_plot(table = table_fung_compar,
                         metadata = metadata_fungi_compar,
                         group_col = "Source",
                         top_n = 20,
                         condition_A = "Rhizosphere",
                         condition_B = "Roots",
                         group_colors = rhizo_roots_colors)
#plot
top_row <- cowplot::plot_grid(a1, a2, a3, labels = c("A", "B", "C"),
                              nrow = 1, rel_widths = c(1.2, 1, 1))
bottom_row <- cowplot::plot_grid(a4, a5, labels = c("D", "E"),
                                 nrow = 1, rel_widths = c(1, 1.1))
diff_abund <- cowplot::plot_grid(top_row, bottom_row, ncol = 1,
                                 rel_heights = c(1, 1.15))

ggsave("plots/differential_abundance.png", diff_abund,
       width = 18, height = 11, units = "in", dpi = 300, bg = "white")

source("1.load_data.R")

#alpha

acm1 <- alpha_hill_corr_plot(table = table_fung )
acm2 <- alpha_hill_corr_plot(table = table_meta, facet_orientation ="vertical", title = NULL)


acms <- cowplot::plot_grid(acm1, acm2,labels = "AUTO", rel_widths = c(1.2,1), ncol = 1)

#ggsave("plots/alpha_cors.png", width = 12, height = 8, units = "in", dpi = 300)


# alpha hill
a1 <- alpha_hill_plot(table = table_fung,
                metadata = metadata_fungi,
                x_col = "Treatment",
                fill_col = "Treatment",
                facet_by = "Source",
                facet_orientation = "horizontal",
                panel_labels = c("D", "E", "F", "G", "H", "I", "J", "K", "L"),
                free_y = TRUE,
                stat = "kruskal.test",
                show_legend = FALSE,
                save_table = FALSE)

a2 <- alpha_hill_plot(table = table_meta,
                metadata = metadata_meta,
                x_col = "Polygon",
                fill_col = "Polygon",
                facet_orientation = "vertical",
                panel_labels = c("D", "E", "F"),
                free_y = TRUE,
                legend_title = "",
                stat = "kruskal.test",
               show_legend = FALSE,
                save_table = FALSE)

#alphas <- cowplot::plot_grid(a1, a2,labels = "AUTO", rel_widths = c(1.2,1))


alphas1 <- cowplot::plot_grid(acm1, a1, rel_heights = c(1,2), ncol = 1)


ggsave("plots/alphas1.png", width = 13, height = 13, units = "in", dpi = 300)

alphas2 <- cowplot::plot_grid(acm2, a2,nrow = 1, rel_widths = c(1.2,1))
ggsave("plots/alphas2.png", width = 12, height = 8, units = "in", dpi = 300)
ggsave("plots/alphas2.png", width = 8, height = 12, units = "in", dpi = 300)


#decay
metadata_fungi_decay <- metadata_fungi %>%
  filter(Source %in% c("Bulk soil", "Rhizosphere"))

ad <- alpha_decay_plot(table = table_fung,
                 metadata = metadata_fungi_decay,
                 cont_var = "pH",
                 group_col = "Source",
                 x_axis_title = "Soil pH")

venn <- venn_plot(table = table_fung,   
          metadata = metadata_fungi,   
          merge_by = "Source",   
          min_prevalence = 0 )


alpha_more<- cowplot::plot_grid(ad, venn, labels = c("", "D"), label_size = 15, label_y = 1, rel_widths = c(1.7,1))

ggsave("plots/alpha_more.png", width = 11, height = 6, units = "in", dpi = 300)

# Comparing packages
devtools::install("../MicroBioMeta", upgrade = FALSE)
# Install latest version: devtools::install("../MicroBioMeta")
library(MicroBioMeta)
suppressPackageStartupMessages({
  library(phyloseq)
  library(rbiom)
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
})

#parameters
N_REP    <- 10
TOP_N    <- 20                                        
PALETA_MBM    <- "colorb"  # paleta de alfa, beta y CCA de MicroBioMeta ("colorb", "safe", "earth", ...)
PALETA_BARRAS <- "safe"   # abundance_bar_plot solo tiene "default", "brewer" y "safe"
FREE_Y_ALFA <- FALSE    # FALSE: una escala por índice (como rbiom); TRUE: una escala por panel
ANALISIS <- c("barras", "alfa", "heatmap", "beta", "cca")   
STATS <- c("permanova")
PERMANOVA_RELATIVA <- FALSE   
MEDIDAS_PHY   <- c("Chao1", "Shannon", "Simpson")    
MEDIDAS_RBIOM <- c("chao1", "shannon", "simpson")
dir_figuras <- file.path(getwd(), "figuras_comparacion")
dir.create(dir_figuras, showWarnings = FALSE, recursive = TRUE)

#load data
bacteria_table <- read.delim(
  system.file("extdata", "table_bacteria.txt", package = "MicroBioMeta"),
  row.names = 1, check.names = FALSE
)

taxonomy_bacteria <- read.delim(
  system.file("extdata", "taxonomy_bacteria.txt", package = "MicroBioMeta"),
  check.names = FALSE
) %>%
  rename(taxonomy = Taxon) %>%
  dplyr::select(-Confidence) %>%
  column_to_rownames(var = "Feature.ID")

bacteria_metadata <- read.delim(
  system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta"),
  check.names = FALSE
) %>%
  filter(Month == "2")

samples_bac <- bacteria_metadata$SAMPLEID[
  bacteria_metadata$SAMPLEID %in% colnames(bacteria_table)]

table_bacteria <- bacteria_table[, samples_bac]

metadata_bacteria <- bacteria_metadata[
  bacteria_metadata$SAMPLEID %in% samples_bac, ]
metadata_bacteria$Treatment <- as.factor(metadata_bacteria$Treatment)

# data MBM
table_bac <- merge_feature_taxonomy(table = table_bacteria,
                                    taxonomy = taxonomy_bacteria)


env_table_bac <- metadata_bacteria %>%
  dplyr::select(SAMPLEID, pH:Arbus_per) %>%
  remove_rownames() %>%
  column_to_rownames(var = "SAMPLEID")

#phyloseq object

otu_bac  <- otu_table(as.matrix(table_bacteria), taxa_are_rows = TRUE)
tax_bac  <- tax_table(as.matrix(taxonomy_bacteria))
meta_bac <- sample_data(data.frame(metadata_bacteria,
                                   row.names = metadata_bacteria$SAMPLEID))

ps_bacteria <- phyloseq(otu_bac, tax_bac, meta_bac)

tax_df <- as.data.frame(as(tax_table(ps_bacteria), "matrix"))
nombre_columna <- colnames(tax_df)[1]

tax_limpia <- tax_df %>%
  separate(!!sym(nombre_columna),
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";\\s*", fill = "right") %>%
  mutate(across(everything(), ~ gsub("^[a-z]__", "", .)))
tax_limpia[tax_limpia == ""] <- NA
rownames(tax_limpia) <- rownames(tax_df)
tax_table(ps_bacteria) <- tax_table(as.matrix(tax_limpia))

# rbiom
data <- as_rbiom(ps_bacteria)

# draw

pintar <- function(x) {
  if (is.function(x)) return(x())
  if (is.data.frame(x) && requireNamespace("gridExtra", quietly = TRUE)) {
    num <- vapply(x, is.numeric, logical(1))
    x[num] <- lapply(x[num], signif, digits = 4)
    grid::grid.newpage()
    gridExtra::grid.table(x)
    return(invisible(NULL))
  }
  print(x)
}

dibujar <- function(x) {
  f <- tempfile(fileext = ".png")
  png(f, width = 1600, height = 1200, res = 200)
  on.exit(dev.off())
  pintar(x)
  invisible(NULL)
}

dibujar_guardar <- function(x, nombre) {
  archivo <- file.path(dir_figuras, paste0(nombre, ".png"))
  png(filename = archivo, width = 1600, height = 1200, res = 200)
  pintar(x)
  dev.off()
  message("Guardado: ", normalizePath(archivo))
  invisible(NULL)
}

# bars
barras_mbm <- function(dib = dibujar) {
  p <- abundance_bar_plot(table = table_bac,
                          metadata = metadata_bacteria,
                          taxonomy_db = "silva",
                          level = "phylum",
                          top_n = 15,
                          x_col = "Type_of_soil",
                          facet_by = "Treatment",
                          add_remained = TRUE,
                          label = "Phylum",
                          x_label_angle = 45,
                          palette = PALETA_BARRAS,
                          taxa_order = "descending",
                          save_table = FALSE)
  dib(p)
}

barras_phy <- function(dib = dibujar) {
  ps_clas <- phyloseq::subset_taxa(ps_bacteria, !is.na(Phylum))
  ps_filo <- phyloseq::tax_glom(ps_clas, taxrank = "Phylum", NArm = FALSE)
  ps_rel <- phyloseq::transform_sample_counts(
  ps_filo, function(x) { total <- sum(x, na.rm = TRUE); if (total > 0) x / total else x })
  meta <- data.frame(sample_data(ps_rel), check.names = FALSE)
  meta$Group <- paste(meta$Type_of_soil, meta$Treatment, sep = "_")
  sample_data(ps_rel)$Group <- meta$Group
  ps_prom <- phyloseq::merge_samples(ps_rel, group = "Group", fun = mean)
  df_meta <- data.frame(Group = sample_names(ps_prom))
  df_meta <- tidyr::separate(df_meta, col = Group,
                             into = c("Type_of_soil", "Treatment"),
                             sep = "_", remove = FALSE)
  rownames(df_meta) <- df_meta$Group
  sample_data(ps_prom) <- phyloseq::sample_data(df_meta)
  ps_prom <- phyloseq::transform_sample_counts(ps_prom, function(x) x / sum(x))

  tax_mat   <- as(phyloseq::tax_table(ps_prom), "matrix")
  asignados <- rownames(tax_mat)[!is.na(tax_mat[, "Phylum"])]
  top15 <- names(sort(phyloseq::taxa_sums(ps_prom)[asignados], decreasing = TRUE))[1:15]
  otros <- !(rownames(tax_mat) %in% top15)
  tax_mat[otros, c("Domain", "Phylum")] <- "Otros"
  phyloseq::tax_table(ps_prom) <- phyloseq::tax_table(tax_mat)

  ps_final <- phyloseq::tax_glom(ps_prom, taxrank = "Phylum", NArm = FALSE)
  p <- phyloseq::plot_bar(ps_final, x = "Type_of_soil", fill = "Phylum",
                          facet_grid = ~Treatment) +
    ggplot2::geom_bar(ggplot2::aes(color = Phylum, fill = Phylum),
                      stat = "identity", position = "stack")
  dib(p)
}

# rbiom: group
barras_rbiom <- function(dib = dibujar) {
  ps_clas <- phyloseq::subset_taxa(ps_bacteria, !is.na(Phylum))

  ps_rel <- phyloseq::transform_sample_counts(
    ps_clas, function(x) { total <- sum(x, na.rm = TRUE); if (total > 0) x / total else x })

  meta <- data.frame(sample_data(ps_rel), check.names = FALSE)
  meta$Group <- paste(meta$Type_of_soil, meta$Treatment, sep = "_")
  sample_data(ps_rel)$Group <- meta$Group

  ps_prom <- phyloseq::merge_samples(ps_rel, group = "Group", fun = mean)

  df_meta <- data.frame(Group = sample_names(ps_prom))
  df_meta <- tidyr::separate(df_meta, col = Group,
                             into = c("Type_of_soil", "Treatment"),
                             sep = "_", remove = FALSE)
  rownames(df_meta) <- df_meta$Group
  sample_data(ps_prom) <- phyloseq::sample_data(df_meta)

  d  <- as_rbiom(ps_prom)
  d3 <- biom_relativize(d)
  p  <- taxa_stacked(d3, rank = "Phylum", taxa = 15,
                     label.by = "Type_of_soil", facet.by = "Treatment")
  dib(p)
}

# alfa

alfa_mbm <- function(dib = dibujar) {
  p <- alpha_diversity_plot(table = table_bac,
                            metadata = metadata_bacteria,
                            x_col = "Treatment",
                            fill_col = "Treatment",
                            facet_by = "Type_of_soil",
                            facet_orientation = "vertical",
                            save_table = FALSE,
                            palette = PALETA_MBM,
                            legend_title = "Treatment",
                            legend_position = "right",
                            aspect_ratio = 1,
                            free_y = FREE_Y_ALFA,
                            panel_labels = rep("", 12))   # sin letras internas (las pone la figura S2)
  dib(p)
}

alfa_phy <- function(dib = dibujar) {
  p <- phyloseq::plot_richness(ps_bacteria,
                               x = "Treatment",
                               color = "Treatment",
                               measures = MEDIDAS_PHY) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = Treatment), alpha = 0.6) +
    ggplot2::facet_grid(variable ~ Type_of_soil, scales = "free_y")
  dib(p)
}

alfa_rbiom <- function(dib = dibujar) {
  p <- rbiom::adiv_boxplot(data,
                           x = "Treatment",
                           facet.by = "Type_of_soil",
                           adiv = MEDIDAS_RBIOM)
  dib(p)
}

#heatmap

heatmap_mbm <- function(dib = dibujar) {
  dib(function() {
    suppressWarnings(
      print(abundance_heatmap_plot(table = table_bac,
                                   metadata = metadata_bacteria,
                                   top_n = TOP_N,
                                   exclude_unclassified = FALSE,
                                   cluster = TRUE,
                                   show_column_names = FALSE,
                                   condition1 = "Treatment",
                                   condition2 = "Type_of_soil",
                                   feature_prefix = "ASV",
                                   draw = FALSE))   # celdas que llenan el panel
    )
  })
}

heatmap_phy <- function(dib = dibujar) {
  ps_nz  <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_bacteria) > 0, ps_bacteria)
  ps_rel <- phyloseq::transform_sample_counts(ps_nz, function(x) x / sum(x) * 100)

  top_asv <- names(sort(phyloseq::taxa_sums(ps_rel), decreasing = TRUE))[
    seq_len(min(TOP_N, phyloseq::ntaxa(ps_rel)))]
  ps_top <- phyloseq::prune_taxa(top_asv, ps_rel)

  p <- phyloseq::plot_heatmap(ps_top,
                              method       = "NMDS",
                              sample.order = phyloseq::sample_names(ps_top),
                              taxa.order   = NULL,
                              trans        = scales::identity_trans())
  dib(function() print(p))
}

heatmap_rbiom <- function(dib = dibujar) {
  p <- rbiom::taxa_heatmap(data,
                           rank   = 7,      # comprobar que es el nivel ASV
                           taxa   = TOP_N,
                           tracks = c("Treatment", "Type_of_soil"),
                           grid   = blues9)
  dib(function() print(p))
}

# beta
beta_mbm <- function(dib = dibujar) {
  p <- beta_ord_plot(table = table_bac,
                     metadata = metadata_bacteria,
                     distance = "bray",
                     ordination = "PCoA",
                     group_col = "Treatment",
                     shape_col = "Type_of_soil",
                     palette = PALETA_MBM)
  dib(p)
}

beta_phy <- function(dib = dibujar) {
  ord <- phyloseq::ordinate(ps_bacteria, method = "PCoA", distance = "bray")
  p <- phyloseq::plot_ordination(ps_bacteria, ord,
                                 type  = "samples",
                                 color = "Treatment",
                                 shape = "Type_of_soil",
                                 title = "PCoA")
  dib(p)
}

# rbiom elipses
beta_rbiom_completo <- function(dib = dibujar) {
  p <- rbiom::bdiv_ord_plot(data,
                            bdiv = "bray",
                            ord = "PCoA",
                            layers = "petm",
                            stat.by = "Treatment",
                            facet.by = NULL,
                            colors = TRUE,
                            shapes = TRUE,
                            tree = NULL,
                            test = "adonis2",
                            seed = 0,
                            permutations = 999)
  dib(p)
}

beta_rbiom_min <- function(dib = dibujar) {
  p <- rbiom::bdiv_ord_plot(data,
                            bdiv = "Bray-Curtis",
                            ord = "PCoA",
                            layers = "p",
                            stat.by = "Treatment",
                            test = "none")
  dib(p)
}

# cca
cca_mbm <- function(dib = dibujar) {
  p <- cca_rda_biplot(table = table_bac,
                      env_data = env_table_bac,
                      metadata = metadata_bacteria,
                      env_vars = c("pH", "TN", "WHC", "EC", "Clay"),
                      analysis = "CCA",
                      show_all_env_vectors = TRUE,
                      group_col = "Type_of_soil",
                      legend_title = "Type_of_soil",
                      scale_arrows = 3)
  dib(p)
}

cca_phy <- function(dib = dibujar) {
  cca_resultado <- phyloseq::ordinate(ps_bacteria, method = "CCA",
                                      formula = ~ pH + TN + WHC + EC + Clay)

  p <- phyloseq::plot_ordination(ps_bacteria, cca_resultado, type = "samples",
                                 color = "Type_of_soil", shape = "Treatment") +
    ggplot2::theme_minimal() +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(title = "CCA")
  flechas <- as.data.frame(cca_resultado$CCA$biplot)
  flechas$Variable <- rownames(flechas)

  p <- p +
    ggplot2::geom_segment(data = flechas,
                          ggplot2::aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
                          arrow = grid::arrow(length = grid::unit(0.2, "cm")),
                          color = "blue", inherit.aes = FALSE) +
    ggplot2::geom_text(data = flechas,
                       ggplot2::aes(x = CCA1, y = CCA2, label = Variable),
                       color = "blue", vjust = -0.5, inherit.aes = FALSE)
  dib(p)
}

#permanova
permanova_mbm <- function(dib = dibujar) {
  set.seed(123)
  p <- beta_test_table(table = table_bac,
                       metadata = metadata_bacteria,
                       formula_str = "Treatment*Type_of_soil",
                       distance = "compositional",
                       test = "permanova",
                       permutations = 999)
  dib(p)
  invisible(p)
}

# rbiom: one  factor (Treatment)
permanova_rbiom_distmat <- function(dib = dibujar) {
  set.seed(123)
  dm     <- rbiom::bdiv_distmat(data, bdiv = "euclidean")
  grupos <- pull(data, "Treatment")
  p <- rbiom::distmat_stats(dm, groups = grupos, test = "adonis2")
  dib(p)
  invisible(p)
}

# register
TAREAS <- list(
  barras  = list(MBM = barras_mbm, phyloseq = barras_phy, rbiom = barras_rbiom),
  alfa    = list(MBM = alfa_mbm,   phyloseq = alfa_phy,   rbiom = alfa_rbiom),
  heatmap = list(MBM = heatmap_mbm, phyloseq = heatmap_phy, rbiom = heatmap_rbiom),
  beta    = list(MBM = beta_mbm,   phyloseq = beta_phy,
                 rbiom_completo = beta_rbiom_completo, rbiom_min = beta_rbiom_min),
  cca     = list(MBM = cca_mbm,    phyloseq = cca_phy)
)

TAREAS_STATS <- list(
  permanova = list(MBM = permanova_mbm,
                   rbiom_distmat = permanova_rbiom_distmat)
)

if ("cca" %in% ANALISIS && !exists("env_table_bac")) {
  message("env_table_bac no existe en la sesión: se omite el CCA.")
  ANALISIS <- setdiff(ANALISIS, "cca")
}

# execution and measure
medir <- function(f, n = N_REP, antes = function() NULL) {
  antes(); f()                                   # calentamiento
  t <- replicate(n, {
    antes()
    unname(system.time(f())["elapsed"])
  })
  c(mediana = median(t), min = min(t), max = max(t))
}

# clear caché
ANTES_RBIOM <- function() NULL
for (a in ANALISIS) {
  for (m in names(TAREAS[[a]])) {
    TAREAS[[a]][[m]](function(x) dibujar_guardar(x, paste0(a, "_", m)))
  }
}
message("Figuras en ", dir_figuras, ". Revísalas antes de dar por buenos los tiempos.")

# measure
res_lista <- list()
for (a in ANALISIS) {
  for (m in names(TAREAS[[a]])) {
    antes <- if (grepl("^rbiom", m)) ANTES_RBIOM else function() NULL
    r <- medir(TAREAS[[a]][[m]], antes = antes)
    res_lista[[paste(a, m, sep = "|")]] <- data.frame(
      analisis = a, metodo = m,
      mediana = r[["mediana"]], min = r[["min"]], max = r[["max"]])
    message(sprintf("%-8s %-15s mediana = %.3f s", a, m, r[["mediana"]]))
  }
}

# analyses

resultados_stats <- list()
for (a in STATS) {
  for (m in names(TAREAS_STATS[[a]])) {
    resultados_stats[[paste(a, m, sep = "|")]] <-
      TAREAS_STATS[[a]][[m]](function(x) dibujar_guardar(x, paste0(a, "_", m)))
    antes <- if (grepl("^rbiom", m)) ANTES_RBIOM else function() NULL
    r <- medir(TAREAS_STATS[[a]][[m]], antes = antes)
    res_lista[[paste(a, m, sep = "|")]] <- data.frame(
      analisis = a, metodo = m,
      mediana = r[["mediana"]], min = r[["min"]], max = r[["max"]])
    message(sprintf("%-9s %-17s mediana = %.3f s", a, m, r[["mediana"]]))
  }
}
if (length(resultados_stats)) {
  message("\n--- Resultados de los tests (compara estadístico, R2 y p) ---")
  for (k in names(resultados_stats)) { message("\n", k); print(resultados_stats[[k]]) }
}

res <- do.call(rbind, res_lista)
rownames(res) <- NULL
res[, c("mediana", "min", "max")] <- round(res[, c("mediana", "min", "max")], 3)
print(res)
write.csv(res, "results_todo.csv", row.names = FALSE)

# check equivalents tasks

n_m <- phyloseq::nsamples(ps_bacteria)
capa_puntos <- function(p) {
  L <- ggplot2::ggplot_build(p)$data
  L[[which(vapply(L, nrow, integer(1)) == n_m)[1]]][, c("x", "y")]
}
comparar <- function(a, b) c(
  eje1 = isTRUE(all.equal(sort(abs(a$x)), sort(abs(b$x)), tolerance = 1e-6)),
  eje2 = isTRUE(all.equal(sort(abs(a$y)), sort(abs(b$y)), tolerance = 1e-6)))
if ("beta" %in% ANALISIS) {
  c_mbm <- capa_puntos(beta_mbm(function(x) x))
  c_phy <- capa_puntos(beta_phy(function(x) x))
  print(comparar(c_mbm, c_phy))                  # deben ser TRUE
}

if ("heatmap" %in% ANALISIS) {
  tmp <- tempfile(fileext = ".txt")
  invisible(suppressWarnings(abundance_heatmap_plot(
    table = table_bac, metadata = metadata_bacteria, top_n = TOP_N,
    exclude_unclassified = FALSE, feature_prefix = "ASV",
    save_table = TRUE, table_filename = tmp, draw = FALSE)))
  ids_mbm <- read.delim(tmp)$OTUID
  ps_nz   <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_bacteria) > 0, ps_bacteria)
  ps_rel  <- phyloseq::transform_sample_counts(ps_nz, function(x) x / sum(x) * 100)
  ids_phy <- names(sort(phyloseq::taxa_sums(ps_rel), decreasing = TRUE))[1:TOP_N]
  message("Heatmap, MBM y phyloseq: mismos ", TOP_N, " ASV = ", setequal(ids_mbm, ids_phy))
}

# preparing data
t_as_rbiom <- median(replicate(N_REP,
                               unname(system.time(as_rbiom(ps_bacteria))["elapsed"])))
message("as_rbiom(ps_bacteria): ", round(t_as_rbiom, 3), " s")

t_ps <- median(replicate(N_REP, unname(system.time({
  ps_tmp <- phyloseq(otu_table(as.matrix(table_bacteria), taxa_are_rows = TRUE),
                     tax_table(as.matrix(taxonomy_bacteria)),
                     sample_data(data.frame(metadata_bacteria,
                                            row.names = metadata_bacteria$SAMPLEID)))
  td <- as.data.frame(as(tax_table(ps_tmp), "matrix"))
  tl <- td %>%
    separate(!!sym(colnames(td)[1]),
             into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
             sep = ";\\s*", fill = "right") %>%
    mutate(across(everything(), ~ gsub("^[a-z]__", "", .)))
  tl[tl == ""] <- NA
  rownames(tl) <- rownames(td)
  tax_table(ps_tmp) <- tax_table(as.matrix(tl))
})["elapsed"])))
message("Construir ps_bacteria + limpiar taxonomía: ", round(t_ps, 3), " s")

# from raw data
#plots
BETA_RBIOM <- "rbiom_completo"  

COLUMNAS <- c(MBM = "MicroBioMeta", phyloseq = "phyloseq", rbiom = "rbiom")

FILAS_S1 <- list(
  "Relative abundance (bar plot)" = c(MBM = "barras_MBM", phyloseq = "barras_phyloseq", rbiom = "barras_rbiom"),
  "Relative abundance (heatmap)"  = c(MBM = "heatmap_MBM", phyloseq = "heatmap_phyloseq", rbiom = "heatmap_rbiom")
)
FILAS_S2 <- list(
  "Alpha diversity"           = c(MBM = "alfa_MBM", phyloseq = "alfa_phyloseq", rbiom = "alfa_rbiom"),
  "Beta diversity ordination" = c(MBM = "beta_MBM", phyloseq = "beta_phyloseq", rbiom = paste0("beta_", BETA_RBIOM)),
  "CCA"                       = c(MBM = "cca_MBM", phyloseq = "cca_phyloseq", rbiom = NA),
  "PERMANOVA"                 = c(MBM = "permanova_MBM", phyloseq = NA, rbiom = "permanova_rbiom_distmat")
)

texto <- function(x, ...) cowplot::ggdraw() + cowplot::draw_label(x, fontface = "bold", size = 14, ...)

componer <- function(filas, nombre, ancho = 15) {
  celdas <- list(texto(""))
  for (m in names(COLUMNAS)) celdas <- c(celdas, list(texto(COLUMNAS[[m]])))
  k <- 0
  for (fila in names(filas)) {
    celdas <- c(celdas, list(texto(fila, angle = 90)))
    for (m in names(COLUMNAS)) {
      archivo <- file.path(dir_figuras, paste0(filas[[fila]][[m]], ".png"))
      if (is.na(filas[[fila]][[m]]) || !file.exists(archivo)) {
        if (!is.na(filas[[fila]][[m]])) message("No existe ", archivo, ": se deja vacío")
        celdas <- c(celdas, list(texto("—", colour = "grey60")))
        next
      }
      k <- k + 1
      celdas <- c(celdas, list(
        cowplot::ggdraw() +
          cowplot::draw_image(archivo, x = 0.04, y = 0, width = 0.96, height = 0.95) +   
          cowplot::draw_label(LETTERS[k], x = 0.01, y = 0.99, hjust = 0, vjust = 1,
                              fontface = "bold", size = 18)))
    }
  }
  w_etiqueta <- 0.08; h_titulo <- 0.08
  fig <- cowplot::plot_grid(plotlist = celdas, ncol = length(COLUMNAS) + 1,
                            rel_widths = c(w_etiqueta, rep(1, length(COLUMNAS))),
                            rel_heights = c(h_titulo, rep(0.75, length(filas))))
  alto <- ancho / (w_etiqueta + length(COLUMNAS)) * (h_titulo + 0.75 * length(filas))
  for (ext in c("pdf", "png")) {
    archivo <- file.path(dir_figuras, paste0(nombre, ".", ext))
    ggsave(archivo, fig, width = ancho, height = alto, dpi = 300, bg = "white")
    message("Guardado: ", normalizePath(archivo))
  }
  invisible(fig)
}

componer(FILAS_S1, "Figure_S1_taxonomic_composition")
componer(FILAS_S2, "Figure_S2_diversity_ordination")

sessionInfo()


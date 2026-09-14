#loading data
library(qiime2R)
library(tidyverse)
library(MicroBioMeta)

#load and pre-processing data of metabarcoding
fungi_table <- read_qza("data/table_fungi.qza")$data %>%
  as.data.frame()

taxonomy_fungi <- read_qza("data/taxonomy_fungi.qza")$data %>% as.data.frame() %>%
  rename(taxonomy = Taxon) %>%
  dplyr::select(-Consensus) %>%
  column_to_rownames(var = "Feature.ID")

fungi_metadata <- read.delim("data/metadata_fungis.txt", check.names = FALSE)

samples_fun <- fungi_metadata$SAMPLEID[
  fungi_metadata$SAMPLEID %in% colnames(fungi_table)
]

table_fungi <- fungi_table[, samples_fun]

metadata_fungi <- fungi_metadata[
  fungi_metadata$SAMPLEID %in% samples_fun,
]

table_fung <- merge_feature_taxonomy(table = table_fungi, 
                                     taxonomy = taxonomy_fungi)



#table metagenomic

table_kraken <- read.delim("data/table_kraken.txt", row.names = 1, check.names = FALSE)

table_kraken[1:3, c(1:2, ncol(table_kraken))]

metadata_meta <- read.delim("data/metadata_metagenomic.txt", check.names = FALSE)

sample_cols <- setdiff(colnames(table_kraken), "taxonomy")

id_metagenome <- as.integer(gsub("kraken_pluspfp_|_report_bracken_species", "", sample_cols))

colnames(table_kraken)[match(sample_cols, colnames(table_kraken))] <-
  metadata_meta$SAMPLEID[match(id_metagenome, metadata_meta$id_metagenome)]

metadata_meta$Polygon <- factor(paste0("Pol", metadata_meta$Poligono), levels = paste0("Pol", 1:6))
metadata_meta$Site <- as.character(metadata_meta$Sitio)
table_meta <- table_kraken[, c(metadata_meta$SAMPLEID, "taxonomy")]



## --- Inputs ---------------------------------------------------------------
# Set the path to your clusters .txt file (tab-delimited, first 2 cols: name/type; rest: Entrez IDs)
clusters_path <- ""   # <-- adjust path if needed

## --- Libraries ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(annotables)        # provides grch38
library(clusterProfiler)

## --- 1) Get DEGs in ENTREZ ID space --------------------------------------
# Use your significance rule (padj < 0.05 here; change if needed)
ensembl_ids <- sub("\\.\\d+$", "", rownames(res))
deg_entrez_ids <- mapIds(org.Hs.eg.db,
                         keys      = ensembl_ids[!is.na(res$padj) & res$padj < 0.05],
                         column    = "ENTREZID",
                         keytype   = "ENSEMBL",
                         multiVals = "first") |> unname()
deg_entrez_ids <- unique(na.omit(deg_entrez_ids))

## --- 2) Protein-coding background (universe) in ENTREZ IDs ----------------
pc_universe_entrez <- unique(grch38$entrez[grch38$biotype == "protein_coding"])
pc_universe_entrez <- pc_universe_entrez[!is.na(pc_universe_entrez)]

## --- 3) Read clusters and build TERM2GENE/TERM2NAME -----------------------
clusters_raw <- read.delim(clusters_path, check.names = FALSE)  # keeps spaces in colnames
colnames(clusters_raw)[1:2] <- c("cluster_name", "cluster_type")

# Pivot all columns from the 3rd to the last into a single 'entrez' column
clusters_long <- clusters_raw %>%
  tidyr::pivot_longer(
    cols = 3:ncol(.),
    names_to  = "gene_col",
    values_to = "entrez"
  ) %>%
  dplyr::filter(!is.na(entrez) & entrez != "") %>%
  dplyr::mutate(entrez = as.character(entrez)) %>%
  dplyr::distinct(cluster_name, entrez)


TERM2GENE <- clusters_long |>
  select(term = cluster_name, gene = entrez)

TERM2NAME <- clusters_raw |>
  distinct(cluster_name) |>
  transmute(term = cluster_name, name = cluster_name)

## --- 4) Run enrichment against your custom clusters -----------------------
# IMPORTANT: 'gene' and 'universe' must be the same ID type (ENTREZ here)
enr_clusters <- enricher(gene       = deg_entrez_ids,
                         TERM2GENE  = TERM2GENE,
                         TERM2NAME  = TERM2NAME,
                         universe   = pc_universe_entrez,
                         pAdjustMethod = "BH",
                         qvalueCutoff  = 0.05)

# Inspect & save
head(as.data.frame(enr_clusters))
write.csv(as.data.frame(enr_clusters), "custom_cluster_enrichment.csv", row.names = FALSE)

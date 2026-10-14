# Huang HEK293 KO panel: substructure-match + GSEA per KO.
#
# All 7 KO glycomeRes share the same (group, toposeq) backbone — only the
# logFC differs. So we run step3-style substructure matching ONCE on the
# shared backbone, build the enzyme×group OR-collapsed gmt ONCE, then loop
# 7 KOs through step4-style fgsea with KO-specific logFC.
#
# Pattern adapted from:
#   src/scripts/step3_build_pk.R           (substructure + class rules + OR collapse)
#   src/scripts/step4_glycome_ORA.R        (fgsea + qvalue)
#   validation/toglycoenz/scripts/02_run_step3_step4.R (per-dataset wrapper)
# Run from the repository root.

suppressMessages({
  library(parallel); library(dplyr); library(tibble); library(igraph)
  library(fgsea); library(qvalue)
})
source("R/utils_formatting.R")

OUTDIR <- "output/hek293"
KOS <- c("MGAT1_KO", "MGAT2_KO", "MGAT4A_4B_KO", "MGAT5_KO",
         "T_KO_MANA1_A2_B1_", "B4GALNT3_4_KO", "SLC35C1_KO")

# ---------- 1. shared triplet DB ----------
triplet_enzdb  <- readRDS("output/triplet_enzdb.rds")
triplet_glygen <- readRDS("output/triplet_glygen.rds")

triplet_enzdb_modi <- triplet_enzdb %>%
  dplyr::select(Enzyme_gene, Reactant, Product, Enzyme_class, dbID = ENZdbID) %>%
  dplyr::distinct(Enzyme_gene, Reactant, Product, .keep_all = TRUE)
triplet_glygen_modi <- triplet_glygen %>%
  dplyr::select(Enzyme_gene, Reactant, Product, Enzyme_class, dbID = GlyGendbID) %>%
  dplyr::distinct(Enzyme_gene, Reactant, Product, .keep_all = TRUE)
triplet <- rbind(triplet_enzdb_modi, triplet_glygen_modi) %>%
  dplyr::distinct(Enzyme_gene, Reactant, Product, .keep_all = TRUE)
cat(sprintf("triplet rows: %d\n", nrow(triplet)))

# ---------- 2. shared (group, toposeq) backbone ----------
# Take any one KO file; all 7 share the same backbone (verified earlier).
gly_template <- readRDS(file.path(OUTDIR, "glycomeRes_MGAT1_KO.rds"))
gly_template <- na.omit(gly_template)

unique_toposeq <- unique(gly_template$toposeq)
unique_groups  <- unique(gly_template$group)
cat(sprintf("unique toposeq: %d  unique groups: %d\n",
            length(unique_toposeq), length(unique_groups)))

# ---------- 3. substructure matching (toposeq × triplet) ----------
cat("Building triplet igraphs (parallel) ...\n")
g_triplet_reactant <- mclapply(triplet$Reactant, function(x) {
  x %>% build_tree_toposeq() %>% tree_to_igraph() %>% canonicalize_igraph()
}, mc.cores = 20)
g_triplet_product <- mclapply(triplet$Product, function(x) {
  x %>% build_tree_toposeq() %>% tree_to_igraph() %>% canonicalize_igraph()
}, mc.cores = 20)

cat("Building glycome igraphs ...\n")
g_glycome <- lapply(unique_toposeq, function(x) {
  x %>% build_tree_toposeq() %>% tree_to_igraph() %>% canonicalize_igraph()
})

cat("Substructure matching ...\n")
t0 <- Sys.time()
mat_react   <- substructure_match_batch(g_triplet_reactant, g_glycome, mc.cores = 20)
mat_product <- substructure_match_batch(g_triplet_product,  g_glycome, mc.cores = 20)
cat(sprintf("elapsed: %.1f sec\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

rownames(mat_react)   <- rownames(mat_product) <- triplet$Enzyme_gene
colnames(mat_react)   <- colnames(mat_product) <- unique_toposeq

# ---------- 4. class-aware AND / AND-NOT logic (same as step3) ----------
res <- matrix(NA, nrow = nrow(mat_react), ncol = ncol(mat_react),
              dimnames = dimnames(mat_react))
is_genz_gt <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == "Glycosyltransferase"
is_glygen  <- grepl("GLYGEN", triplet$dbID)
row_idx <- is_genz_gt | is_glygen
res[row_idx, ] <- mat_react[row_idx, ] & mat_product[row_idx, ]
row_idx <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == ""
res[row_idx, ] <- mat_react[row_idx, ] & mat_product[row_idx, ]
row_idx <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == "Glycosidase"
res[row_idx, ] <- mat_product[row_idx, ] & !mat_react[row_idx, ]
row_idx <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == "N_linked"
res[row_idx, ] <- mat_react[row_idx, ] & mat_product[row_idx, ]

# Enzyme dedup: OR across rows with the same Enzyme_gene
enz_names <- unique(rownames(res))
res_summarized <- t(sapply(enz_names, function(e) {
  rows <- res[rownames(res) == e, , drop = FALSE]
  if (nrow(rows) == 1) return(rows[1, ])
  apply(rows, 2, function(x) {
    if (all(is.na(x))) return(NA)
    any(x, na.rm = TRUE)
  })
}))
rownames(res_summarized) <- enz_names

# Cleaning: drop all-FALSE / all-NA enzymes and full-matching enzymes
out_all_true <- rowSums(res_summarized, na.rm = TRUE) == ncol(res_summarized)
out_all_false <- rowSums(res_summarized, na.rm = TRUE) == 0
res_clean <- res_summarized[!(out_all_true | out_all_false), , drop = FALSE]
cat(sprintf("res_clean (toposeq level): %d enzymes × %d toposeq\n",
            nrow(res_clean), ncol(res_clean)))
saveRDS(res_clean, file.path(OUTDIR, "res_clean.rds"))

# ---------- 5. OR-collapse to group level (per spec §D2, step3:198-244) ----
res_clean_group <- matrix(
  NA, nrow = nrow(res_clean), ncol = length(unique_groups),
  dimnames = list(rownames(res_clean), unique_groups)
)
for (i in seq_along(unique_groups)) {
  g <- unique_groups[i]
  target_toposeq <- intersect(
    gly_template$toposeq[gly_template$group == g],
    colnames(res_clean)
  )
  if (length(target_toposeq) == 0) next
  if (length(target_toposeq) == 1) {
    res_clean_group[, g] <- res_clean[, target_toposeq]
  } else {
    sub <- res_clean[, target_toposeq, drop = FALSE]
    res_clean_group[, g] <- apply(sub, 1, function(x) {
      if (all(is.na(x))) return(NA)
      as.integer(any(x, na.rm = TRUE))
    })
  }
}

# Cleaning at group level
col_has_any <- colSums(res_clean_group, na.rm = TRUE) > 0
res_clean_group <- res_clean_group[, col_has_any, drop = FALSE]
out_all_true  <- rowSums(res_clean_group, na.rm = TRUE) == ncol(res_clean_group)
out_all_false <- rowSums(res_clean_group, na.rm = TRUE) == 0
res_clean_group <- res_clean_group[!(out_all_true | out_all_false), , drop = FALSE]
cat(sprintf("res_clean_group: %d enzymes × %d groups\n",
            nrow(res_clean_group), ncol(res_clean_group)))
saveRDS(res_clean_group, file.path(OUTDIR, "res_clean_group.rds"))

# ---------- 6. enzyme pattern grouping (REVIGO-like) ----------
mat_numeric <- res_clean_group * 1
mat_numeric[is.na(mat_numeric)] <- -1
row_profiles <- apply(mat_numeric, 1, paste, collapse = "")
enzyme_pattern_groups <- split(names(row_profiles), row_profiles)

extract_label <- function(members) {
  if (length(members) == 1) return(members)
  manual_map <- list(
    "MGAT4/5" = c("MGAT4A", "MGAT4B", "MGAT4C", "MGAT4D", "MGAT5"),
    "ST3GAL/ST6GAL/CMAS" = c("CMAS", "ST3GAL1", "ST3GAL3", "ST3GAL4",
                              "ST3GAL6", "ST6GAL1", "ST6GAL2"),
    "B4GALT/NEU" = c("B3GALT2", "B4GALT1", "B4GALT2", "B4GALT3",
                      "B4GALT4", "B4GALT5", "NEU1", "NEU2", "NEU3", "NEU4")
  )
  ms <- sort(members)
  for (lab in names(manual_map)) if (identical(ms, sort(manual_map[[lab]]))) return(lab)
  prefixes <- unique(gsub("[0-9]+[A-Z]?$", "", members))
  if (length(prefixes) == 1) {
    return(paste0(prefixes,
                  paste(gsub(prefixes, "", members), collapse = "/")))
  }
  paste(members, collapse = "/")
}

group_labels <- sapply(enzyme_pattern_groups, extract_label)
names(enzyme_pattern_groups) <- group_labels

# Collapse to one row per pattern group
res_grouped <- t(sapply(enzyme_pattern_groups, function(members) {
  res_clean_group[members[1], ]
}))
rownames(res_grouped) <- names(enzyme_pattern_groups)
gmt_enzyme_grouped <- apply(res_grouped, 1, function(row)
  colnames(res_grouped)[which(row & !is.na(row))])
gmt_enzyme_grouped <- gmt_enzyme_grouped[lengths(gmt_enzyme_grouped) > 0]
cat(sprintf("enzyme pattern groups: %d (was %d single enzymes)\n",
            length(gmt_enzyme_grouped), nrow(res_clean_group)))
saveRDS(gmt_enzyme_grouped, file.path(OUTDIR, "gmt_enzyme_grouped.rds"))
saveRDS(enzyme_pattern_groups, file.path(OUTDIR, "enzyme_pattern_groups.rds"))

# ---------- 7. GSEA per KO ----------
run_gsea_for_ko <- function(ko) {
  gly <- readRDS(file.path(OUTDIR, sprintf("glycomeRes_%s.rds", ko)))
  gly <- na.omit(gly)
  # one logFC per group (already constant within group; take first)
  lfc <- gly %>%
    dplyr::distinct(group, .keep_all = TRUE) %>%
    dplyr::select(group, logFC) %>%
    tibble::deframe()
  lfc <- lfc[is.finite(lfc)]
  lfc <- sort(lfc)
  set.seed(42)
  res <- fgsea(pathways = gmt_enzyme_grouped,
               stats    = lfc,
               minSize  = 1,
               maxSize  = Inf,
               eps      = 0,
               nPermSimple = 50000)
  res <- as.data.frame(res)
  res$qvalue <- tryCatch(
    qvalue(res$pval, pi0.method = "bootstrap")$qvalues,
    error = function(e) res$padj
  )
  res <- res %>% dplyr::arrange(pval)
  res$ko <- ko
  res
}

resGSEA_by_ko <- list()
for (ko in KOS) {
  cat(sprintf("\n--- GSEA: %s ---\n", ko))
  r <- run_gsea_for_ko(ko)
  resGSEA_by_ko[[ko]] <- r
  cat(sprintf("  pathways: %d; q<0.05: %d; q<0.10: %d; pval<0.05: %d\n",
              nrow(r),
              sum(r$qvalue < 0.05, na.rm = TRUE),
              sum(r$qvalue < 0.10, na.rm = TRUE),
              sum(r$pval   < 0.05, na.rm = TRUE)))
  cat("  top 5:\n")
  print(r %>% head(5) %>%
          dplyr::select(pathway, NES, pval, qvalue, size))
}

saveRDS(resGSEA_by_ko, file.path(OUTDIR, "resGSEA_by_ko.rds"))
cat("\nDone. Saved output/hek293/resGSEA_by_ko.rds\n")

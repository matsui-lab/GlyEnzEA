# =============================================================================
# step3b_glycodecipher_mapping.R
# Full GlycoDecipher library mapping against enzyme reaction knowledge
#
# Builds the same enzyme x glycan binary matrix as step3_build_pk.R,
# but for the entire GlycoDecipher library (~10K structures) instead of
# only the Kizuka MCP 2019 validation set (97 structures).
# =============================================================================

library(dplyr)
library(readr)
library(igraph)
library(parallel)

# --- Working directory handling ---
# When run from scripts/ folder, move up to project root
# Run from the repository root.
cat("Working directory:", getwd(), "\n\n")

source("R/utils_formatting.R")


# =============================================================================
# Section 1: Data loading
# =============================================================================

cat("=== Section 1: Loading data ===\n")

glycan_lib <- read_csv(
  "data/GlycanLibrary/glycanDatabase_for_GlycoDecipher.csv",
  show_col_types = FALSE
)
cat("GlycoDecipher library loaded:", nrow(glycan_lib), "entries\n")

triplet <- readRDS("output/triplet.rds")
cat("Triplet data loaded:", nrow(triplet), "enzyme-reaction pairs\n")
cat("  Unique enzymes:", length(unique(triplet$Enzyme_gene)), "\n")
cat("  Enzyme classes:", paste(unique(triplet$Enzyme_class), collapse = ", "), "\n\n")


# =============================================================================
# Section 2: Convert GlycoDecipher WURCS -> topoSeq
# =============================================================================

cat("=== Section 2: WURCS -> topoSeq conversion ===\n")

n_cores <- max(1, detectCores() - 1)
cat("Using", n_cores, "cores for parallel processing\n")

lib_wurcs <- glycan_lib$StructureInformation
names(lib_wurcs) <- glycan_lib$GlycanID

cat("Converting", length(lib_wurcs), "WURCS strings to topoSeq...\n")

lib_toposeq <- mclapply(lib_wurcs, function(x) {
  tryCatch({
    x %>%
      build_tree_wurcs() %>%
      modify_tree_wurcs() %>%
      tree_to_igraph() %>%
      canonicalize_igraph() %>%
      igraph_to_toposeq()
  }, error = function(e) NA_character_)
}, mc.cores = n_cores)

lib_toposeq <- unlist(lib_toposeq)
lib_toposeq_valid <- lib_toposeq[!is.na(lib_toposeq)]

cat("topoSeq conversion complete:\n")
cat("  Total entries:", length(lib_toposeq), "\n")
cat("  Valid topoSeqs:", length(lib_toposeq_valid), "\n")
cat("  Unique topoSeqs:", length(unique(lib_toposeq_valid)), "\n\n")

# Work with unique topoSeqs to avoid redundant matching
unique_toposeqs <- unique(lib_toposeq_valid)
cat("Working with", length(unique_toposeqs), "unique topoSeqs\n\n")


# =============================================================================
# Section 3: Convert topoSeqs to igraph objects
# =============================================================================

cat("=== Section 3: Building igraph objects ===\n")

# Convert library glycans
cat("Converting", length(unique_toposeqs), "library topoSeqs to igraph...\n")
g_library <- mclapply(unique_toposeqs, function(x) {
  tryCatch({
    x %>% build_tree_toposeq() %>% tree_to_igraph() %>% canonicalize_igraph()
  }, error = function(e) NULL)
}, mc.cores = n_cores)

# Remove failed conversions
valid_idx <- !sapply(g_library, is.null)
g_library <- g_library[valid_idx]
unique_toposeqs <- unique_toposeqs[valid_idx]
cat("Valid igraph objects:", length(g_library), "\n")

# Convert triplet reactants and products
cat("Converting triplet reactants and products to igraph...\n")
g_triplet_reactant <- mclapply(triplet$Reactant, function(x) {
  x %>% build_tree_toposeq() %>% tree_to_igraph() %>% canonicalize_igraph()
}, mc.cores = n_cores)

g_triplet_product <- mclapply(triplet$Product, function(x) {
  x %>% build_tree_toposeq() %>% tree_to_igraph() %>% canonicalize_igraph()
}, mc.cores = n_cores)

cat("Triplet igraph conversion complete\n\n")


# =============================================================================
# Section 4: Substructure matching
# =============================================================================

cat("=== Section 4: Substructure matching ===\n")
cat("Matrix dimensions:", nrow(triplet), "triplets x", length(g_library), "glycans\n")

# Batch substructure matching with progress reporting
# Process in chunks to report progress
chunk_size <- 1000
n_glycans <- length(g_library)
n_chunks <- ceiling(n_glycans / chunk_size)

cat("Matching reactants against library (", n_chunks, " chunks)...\n")
t_start <- Sys.time()

mat_triplet_reactant <- substructure_match_batch(
  g_triplet_reactant, g_library, mc.cores = n_cores
)

t_reactant <- Sys.time()
cat("  Reactant matching done in",
    round(difftime(t_reactant, t_start, units = "mins"), 1), "minutes\n")

cat("Matching products against library...\n")
mat_triplet_product <- substructure_match_batch(
  g_triplet_product, g_library, mc.cores = n_cores
)

t_product <- Sys.time()
cat("  Product matching done in",
    round(difftime(t_product, t_reactant, units = "mins"), 1), "minutes\n")

# Label rows and columns
rownames(mat_triplet_reactant) <- triplet$Enzyme_gene
colnames(mat_triplet_reactant) <- unique_toposeqs
rownames(mat_triplet_product) <- triplet$Enzyme_gene
colnames(mat_triplet_product) <- unique_toposeqs


# =============================================================================
# Section 4b: Apply GT/GH rules
# =============================================================================

cat("\n=== Section 4b: Applying enzyme class rules ===\n")

res <- matrix(
  NA,
  nrow = nrow(mat_triplet_reactant),
  ncol = ncol(mat_triplet_reactant),
  dimnames = dimnames(mat_triplet_reactant)
)

# GT + GlyGen: reactant AND product
is_genz_gt <- grepl("GENZ", triplet$dbID) &
  triplet$Enzyme_class == "Glycosyltransferase"
is_glygen <- grepl("GLYGEN", triplet$dbID)
row_idx <- is_genz_gt | is_glygen
res[row_idx, ] <- mat_triplet_reactant[row_idx, ] & mat_triplet_product[row_idx, ]

# Empty class: same as GT
row_idx <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == ""
res[row_idx, ] <- mat_triplet_reactant[row_idx, ] & mat_triplet_product[row_idx, ]

# GH: product AND NOT reactant
row_idx <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == "Glycosidase"
res[row_idx, ] <- mat_triplet_product[row_idx, ] & !mat_triplet_reactant[row_idx, ]

# N_linked: same as GT
row_idx <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == "N_linked"
res[row_idx, ] <- mat_triplet_reactant[row_idx, ] & mat_triplet_product[row_idx, ]

cat("Rule application complete\n")
cat("  Non-NA entries:", sum(!is.na(res)), "of", length(res), "\n")
cat("  TRUE entries:", sum(res, na.rm = TRUE), "\n\n")


# =============================================================================
# Section 5: Aggregate by enzyme and clean
# =============================================================================

cat("=== Section 5: Aggregation and cleaning ===\n")

# Aggregate duplicate enzyme rows (OR logic)
enz <- unique(rownames(res))
for (i in seq_along(enz)) {
  e <- enz[i]
  idx <- rownames(res) == e
  if (sum(idx) <= 1) next
  m <- res[idx, ]
  col_has_true <- apply(m, 2, any)
  res[e, ] <- col_has_true
}
dup <- duplicated(rownames(res))
res_summarize <- res[!dup, ]

cat("After enzyme aggregation:", nrow(res_summarize), "enzymes x",
    ncol(res_summarize), "glycans\n")

# Remove enzymes that map to ALL glycans (uninformative)
out_row <- rowSums(res_summarize, na.rm = TRUE) == ncol(res_summarize)
res_clean <- res_summarize[!out_row, , drop = FALSE]
cat("  Removed", sum(out_row), "enzymes mapping to all glycans\n")

# Remove enzymes that map to NO glycans
out_row <- rowSums(res_clean, na.rm = TRUE) == 0
res_clean <- res_clean[!out_row, , drop = FALSE]
cat("  Removed", sum(out_row), "enzymes mapping to no glycans\n")

cat("Final matrix:", nrow(res_clean), "enzymes x", ncol(res_clean), "glycans\n\n")


# =============================================================================
# Section 6: Summary statistics
# =============================================================================

cat("=== Section 6: Summary statistics ===\n")
cat("Total library glycans (valid topoSeq):", length(lib_toposeq_valid), "\n")
cat("Unique topoSeqs in library:", length(unique_toposeqs), "\n")
cat("Enzymes in final matrix:", nrow(res_clean), "\n")
cat("Glycans with at least one enzyme:", sum(colSums(res_clean, na.rm = TRUE) > 0), "\n")

mappings_per_enzyme <- rowSums(res_clean, na.rm = TRUE)
cat("Mappings per enzyme: min =", min(mappings_per_enzyme),
    ", median =", median(mappings_per_enzyme),
    ", max =", max(mappings_per_enzyme), "\n")

mappings_per_glycan <- colSums(res_clean, na.rm = TRUE)
cat("Enzymes per glycan: min =", min(mappings_per_glycan),
    ", median =", median(mappings_per_glycan),
    ", max =", max(mappings_per_glycan), "\n\n")

# Save results
saveRDS(res_clean, "output/res_clean_glycodecipher.rds")
cat("Saved: output/res_clean_glycodecipher.rds\n")

# Also save the topoSeq mapping for reference
lib_toposeq_df <- data.frame(
  GlycanID = names(lib_toposeq),
  topoSeq = unname(lib_toposeq),
  stringsAsFactors = FALSE
)
saveRDS(lib_toposeq_df, "output/glycodecipher_toposeq_map.rds")
cat("Saved: output/glycodecipher_toposeq_map.rds\n\n")


# =============================================================================
# Section 7: Comparison with Kizuka subset
# =============================================================================

cat("=== Section 7: Comparison with Kizuka subset ===\n")

if (file.exists("output/res_clean.rds")) {
  res_kizuka <- readRDS("output/res_clean.rds")
  cat("Kizuka res_clean:", nrow(res_kizuka), "enzymes x",
      ncol(res_kizuka), "glycans\n")
  cat("GlycoDecipher res_clean:", nrow(res_clean), "enzymes x",
      ncol(res_clean), "glycans\n\n")

  # Check enzyme overlap
  enz_kizuka <- rownames(res_kizuka)
  enz_gd <- rownames(res_clean)
  cat("Enzymes in Kizuka:", length(enz_kizuka), "\n")
  cat("Enzymes in GlycoDecipher:", length(enz_gd), "\n")
  cat("Overlap:", length(intersect(enz_kizuka, enz_gd)), "\n")
  cat("Kizuka-only:", length(setdiff(enz_kizuka, enz_gd)), "\n")
  cat("GlycoDecipher-only:", length(setdiff(enz_gd, enz_kizuka)), "\n\n")

  # Check glycan (topoSeq) overlap
  ts_kizuka <- colnames(res_kizuka)
  ts_gd <- colnames(res_clean)
  cat("topoSeqs in Kizuka:", length(ts_kizuka), "\n")
  cat("topoSeqs in GlycoDecipher:", length(ts_gd), "\n")
  cat("Overlap:", length(intersect(ts_kizuka, ts_gd)), "\n")
  cat("Kizuka subset coverage:",
      round(length(intersect(ts_kizuka, ts_gd)) / length(ts_kizuka) * 100, 1),
      "%\n\n")

  # For overlapping enzymes and glycans, check consistency
  common_enz <- intersect(enz_kizuka, enz_gd)
  common_ts <- intersect(ts_kizuka, ts_gd)
  if (length(common_enz) > 0 && length(common_ts) > 0) {
    sub_kizuka <- res_kizuka[common_enz, common_ts, drop = FALSE]
    sub_gd <- res_clean[common_enz, common_ts, drop = FALSE]
    n_agree <- sum(sub_kizuka == sub_gd, na.rm = TRUE)
    n_total <- sum(!is.na(sub_kizuka) & !is.na(sub_gd))
    cat("Consistency on overlapping entries:", n_agree, "/", n_total,
        "(", round(n_agree / n_total * 100, 1), "%)\n")
  }
} else {
  cat("output/res_clean.rds not found; skipping Kizuka comparison\n")
}

t_end <- Sys.time()
cat("\n=== Total runtime:",
    round(difftime(t_end, t_start, units = "mins"), 1), "minutes ===\n")
cat("=== Processing complete ===\n")

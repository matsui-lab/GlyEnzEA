# Compare regenerated HEK293 / robustness outputs with the frozen copies.
fz <- "provenance/frozen_outputs"
a <- readRDS("output/hek293/resGSEA_by_ko.rds"); b <- readRDS(file.path(fz, "hek293/resGSEA_by_ko.rds"))
for (k in names(b)) { x <- as.data.frame(a[[k]]); y <- as.data.frame(b[[k]]); x <- x[order(x$pathway), ]; y <- y[order(y$pathway), ]
  cat(sprintf("%-20s sets %d/%d  NES max|d| %.1e  pval max|d| %.1e  q max|d| %.1e\n", k, nrow(x), nrow(y),
              max(abs(x$NES - y$NES)), max(abs(x$pval - y$pval)), max(abs(x$qvalue - y$qvalue)))) }
for (f in c("hek293/res_clean.rds", "hek293/res_clean_group.rds", "hek293/gmt_enzyme_grouped.rds", "hek293/enzyme_pattern_groups.rds", "hek293/rowmeta.rds")) {
  x <- readRDS(file.path("output", f)); y <- readRDS(file.path(fz, f))
  cat(sprintf("%-40s %s\n", f, if (identical(x, y)) "IDENTICAL" else if (isTRUE(all.equal(unname(x), unname(y)))) "same values (names/attrs differ)" else "DIFFERENT")) }
for (f in c("gsea_network_null.tsv", "TableY_minsize_sensitivity.tsv", "FigureS_three_null_MGAT3.tsv", "FigureS_nperm_sensitivity.tsv", "TableSx_jaccard_matrix.tsv", "hek293/suppT_gsea_all.tsv")) {
  x <- readLines(file.path("output", f)); y <- readLines(file.path(fz, f))
  cat(sprintf("%-40s %s\n", f, if (identical(x, y)) "IDENTICAL" else sprintf("differs in %d of %d lines", sum(x != y[seq_along(x)]), length(x)))) }

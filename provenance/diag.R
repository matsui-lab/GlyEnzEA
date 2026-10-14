d <- function(f) { a <- readRDS(file.path("output", f)); b <- readRDS(file.path("provenance/frozen_outputs", f))
  cat("\n##", f, "| class:", class(a)[1], "vs", class(b)[1], "| dim:", paste(dim(a), collapse="x"), "vs", paste(dim(b), collapse="x"), "\n")
  ae <- all.equal(as.data.frame(a), as.data.frame(b)); if (isTRUE(ae)) cat("  as.data.frame: all.equal TRUE\n") else print(head(ae, 6))
  if (is.matrix(a)) { cat("  same dimnames:", identical(dimnames(a), dimnames(b)), " same values:", isTRUE(all.equal(unname(a), unname(b))), "\n")
    if (!identical(dimnames(a), dimnames(b))) { cat("  rows only in new:", setdiff(rownames(a), rownames(b)), "| cols diff n:", length(setdiff(colnames(a), colnames(b))), "\n") } }
}
for (f in c("triplet_enzdb.rds","triplet_glygen.rds","glycomeRes.rds","triplet.rds","res_clean.rds","glycomeRes_group.rds")) d(f)

# Compare the regenerated Glyco-Decipher library mapping (step3b) with the frozen copy.
fz <- "provenance/frozen_outputs"
a <- readRDS("output/res_clean_glycodecipher.rds"); b <- readRDS(file.path(fz, "res_clean_glycodecipher.rds"))
cat("dims:", paste(dim(a), collapse = "x"), "vs", paste(dim(b), collapse = "x"), "\n")
cat("rows identical:", identical(rownames(a), rownames(b)), " values identical (positional):", identical(unname(a), unname(b)), "\n")
pa <- sort(apply(a, 2, function(v) paste(as.integer(v), collapse = ""))); pb <- sort(apply(b, 2, function(v) paste(as.integer(v), collapse = "")))
cat("column patterns identical as multiset:", identical(unname(pa), unname(pb)), "\n")
cat("topology names differing:", sum(!(colnames(a) %in% colnames(b))), "of", ncol(a), "\n")
m <- readRDS("output/glycodecipher_toposeq_map.rds"); n <- readRDS(file.path(fz, "glycodecipher_toposeq_map.rds"))
cat("toposeq map rows:", nrow(m), "vs", nrow(n), " converted (non-NA):", sum(!is.na(m$toposeq)), "vs", sum(!is.na(n$toposeq)), "\n")

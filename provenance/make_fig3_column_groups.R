# One-off (2026-09-28): derive the Figure 3 column groups (G1..G29) as enzyme-mapping patterns,
# so the figure script does not depend on the exact topology strings (canonical strings can change
# between versions of R/utils_formatting.R while the mapping matrix stays the same).
# Inputs: the frozen res_clean.rds and Figure2C_glycan_legend.csv from the April 2026 run.
rc <- readRDS("provenance/frozen_outputs/res_clean.rds")
gl <- read.csv("provenance/frozen_outputs/Figure2C_glycan_legend.csv", stringsAsFactors = FALSE)
pat <- function(cols) { p <- unique(apply(rc[, cols, drop = FALSE], 2, function(v) paste(as.integer(v), collapse = ""))); stopifnot(length(p) == 1); p }
out <- data.frame(GroupLabel = gl$GroupLabel,
                  Pattern = vapply(strsplit(gl$Members, " \\| "), pat, ""),
                  NTopologies = vapply(strsplit(gl$Members, " \\| "), length, 1L), stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(out$Pattern), sum(out$NTopologies) == 97)
attr(out, "enzyme_order") <- rownames(rc)
write.csv(out, "data/figure3_column_groups.csv", row.names = FALSE)
writeLines(rownames(rc), "data/figure3_enzyme_order.txt")
cat("wrote data/figure3_column_groups.csv (", nrow(out), "groups ) and data/figure3_enzyme_order.txt\n")

# Compare regenerated outputs with the frozen copies (provenance/frozen_outputs).
cmp <- function(f, fun = identical) {
  a <- readRDS(file.path("output", f)); b <- readRDS(file.path("provenance/frozen_outputs", f))
  cat(sprintf("%-40s %s\n", f, if (fun(a, b)) "IDENTICAL" else "DIFFERENT"))
}
for (f in c("triplet_enzdb.rds", "triplet_glygen.rds", "glycomeRes.rds", "triplet.rds", "res_clean.rds",
            "res_clean_group.rds", "glycomeRes_group.rds", "glycan2enzMap_enzyme_grouped.rds",
            "enzyme_pattern_groups_list.rds", "resORA_down.rds")) cmp(f)
g <- as.data.frame(readRDS("output/resGSEA.rds")); h <- as.data.frame(readRDS("provenance/frozen_outputs/resGSEA.rds"))
g <- g[order(g$pathway), ]; h <- h[order(h$pathway), ]
cat(sprintf("%-40s NES max|diff| = %.2e, pval max|diff| = %.2e, qvalue max|diff| = %.2e\n", "resGSEA.rds",
            max(abs(g$NES - h$NES)), max(abs(g$pval - h$pval)), max(abs(g$qvalue - h$qvalue))))
print(g[g$pathway == "MGAT3", c("pathway", "NES", "pval", "padj", "qvalue")])

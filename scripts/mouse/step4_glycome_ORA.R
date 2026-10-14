library(fgsea)
library(qvalue)
library(dplyr)

# Load enzyme-grouped GMT (non-redundant enzyme pattern groups x 33 glycan groups)
# Enzyme grouping reduces redundant hypothesis tests (analogous to REVIGO for GO terms).
# Glycan groups remain individual as each is an independently quantified molecular species.
gmt_enzyme_grouped <- readRDS("output/glycan2enzMap_enzyme_grouped.rds")
enzyme_pattern_groups <- readRDS("output/enzyme_pattern_groups_list.rds")
glycomeRes_group <- readRDS("output/glycomeRes_group.rds")   # 33 ranked observed glycan species (step3)


## ------------------------------------------------------------
## ORA (Fisher's exact test)
## ------------------------------------------------------------

sigGlycan <- glycomeRes_group$group[glycomeRes_group$diff %in% c("increased", "lost")] %>% unique()
sigGlycan_up <- glycomeRes_group$group[glycomeRes_group$diff == "increased"] %>% unique()
sigGlycan_down <- glycomeRes_group$group[glycomeRes_group$diff == "lost"] %>% unique()

runORA_fisher <- function(sigGlycan, gmt, all_glycan) {
  background <- unique(c(unlist(gmt), all_glycan))
  N <- length(background)
  K <- length(intersect(sigGlycan, background))
  results <- do.call(rbind, lapply(names(gmt), function(enz) {
    pathway_members <- gmt[[enz]]
    n <- length(intersect(pathway_members, background))
    k <- length(intersect(pathway_members, sigGlycan))
    mat <- matrix(c(max(k, 0), max(K - k, 0),
                     max(n - k, 0), max(N - n - K + k, 0)), nrow = 2)
    ft <- fisher.test(mat, alternative = "greater")
    expected <- (K * n) / N
    data.frame(
      pathway = enz, k = k, n = n, K = K, N = N,
      expected = round(expected, 2),
      pvalue = ft$p.value,
      OR = as.numeric(ft$estimate),
      fold_enrichment = ifelse(expected > 0, round(k / expected, 2), Inf),
      stringsAsFactors = FALSE
    )
  }))
  results$padj <- p.adjust(results$pvalue, method = "BH")
  results$qvalue <- tryCatch(
    qvalue(results$pvalue, pi0.method = "bootstrap")$qvalues,
    error = function(e) p.adjust(results$pvalue, method = "BH")
  )
  results
}

resORA_up <- runORA_fisher(sigGlycan_up, gmt_enzyme_grouped, glycomeRes_group$group)
resORA_down <- runORA_fisher(sigGlycan_down, gmt_enzyme_grouped, glycomeRes_group$group)

resORA_up %>% filter(pvalue < 0.05)
resORA_down %>% filter(pvalue < 0.05)


## ------------------------------------------------------------
## GSEA (enzyme-grouped pathways)
## ------------------------------------------------------------

lfc_rank <- setNames(glycomeRes_group$logFC, glycomeRes_group$group)
lfc_rank <- sort(lfc_rank)

set.seed(42)
resGSEA <- fgsea(pathways = gmt_enzyme_grouped,
                 stats    = lfc_rank,
                 minSize  = 1,
                 maxSize  = 500,
                 nPermSimple = 50000)

# Storey's q-value (more appropriate FDR than BH when pi0 < 1)
resGSEA$qvalue <- tryCatch(
  qvalue(resGSEA$pval, pi0.method = "bootstrap")$qvalues,
  error = function(e) resGSEA$padj  # fallback to BH
)

resGSEA <- resGSEA %>% arrange(pval)
resGSEA

(resGSEA$pval < 0.05) %>% table()
(resGSEA$qvalue < 0.05) %>% table()
(resGSEA$qvalue < 0.1) %>% table()
(resGSEA$qvalue < 0.25) %>% table()

resGSEA$pathway[resGSEA$qvalue < 0.05]
resGSEA$pathway[resGSEA$qvalue < 0.25]

enriched <- resGSEA %>%
  filter(pval < 0.05) %>%
  select(pathway, NES, pval, padj, qvalue, size)

## ------------------------------------------------------------
## Save (added for the paper repository; the original script kept these in the session)
## ------------------------------------------------------------
saveRDS(resGSEA,    "output/resGSEA.rds")
saveRDS(resORA_down, "output/resORA_down.rds")
saveRDS(resORA_up,   "output/resORA_up.rds")
cat("saved output/resGSEA.rds, resORA_down.rds, resORA_up.rds\n")

## ------------------------------------------------------------
## GSEA network null: degree-preserving bipartite randomization
##
## Tests whether GSEA enrichment of enzyme-grouped pathways depends
## on the specific enzyme<->glycan topology, or can be explained by
## the marginal degree distributions alone.
##
## Null model: curveball algorithm (Strona et al., Nat Commun 2014)
## preserves both row degree (# enzymes/glycan) and column degree
## (# glycans/enzyme), shuffling only edge endpoints.
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(fgsea)
  library(qvalue)
  library(dplyr)
})

set.seed(42)

## ---- Inputs --------------------------------------------------
gmt_enzyme_grouped <- readRDS("output/glycan2enzMap_enzyme_grouped.rds")
glycomeRes_group   <- readRDS("output/glycomeRes_group.rds")

lfc_rank <- setNames(glycomeRes_group$logFC, glycomeRes_group$group)
## sort decreasing to match fgsea's internal convention (positive ES =
## enrichment among high-stat glycans). calcGseaStat respects given order,
## so this is the convention that aligns ES signs with fgsea_NES.
lfc_rank <- sort(lfc_rank, decreasing = TRUE)

## ---- Observed fgsea (recompute so we share the same stats) ----
resGSEA_obs <- fgsea(pathways = gmt_enzyme_grouped,
                     stats    = lfc_rank,
                     minSize  = 1,
                     maxSize  = 500,
                     nPermSimple = 50000)
resGSEA_obs$qvalue <- tryCatch(
  qvalue::qvalue(resGSEA_obs$pval, pi0.method = "bootstrap")$qvalues,
  error = function(e) resGSEA_obs$padj
)

## ---- Bipartite matrix (glycan rows x enzyme cols) -------------
universe <- names(lfc_rank)
M <- sapply(gmt_enzyme_grouped, function(s) as.integer(universe %in% s))
rownames(M) <- universe
stopifnot(all(colSums(M) == sapply(gmt_enzyme_grouped, length)))

## ---- Curveball randomizer ------------------------------------
## One curveball "trade": pick two rows; among columns where exactly one
## row has a 1, randomly redistribute. Preserves all row & column sums.
curveball_swap <- function(M, n_trades = 5L * sum(M)) {
  nr <- nrow(M)
  ## represent rows as integer position vectors of their 1-columns
  rows <- lapply(seq_len(nr), function(i) which(M[i, ] == 1L))
  for (t in seq_len(n_trades)) {
    ij <- sample.int(nr, 2L)
    A <- rows[[ij[1]]]; B <- rows[[ij[2]]]
    common <- intersect(A, B)
    onlyA  <- setdiff(A, common)
    onlyB  <- setdiff(B, common)
    L <- length(onlyA) + length(onlyB)
    if (L == 0L || length(onlyA) == 0L || length(onlyB) == 0L) next
    pool <- c(onlyA, onlyB)
    newA <- sample(pool, length(onlyA))
    newB <- setdiff(pool, newA)
    rows[[ij[1]]] <- sort(c(common, newA))
    rows[[ij[2]]] <- sort(c(common, newB))
  }
  out <- matrix(0L, nrow = nr, ncol = ncol(M),
                dimnames = dimnames(M))
  for (i in seq_len(nr)) out[i, rows[[i]]] <- 1L
  out
}

## sanity: degrees preserved
.M1 <- curveball_swap(M, n_trades = 100L)
stopifnot(identical(rowSums(.M1), rowSums(M)))
stopifnot(identical(colSums(.M1), colSums(M)))
cat("[curveball] degree preservation: OK\n")

## ---- Network null: B randomizations, ES per enzyme ----------
B <- 1000L
enz_names <- colnames(M)

## observed ES via fgsea internal helper
es_one <- function(set_glycans) {
  idx <- match(set_glycans, names(lfc_rank))
  idx <- idx[!is.na(idx)]
  if (length(idx) < 1L) return(NA_real_)
  fgsea::calcGseaStat(stats = lfc_rank, selectedStats = idx,
                      gseaParam = 1, returnAllExtremes = FALSE)
}

ES_obs <- vapply(gmt_enzyme_grouped, es_one, numeric(1))

ES_null <- matrix(NA_real_, nrow = length(enz_names), ncol = B,
                  dimnames = list(enz_names, paste0("perm", seq_len(B))))

t0 <- Sys.time()
for (b in seq_len(B)) {
  Mb <- curveball_swap(M)
  for (e in seq_along(enz_names)) {
    set_b <- rownames(Mb)[Mb[, e] == 1L]
    ES_null[e, b] <- es_one(set_b)
  }
  if (b %% 100 == 0) cat(sprintf("[null] %d/%d (%.1f min)\n", b, B,
                                 as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---- Empirical p-values --------------------------------------
## Sign-aware: if ES_obs > 0, p = P(ES_null >= ES_obs); if < 0, P(ES_null <= ES_obs).
## (Same one-sided semantics as fgsea NES.)
emp_p <- vapply(enz_names, function(e) {
  obs <- ES_obs[e]
  if (is.na(obs)) return(NA_real_)
  v <- ES_null[e, ]; v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_real_)
  if (obs >= 0) (sum(v >= obs) + 1) / (length(v) + 1)
  else          (sum(v <= obs) + 1) / (length(v) + 1)
}, numeric(1))

## network-NES: ES_obs normalized by mean(|ES_null| with matching sign)
emp_nes <- vapply(enz_names, function(e) {
  obs <- ES_obs[e]
  if (is.na(obs)) return(NA_real_)
  v <- ES_null[e, ]; v <- v[!is.na(v)]
  same <- if (obs >= 0) v[v >= 0] else v[v < 0]
  if (length(same) == 0L) return(NA_real_)
  obs / mean(abs(same))
}, numeric(1))

emp_q <- tryCatch(
  qvalue::qvalue(emp_p, pi0.method = "bootstrap")$qvalues,
  error = function(e) p.adjust(emp_p, method = "BH")
)

## ---- Comparison table ----------------------------------------
cmp <- data.frame(
  pathway     = enz_names,
  size        = colSums(M)[enz_names],
  ES_obs      = ES_obs[enz_names],
  fgsea_NES   = resGSEA_obs$NES[match(enz_names, resGSEA_obs$pathway)],
  fgsea_pval  = resGSEA_obs$pval[match(enz_names, resGSEA_obs$pathway)],
  fgsea_qval  = resGSEA_obs$qvalue[match(enz_names, resGSEA_obs$pathway)],
  netnull_NES = emp_nes[enz_names],
  netnull_pval = emp_p[enz_names],
  netnull_qval = emp_q[enz_names],
  stringsAsFactors = FALSE,
  row.names = NULL
) |> dplyr::arrange(fgsea_pval)

## ---- Save & print --------------------------------------------
saveRDS(list(comparison = cmp,
             ES_null    = ES_null,
             ES_obs     = ES_obs,
             B          = B,
             seed       = 42L,
             fgsea_full = resGSEA_obs),
        "output/gsea_network_null.rds")

write.table(cmp, "output/gsea_network_null.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== fgsea (score permutation) vs network null (degree-preserving) ===\n")
print(cmp, digits = 3)

cat("\nCalls at fgsea p<0.05 :", sum(cmp$fgsea_pval < 0.05, na.rm = TRUE), "\n")
cat("Calls at netnull p<0.05:", sum(cmp$netnull_pval < 0.05, na.rm = TRUE), "\n")
cat("Calls at fgsea q<0.25  :", sum(cmp$fgsea_qval < 0.25, na.rm = TRUE), "\n")
cat("Calls at netnull q<0.25:", sum(cmp$netnull_qval < 0.25, na.rm = TRUE), "\n")

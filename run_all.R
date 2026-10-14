#!/usr/bin/env Rscript
# Entry point: reproduce the GlyEnzEA paper analyses, figures and supplementary tables from the
# files in this repository. Run from the repository root:  Rscript run_all.R [stage ...]
# Stages (default = all):
#   knowledge   step1: reaction triplets from GlycoEnzDB and GlyGen (uses cache/ for GlyTouCan->WURCS)
#   mouse       step2-4: Mgat3-KO mouse brain glycome -> mapping -> enrichment (+ curveball null)
#   robustness  minSize / nPermSimple / three-null comparisons, Jaccard matrix, Table S1
#   hek293      Huang 2021 seven knockouts (00-02 + supplementary tables)
#   library     step3b: Glyco-Decipher library mapping (slow; numbers in Results 3.1)
#   figures     Fig. 3-5 into paper/figs/
# Everything is offline: API responses are frozen in cache/. Set GLYENZEA_FETCH=1 to allow new API calls.
stages <- list(
  knowledge  = c("scripts/mouse/step1_ENZDB_to_graph.R", "scripts/mouse/step1_GlyGen_to_graph.R"),
  mouse      = c("scripts/mouse/step2_glycome.R", "scripts/mouse/step3_build_pk.R", "scripts/mouse/step4_glycome_ORA.R",
                 "scripts/mouse/step4b_gsea_network_null.R"),
  robustness = c("scripts/robustness/minsize_sensitivity.R", "scripts/robustness/three_null_comparison.R",
                 "scripts/robustness/figS_stat_robustness.R", "scripts/robustness/figS3_jaccard_heatmap.R",
                 "scripts/robustness/regen_TableS1.R"),
  hek293     = c("scripts/hek293/00_verify_toposeq.R", "scripts/hek293/01_adapter.R",
                 "scripts/hek293/02_run_step3_step4.R", "scripts/hek293/build_supp_data.R"),
  library    = c("scripts/mouse/step3b_glycodecipher_mapping.R"),
  figures    = c("scripts/figures/fig3_matrix.R", "scripts/figures/fig4_mgat3.R", "scripts/figures/fig5_hek293.R"))
args <- commandArgs(trailingOnly = TRUE)
run <- if (length(args) == 0) names(stages) else args
stopifnot(all(run %in% names(stages)), file.exists("R/utils_formatting.R"))
dir.create("logs", showWarnings = FALSE)
for (st in run) for (s in stages[[st]]) {
  cat(sprintf("[%s] %s ... ", st, s)); t0 <- Sys.time()
  log <- file.path("logs", paste0(sub("\\.R$", "", basename(s)), ".log"))
  rc <- system2("Rscript", s, stdout = log, stderr = log)
  cat(sprintf("%s (%.0f s)\n", if (rc == 0) "ok" else "FAILED", as.numeric(Sys.time() - t0, units = "secs")))
  if (rc != 0) stop("see ", log)
}
cat("done\n")

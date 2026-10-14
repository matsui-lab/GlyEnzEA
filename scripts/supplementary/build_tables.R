#!/usr/bin/env Rscript
# Supplementary Tables S1 (mouse, 16 glycan sets) and S2 (HEK293, 7 perturbations x 22 sets) from output/.
# Run from the repository root after `Rscript run_all.R mouse robustness hek293`.
suppressPackageStartupMessages({ library(dplyr); library(qvalue) })
out <- "paper/supplementary/build"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- Table S1: composition + fgsea + Storey/BH + curveball null + minSize sensitivity ----
grp  <- read.csv("output/TableS1_enzyme_pattern_groups.csv", stringsAsFactors = FALSE)
res  <- as.data.frame(readRDS("output/resGSEA.rds")); res$leadingEdge <- NULL
st   <- tryCatch(qvalue(res$pval, pi0.method = "bootstrap")$qvalues, error = function(e) NULL)
res$storey_q <- if (is.null(st)) NA_real_ else st
res$BH       <- p.adjust(res$pval, "BH")
res$method   <- if (is.null(st)) "BH*" else "Storey"
stopifnot(isTRUE(all.equal(res$storey_q, res$qvalue)))            # same as the value used in the paper
nn   <- read.delim("output/gsea_network_null.tsv", check.names = FALSE)
ms   <- read.delim("output/TableY_minsize_sensitivity.tsv", check.names = FALSE)
names(ms) <- c("pathway", "P_minSize1", "P_minSize2", "P_minSize3")[seq_len(ncol(ms))]
t1 <- grp %>% rename(pathway = Pattern_group) %>%
  left_join(res %>% select(pathway, size, NES, nominal_P = pval, storey_q, BH_adjusted_P = BH, method), by = "pathway") %>%
  left_join(nn %>% select(pathway, curveball_NES = netnull_NES, curveball_P = netnull_pval), by = "pathway") %>%
  left_join(ms, by = "pathway") %>%
  arrange(nominal_P) %>%
  transmute(`Glycan set` = pathway, `Member enzymes` = Member_enzymes, Class = Enzyme_class,
            `Set size` = size, NES = round(NES, 2), `Nominal P` = signif(nominal_P, 2), `Storey q` = signif(storey_q, 2),
            `BH P` = signif(BH_adjusted_P, 2), Correction = method,
            `Curveball NES` = round(curveball_NES, 2), `Curveball P` = signif(curveball_P, 2),
            `P (minSize 2)` = signif(P_minSize2, 2), `P (minSize 3)` = signif(P_minSize3, 2))
stopifnot(nrow(t1) == 16, !any(is.na(t1$NES)))
write.csv(t1, file.path(out, "TableS1_mouse_sets.csv"), row.names = FALSE)

# ---- Table S2: HEK293, all perturbations; which correction was used per perturbation ----
# resGSEA_by_ko.rds holds the exact fgsea output per perturbation with padj (BH) and qvalue (Storey q, or padj
# when pi0 estimation failed: 02_run_step3_step4.R). A fallback is therefore qvalue identical to padj.
rb <- readRDS("output/hek293/resGSEA_by_ko.rds")
h  <- read.delim("output/hek293/suppT_gsea_all.tsv", check.names = FALSE)      # KO labels and target flags as reported
t2 <- bind_rows(lapply(rb, function(d) {
  d <- as.data.frame(d); d$leadingEdge <- NULL
  d$method <- if (isTRUE(all.equal(d$qvalue, d$padj))) "BH*" else "Storey"
  d
}))
ko_label <- c(MGAT1_KO = "MGAT1-KO", MGAT2_KO = "MGAT2-KO", MGAT4A_4B_KO = "MGAT4A&4B-KO", MGAT5_KO = "MGAT5-KO",
              T_KO_MANA1_A2_B1_ = "T-KO (MAN1A1&A2&B1)", B4GALNT3_4_KO = "B4GALNT3&4-KO", SLC35C1_KO = "SLC35C1-KO")
t2$KO <- unname(ko_label[t2$ko]); stopifnot(!any(is.na(t2$KO)), setequal(unique(t2$KO), unique(h$KO)))
t2 <- t2 %>% left_join(h %>% select(KO, pathway, is_target), by = c("KO", "pathway")) %>% arrange(KO, pval) %>%
  transmute(Perturbation = KO, `Glycan set` = pathway, `Target set` = ifelse(is_target, "yes", ""), `Set size` = size,
            NES = round(NES, 2), `Nominal P` = signif(pval, 2), q = signif(qvalue, 2), `BH P` = signif(padj, 2), Correction = method)
stopifnot(nrow(t2) == 7 * 22, !any(is.na(t2$`Target set`)))
write.csv(t2, file.path(out, "TableS2_hek293_all.csv"), row.names = FALSE)
cat("Table S1:", nrow(t1), "rows; Table S2:", nrow(t2), "rows\n")
print(unique(t2[, c("Perturbation", "Correction used")]), row.names = FALSE)

#!/usr/bin/env Rscript
# Figure 5 (manuscript v7): seven HEK293 perturbations (Huang et al. 2021).
# Data: output/hek293/resGSEA_by_ko.rds (02_run_step3_step4.R; 22 glycan sets per KO).
# Panel content = the original fig_figure4_7panel.R (seven lowest-q sets plus the target);
# layout: 3 columns x 3 rows (A-G + shared legend), rank/q headings; 2 columns until 2026-10-08.
# Output: paper/figs/fig5.pdf (700 x 358 pt).
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(data.table); library(grid); library(cowplot) })

cp <- function(a, b) { m <- min(nchar(a), nchar(b)); k <- 0L
  while (k < m && substr(a, k+1, k+1) == substr(b, k+1, k+1)) k <- k + 1L; substr(a, 1, k) }
collapse_label <- function(s) {
  toks <- strsplit(s, "/", fixed = TRUE)[[1]]
  if (length(toks) < 2 || !all(grepl("^[A-Z][A-Z0-9]*[0-9]$", toks))) return(s)
  groups <- list(toks[1])
  for (i in 2:length(toks)) { g <- groups[[length(groups)]]
    if (nchar(cp(g[1], toks[i])) >= 3) groups[[length(groups)]] <- c(g, toks[i]) else groups[[length(groups) + 1]] <- toks[i] }
  paste(vapply(groups, function(g) { if (length(g) == 1) return(g[1]); stem <- Reduce(cp, g)
    paste(c(g[1], substring(g[-1], nchar(stem) + 1)), collapse = "/") }, ""), collapse = "/")
}
wrap_slash <- function(s, maxw = 16) {
  if (nchar(s) <= maxw) return(s)
  toks <- strsplit(s, "(?<=/)", perl = TRUE)[[1]]; if (length(toks) < 2) return(s)
  cum <- cumsum(nchar(toks)); i <- min(which.min(abs(cum - sum(nchar(toks))/2)), length(toks) - 1)
  paste0(paste(toks[1:i], collapse = ""), "\n", paste(toks[(i+1):length(toks)], collapse = ""))
}

r <- readRDS("output/hek293/resGSEA_by_ko.rds")
kos <- c("MGAT1_KO", "MGAT2_KO", "MGAT4A_4B_KO", "MGAT5_KO", "T_KO_MANA1_A2_B1_", "B4GALNT3_4_KO", "SLC35C1_KO")
ko_lab <- c(MGAT1_KO = "MGAT1 KO", MGAT2_KO = "MGAT2 KO", MGAT4A_4B_KO = "MGAT4A/B KO", MGAT5_KO = "MGAT5 KO",
            T_KO_MANA1_A2_B1_ = "MAN1A1/A2/B1 KO", B4GALNT3_4_KO = "B4GALNT3/4 KO", SLC35C1_KO = "SLC35C1 KO")
target <- c(MGAT1_KO = "MGAT1/2/5B", MGAT2_KO = "MGAT1/2/5B", MGAT4A_4B_KO = "MGAT4/5", MGAT5_KO = "MGAT4/5",
            T_KO_MANA1_A2_B1_ = "MAN1A1/MAN1A2/MAN1C1", B4GALNT3_4_KO = "B3GALNT2/B4GALNT3/B4GALNT4", SLC35C1_KO = "FUT8")
TOPN <- 7; tcol <- "#1F4E9C"

info <- rbindlist(lapply(kos, function(k) { x <- as.data.frame(r[[k]]); stopifnot(nrow(x) == 22)
  data.table(ko = k, rank = rank(x$pval, ties.method = "min")[x$pathway == target[[k]]],
             q = x$qvalue[x$pathway == target[[k]]], method = if (all(abs(x$qvalue - x$padj) < 1e-12)) "BH" else "Storey") }))
info[, letter := LETTERS[seq_len(.N)]]
print(info)

panel <- function(k, letter) {
  x <- as.data.frame(r[[k]]); x <- x[order(x$qvalue), ]
  keep <- union(head(x$pathway, TOPN), target[[k]]); x <- x[x$pathway %in% keep, ]
  d <- data.table(pathway = x$pathway, NES = x$NES, q = x$qvalue, is_target = x$pathway == target[[k]])
  d[, star := ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", ifelse(q < 0.10, "+", ""))))]
  abbr <- c("ST3GAL3/ST3GAL4/ST3GAL6/ST6GAL1/ST6GAL2/ST3GAL1" = "ST3GAL/ST6GAL",
            "B3GALNT2/B4GALNT3/B4GALNT4" = "B3GALNT2/\nB4GALNT3/4")   # abbreviations defined in the caption
  d[, disp := ifelse(pathway %in% names(abbr), abbr[pathway], vapply(vapply(pathway, collapse_label, ""), wrap_slash, ""))]
  d <- d[order(NES)]; d[, row := factor(pathway, levels = pathway)]
  i <- info[ko == k]
  ggplot(d, aes(NES, row, fill = NES)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_col(width = 0.68, colour = "grey60", linewidth = 0.2) +
    geom_text(aes(label = star, x = NES + ifelse(NES >= 0, 0.1, -0.1), hjust = ifelse(NES >= 0, 0, 1)),
              vjust = 0.75, size = 2.2, colour = "grey15") +
    geom_text(aes(x = -2.75, y = row, label = disp), fontface = ifelse(d$is_target, "bold", "plain"),
              colour = ifelse(d$is_target, tcol, "black"), hjust = 1, size = 1.8, lineheight = 0.8) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2.4, 2.4),
                         oob = scales::squish, name = "NES") +
    scale_x_continuous(breaks = -2:2) +
    coord_cartesian(xlim = c(-2.6, 2.6), clip = "off") +
    labs(x = "NES", y = NULL, title = bquote(bold(.(letter))~~.(ko_lab[[k]])),
         subtitle = sprintf("Target: rank %d/22, q = %.3f", i$rank, i$q)) +
    theme_minimal(base_size = 6.5) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"),
          axis.text.y = element_blank(), axis.text.x = element_text(size = 5.5), axis.title.x = element_text(size = 6),
          plot.title = element_text(size = 7, margin = margin(b = 0)), plot.subtitle = element_text(size = 6, colour = "grey25", margin = margin(b = 2)),
          legend.position = "none", plot.margin = margin(3, 6, 2, 40))
}
plots <- Map(panel, kos, info$letter)

# shared legend: colour bar + text
leg_plot <- plots[[1]] + theme(legend.position = "right", legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
                               legend.key.height = unit(9, "pt"), legend.key.width = unit(5, "pt"))
cbar <- cowplot::get_plot_component(leg_plot, "guide-box-right", return_all = TRUE)
storey <- paste(info[method == "Storey", letter], collapse = ", ")
bhs <- info[method == "BH", letter]
txt <- c("Bold labels: target glycan sets",
         "Rank: nominal P-value rank among 22 glycan sets",
         sprintf("q: Storey (%s); BH fallback (%s)", storey, paste(bhs, collapse = ", ")),
         "* q < 0.05   ** q < 0.01   *** q < 0.001   + q < 0.10",
         "All seven perturbations are shown")
tg <- textGrob(paste(txt, collapse = "\n"), x = 0.02, y = 0.95, just = c("left", "top"), gp = gpar(fontsize = 5.6, lineheight = 1.25))
leg <- grobTree(grobTree(cbar, vp = viewport(x = 0.1, y = 0.55, width = 0.15, height = 0.6)),
                grobTree(tg, vp = viewport(x = 0.62, y = 0.5, width = 0.72, height = 0.9)))

fig <- wrap_plots(c(plots, list(wrap_elements(leg))), ncol = 3, byrow = TRUE)
dir.create("paper/figs", showWarnings = FALSE)
ggsave("paper/figs/fig5.pdf", fig, width = 700/72, height = 358/72, device = pdf, useDingbats = FALSE)
cat("saved paper/figs/fig5.pdf\n")

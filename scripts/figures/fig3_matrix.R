#!/usr/bin/env Rscript
# Figure 3 (manuscript v7): 16 distinct enzyme-centered glycan sets x 29 display groups.
# Data: output/res_clean.rds (48 enzymes x 97 candidate topologies, step3_build_pk.R),
#       output/enzyme_pattern_groups_list.rds (16 enzyme groups with identical rows),
#       data/figure3_column_groups.csv (29 column groups as enzyme-mapping patterns with the G labels
#       of the typeset figure; derived once from the April 2026 run, see provenance/make_fig3_column_groups.R).
# Row/column order = the clustering used in Figure2.R panel C (binary distance, ward.D2,
# ComplexHeatmap default dendrogram reordering). Layout follows the 2026-09-16 typeset figure.
# Output: paper/figs/fig3.pdf (280 x 277 pt). ASCII-only text (pdf device, no cairo).
suppressPackageStartupMessages({
  library(ggplot2); library(ComplexHeatmap); library(patchwork); library(cowplot)
})
source("R/utils_formatting.R")

res_clean <- readRDS("output/res_clean.rds")
pg        <- readRDS("output/enzyme_pattern_groups_list.rds")
cg        <- read.csv("data/figure3_column_groups.csv", stringsAsFactors = FALSE, colClasses = c(Pattern = "character"))   # G label -> enzyme-mapping pattern
enz_order <- readLines("data/figure3_enzyme_order.txt")
stopifnot(setequal(enz_order, rownames(res_clean)))
res_clean <- res_clean[enz_order, ]

# ---- collapse rows (enzyme groups with identical rows) ----
stopifnot(all(vapply(pg, function(m) all(apply(res_clean[m, , drop = FALSE], 2, function(v) length(unique(v)) == 1)), TRUE)))
# ---- collapse columns by enzyme-mapping pattern; label with the G numbers of the typeset figure ----
col_pat <- apply(res_clean, 2, function(v) paste(as.integer(v), collapse = ""))
stopifnot(setequal(unique(col_pat), cg$Pattern))
col_members <- lapply(cg$Pattern, function(p) names(col_pat)[col_pat == p])
stopifnot(all(lengths(col_members) == cg$NTopologies))
mat <- sapply(col_members, function(m) res_clean[vapply(pg, `[`, "", 1), m[1]])
mat <- matrix(as.numeric(mat), nrow = length(pg))
rownames(mat) <- vapply(pg, extract_enzyme_group_label, "")
colnames(mat) <- cg$GroupLabel
n_enz <- vapply(pg, length, 1L)
stopifnot(nrow(mat) == 16, ncol(mat) == 29, sum(n_enz) == 48)

# ---- ordering (same call as Figure2.R) ----
ht <- draw(Heatmap(mat, clustering_distance_rows = "binary", clustering_distance_columns = "binary",
                   clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                   show_heatmap_legend = FALSE))
ro <- row_order(ht); co <- column_order(ht)
row_lab <- rownames(mat)[ro]; col_lab <- colnames(mat)[co]
cat("row order:", paste(row_lab, collapse = ", "), "\n")
cat("col order:", paste(col_lab, collapse = " "), "\n")

short <- function(s) {   # MAN1A1/MAN1A2/MAN1C1 -> MAN1A1/A2/C1
  toks <- strsplit(s, "/", fixed = TRUE)[[1]]
  if (length(toks) < 2 || !all(grepl("^[A-Z][A-Z0-9]*[0-9]$", toks))) return(s)
  cp <- function(a, b) { k <- 0; while (k < min(nchar(a), nchar(b)) && substr(a, k+1, k+1) == substr(b, k+1, k+1)) k <- k + 1; substr(a, 1, k) }
  stem <- Reduce(cp, toks)
  if (nchar(stem) >= 3) paste(c(toks[1], substring(toks[-1], nchar(stem) + 1)), collapse = "/") else s
}
hl <- c("MGAT3", "FUT8")
lab_txt <- vapply(row_lab, short, ""); lab_txt[lab_txt == "ST3GAL/ST6GAL/CMAS"] <- "ST3GAL/ST6GAL\n/CMAS"
lab_df <- data.frame(row = factor(row_lab, levels = rev(row_lab)), lab = lab_txt,
                     col = ifelse(row_lab %in% hl, "#1F5FA8", "black"), face = ifelse(row_lab %in% hl, "bold", "plain"))
long <- data.frame(row = factor(rep(row_lab, times = ncol(mat)), levels = rev(row_lab)),
                   col = factor(rep(col_lab, each = nrow(mat)), levels = col_lab),
                   mapped = as.vector(mat[ro, co]) == 1)
bars <- data.frame(row = factor(row_lab, levels = rev(row_lab)), n = n_enz[ro],
                   cls = ifelse(n_enz[ro] == 1, "1 enzyme", "2-10 enzymes"))

col_map <- c(`TRUE` = "#4A2C2A", `FALSE` = "#F4F1DC")
bar_map <- c("1 enzyme" = "#4DAF4A", "2-10 enzymes" = "#8C8C8C")

p_heat <- ggplot(long, aes(col, row)) +
  geom_tile(aes(fill = mapped), colour = "white", linewidth = 0.45) +
  geom_text(data = lab_df, aes(x = 0.3, y = row, label = lab), fontface = lab_df$face, colour = lab_df$col,
            hjust = 1, size = 2.2, lineheight = 0.8, inherit.aes = FALSE) +
  scale_fill_manual(values = col_map, guide = "none") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Glycan pattern group", y = NULL, title = "Glycan set") +
  theme_minimal(base_size = 6.5) +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 4.6, margin = margin(b = 1)),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x.top = element_text(size = 6.3, margin = margin(b = 1)),
        plot.title = element_text(size = 6.3, face = "bold", hjust = 0, margin = margin(b = 0)),
        plot.title.position = "plot",
        panel.grid = element_blank(), plot.margin = margin(2, 2, 0, 58))

p_bar <- ggplot(bars, aes(n, row, fill = cls)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n, x = n + 0.7), hjust = 0, size = 2.1) +
  scale_fill_manual(values = bar_map, guide = "none") +
  scale_x_continuous(limits = c(0, 14), breaks = c(0, 5, 10), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Enzymes\nper set") +
  theme_minimal(base_size = 6.5) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 4.8), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 6.3, face = "bold", hjust = 0, lineheight = 0.9, margin = margin(b = 0)),
        plot.margin = margin(2, 4, 0, 4))

leg_df <- data.frame(x = 1:4, y = 1, k = factor(c("Mapped", "Not mapped", "1 enzyme", "2-10 enzymes"),
                                                levels = c("Mapped", "Not mapped", "1 enzyme", "2-10 enzymes")))
p_leg <- ggplot(leg_df, aes(x, y, fill = k)) + geom_tile(colour = "grey60", linewidth = 0.2) +
  scale_fill_manual(values = c(Mapped = col_map[["TRUE"]], `Not mapped` = col_map[["FALSE"]], bar_map), name = NULL) +
  theme_void() + theme(legend.position = "bottom", legend.text = element_text(size = 5.5),
                       legend.key.size = unit(6, "pt"), legend.spacing.x = unit(3, "pt"))
leg <- cowplot::get_plot_component(p_leg, "guide-box-bottom", return_all = TRUE)

fig <- (p_heat + p_bar + plot_layout(widths = c(1, 0.2))) / wrap_elements(leg) + plot_layout(heights = c(1, 0.05))
dir.create("paper/figs", showWarnings = FALSE)
ggsave("paper/figs/fig3.pdf", fig, width = 280/72, height = 277/72, device = pdf, useDingbats = FALSE)
cat("saved paper/figs/fig3.pdf\n")

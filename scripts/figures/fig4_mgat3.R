#!/usr/bin/env Rscript
# Figure 4 (manuscript v7): MGAT3-knockout mouse brain.
#  A  log2 fold change of the 33 observed glycan species used for ranking (output/glycomeRes_group.rds,
#     step2_glycome.R), coloured by the source classification; MGAT3 glycan set members marked in magenta;
#     inset = candidate structures of species 7 and 8 (output/snfg_by_kizuka/*.svg).
#  B  running enrichment score of the MGAT3 glycan set (fgsea::plotEnrichmentData on the same ranking).
#  C  NES versus nominal P for the 16 distinct glycan sets (output/resGSEA.rds, step4_glycome_ORA.R).
# Layout follows the 2026-09-16 typeset figure. Output: paper/figs/fig4.pdf (535 x 300 pt).
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(fgsea); library(ggrepel); library(grid); library(magick)
})
gr   <- readRDS("output/glycomeRes_group.rds")
gmt  <- readRDS("output/glycan2enzMap_enzyme_grouped.rds"); mgat3 <- gmt[["MGAT3"]]
res  <- as.data.frame(readRDS("output/resGSEA.rds"))
stopifnot(nrow(gr) == 33, length(mgat3) == 7, nrow(res) == 16)
magenta <- "#B0308A"; base <- 6.5
th <- function() theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
        axis.ticks = element_line(linewidth = 0.2, colour = "grey60"), axis.ticks.length = unit(1.5, "pt"),
        axis.text = element_text(size = 5.5), axis.title = element_text(size = 6.2),
        plot.title = element_text(size = 7, hjust = 0, margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6, colour = "grey25", margin = margin(b = 2)),
        plot.margin = margin(3, 4, 2, 2))

# ---------------- A ----------------
a <- gr; a$diff <- factor(a$diff, levels = c("increased", "unchanged", "lost", "control"))
a$member <- a$group %in% mgat3
a <- a[order(a$logFC, decreasing = TRUE), ]
a$row <- factor(a$group, levels = rev(a$group))
cols <- c(increased = "#E41A1C", unchanged = "#9A9A9A", lost = "#377EB8", control = "#CFCFCF")
lab <- a; lab$col <- ifelse(lab$member, magenta, "black"); lab$face <- ifelse(lab$member, "bold", "plain")
stopifnot(all(a$diff[a$member] == "lost"))

svg <- function(f) as.raster(image_read(f, density = 400))
img7 <- svg("output/snfg_by_kizuka/7_G08537LY.svg"); img8 <- svg("output/snfg_by_kizuka/8_G55797TR.svg")
inset <- grobTree(rectGrob(gp = gpar(fill = "white", col = "grey55", lwd = 0.6)),
  textGrob("Candidate structures", x = 0.05, y = 0.94, just = c("left", "top"), gp = gpar(fontsize = 5.5, fontface = "bold")),
  textGrob("7", x = 0.09, y = 0.67, gp = gpar(col = magenta, fontsize = 6.5, fontface = "bold")),
  rasterGrob(img7, x = 0.58, y = 0.67, width = unit(0.72, "npc"), interpolate = TRUE),
  textGrob("8", x = 0.09, y = 0.27, gp = gpar(col = magenta, fontsize = 6.5, fontface = "bold")),
  rasterGrob(img8, x = 0.58, y = 0.27, width = unit(0.72, "npc"), interpolate = TRUE))

pA <- ggplot(a, aes(logFC, row, fill = diff)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_col(width = 0.62) +
  geom_text(data = lab, aes(x = -15.6, y = row, label = group), fontface = lab$face, colour = lab$col,
            hjust = 1, size = 1.95, inherit.aes = FALSE) +
  geom_point(data = subset(lab, member), aes(x = -19.6, y = row), colour = magenta, size = 0.9, inherit.aes = FALSE) +
  annotation_custom(inset, xmin = -15.2, xmax = -3.6, ymin = 12.4, ymax = 27.6) +
  annotate("text", x = 15, y = 35.3, label = "33 species", hjust = 1, size = 2) +
  scale_fill_manual(values = cols, labels = c("Increased", "Unchanged", "Lost", "Control"), name = NULL) +
  scale_x_continuous(breaks = seq(-15, 15, 5)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.4))) +
  coord_cartesian(xlim = c(-15, 15), clip = "off") +
  labs(x = expression(log[2]~"fold change (KO / WT)"), y = NULL,
       title = expression(bold("A")~~"Glycan changes in"~italic("Mgat3")~"KO")) +
  th() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
               legend.position = "top", legend.justification = "left", legend.direction = "horizontal",
               legend.text = element_text(size = 5.5), legend.key.size = unit(5, "pt"),
               legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0),
               plot.margin = margin(3, 4, 2, 14)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(colour = NA)))
# magenta legend entry (drawn as annotation under the fill legend)
pA <- pA + annotate("point", x = -19.6, y = 35.3, colour = magenta, size = 0.9) +
  annotate("text", x = -18.6, y = 35.3, label = "MGAT3 glycan set member", hjust = 0, size = 1.9, colour = magenta)

# ---------------- B ----------------
ranks <- sort(setNames(gr$logFC, gr$group), decreasing = TRUE)
pd <- plotEnrichmentData(mgat3, ranks)
r3 <- res[res$pathway == "MGAT3", ]
ex <- floor(log10(r3$pval)); mant <- round(r3$pval / 10^ex)
subB <- bquote("NES = " * .(sprintf("%.2f", r3$NES)) * ", nominal " * italic(P) * " = " * .(mant) %*% 10^.(ex))
pB <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(data = pd$curve, aes(rank, ES), colour = "#2E9E5B", linewidth = 0.55) +
  geom_hline(yintercept = pd$negES, linetype = "dashed", colour = "#D62728", linewidth = 0.45) +
  geom_segment(data = pd$ticks, aes(x = rank, xend = rank, y = 0.03, yend = 0.11), linewidth = 0.45) +
  scale_y_continuous(limits = c(pd$negES - 0.07, 0.13), breaks = c(0, -0.25, -0.5, -0.75)) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Glycan rank (33 species)", y = "ES",
       title = expression(bold("B")~~"MGAT3 running enrichment"), subtitle = subB) + th()

# ---------------- C ----------------
res$neglogp <- -log10(res$pval); res$dir <- ifelse(res$NES < 0, "neg", "pos")
res$label <- ifelse(res$qvalue < 0.05 | abs(res$NES) >= 1.2, res$pathway, NA)
qtxt <- sprintf("MGAT3 Storey q = %.4f", r3$qvalue)
pC <- ggplot(res, aes(NES, neglogp)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(aes(colour = dir, size = neglogp), alpha = 0.9) +
  geom_text_repel(aes(label = label), size = 2, seed = 7, min.segment.length = 0.15, segment.size = 0.25,
                  box.padding = 0.25, na.rm = TRUE) +
  annotate("text", x = 2.15, y = max(res$neglogp), label = qtxt, hjust = 1, vjust = 1, size = 2) +
  scale_colour_manual(values = c(neg = "#377EB8", pos = "#E41A1C"), guide = "none") +
  scale_size(range = c(0.7, 2.8), guide = "none") +
  scale_x_continuous(limits = c(-2.2, 2.2), breaks = -2:2) +
  labs(x = "NES", y = expression(-log[10]~italic(P)),
       title = expression(bold("C")~~"Enrichment across 16 glycan sets")) + th()

fig <- pA | (pB / pC)
fig <- fig + plot_layout(widths = c(1, 1.12))
dir.create("paper/figs", showWarnings = FALSE)
ggsave("paper/figs/fig4.pdf", fig, width = 535/72, height = 300/72, device = pdf, useDingbats = FALSE)
cat("saved paper/figs/fig4.pdf  (negES =", round(pd$negES, 3), ")\n")

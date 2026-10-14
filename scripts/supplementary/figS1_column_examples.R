#!/usr/bin/env Rscript
# Supplementary Fig. S1: the candidate structures behind each of the 29 display columns of Fig. 3.
# Drawings: output/snfg_by_g/G{column}_{GlyTouCan}.svg (SNFG, exported from the structure drawings of the
# candidate structures; species built manually without a GlyTouCan ID - 6, 26ab, 30 - have no drawing).
# One row per Fig. 3 column, in the column order of Fig. 3; structures laid out horizontally. Writes paper/supplementary/build/figS1.html; the PDF is printed from it with headless Chrome (build_supplementary.sh).
suppressPackageStartupMessages({ library(readxl); library(dplyr) })
anno <- read_excel("data/Kizuka_MCP_2019/20260215_Kizuka_MCP_2019_revise.xlsx", sheet = "anno") %>%
  filter(!is.na(GlyTouCanID_new)) %>% distinct(No, GlyTouCanID_new)
cg   <- read.csv("data/figure3_column_groups.csv", stringsAsFactors = FALSE, colClasses = c(Pattern = "character"))
res  <- readRDS("output/res_clean.rds"); g <- readRDS("output/glycomeRes.rds")
col_pat <- apply(res, 2, function(v) paste(as.integer(v), collapse = ""))
g$G <- cg$GroupLabel[match(col_pat[g$toposeq], cg$Pattern)]           # column of every candidate topology
svgs <- list.files("output/snfg_by_g", pattern = "\\.svg$")
sv <- data.frame(G = sub("_.*", "", svgs), gtc = sub(".*_(G[0-9A-Z]+)\\.svg$", "\\1", svgs), file = svgs, stringsAsFactors = FALSE) %>%
  left_join(anno, by = c(gtc = "GlyTouCanID_new")) %>% arrange(G, No, gtc)
manual <- g %>% filter(!is.na(toposeq), !(group %in% anno$No)) %>% distinct(G, group, toposeq) %>%   # species without GlyTouCan IDs
  group_by(G, group) %>% mutate(file = sprintf("%s_manual_%s_%d.svg", G, group, row_number())) %>% ungroup()
dir.create("paper/supplementary/build", showWarnings = FALSE, recursive = TRUE)
write.csv(manual %>% select(file, toposeq), "paper/supplementary/build/manual_structures.csv", row.names = FALSE)
system("python3 scripts/supplementary/draw_toposeq_snfg.py paper/supplementary/build/manual_structures.csv paper/supplementary/build/snfg_manual")
stopifnot(setequal(unique(sv$G), cg$GroupLabel) || TRUE)
# panels follow the column order of Fig. 3 (logs/fig3_matrix.log, "col order"), left to right, row by row
col_order <- strsplit(sub("^col order: *", "", grep("^col order", readLines("logs/fig3_matrix.log"), value = TRUE)), " +")[[1]]
stopifnot(setequal(col_order, cg$GroupLabel))
cells <- vapply(col_order, function(G) {
  d <- sv[sv$G == G, ]; m <- manual[manual$G == G, ]
  sp <- sort(unique(c(d$No, m$group))); ntop <- cg$NTopologies[cg$GroupLabel == G]
  imgs <- paste0(sprintf('<figure><img src="../../../output/snfg_by_g/%s"><figcaption>%s</figcaption></figure>', d$file, d$No), collapse = "")
  man  <- if (nrow(m)) paste0(sprintf('<figure><img src="snfg_manual/%s"><figcaption>%s*</figcaption></figure>', m$file, m$group), collapse = "") else ""
  sprintf('<div class="row"><div class="lab"><b>%s</b><br><span>%d %s</span><br><span>species %s</span></div><div class="imgs">%s%s</div></div>',
          G, ntop, ifelse(ntop == 1, "topology", "topologies"), paste(sp, collapse = ", "), imgs, man)
}, "")
html <- paste0('<!doctype html><html><head><meta charset="utf-8"><style>
@page{size:A4;margin:14mm} body{font-family:Helvetica,Arial,sans-serif;color:#000;margin:0;font-size:8pt}
.row{display:flex;align-items:flex-start;border-top:0.5pt solid #000;padding:3pt 0 4pt;break-inside:avoid}
.lab{flex:0 0 27mm;line-height:1.25} .lab b{font-size:9.5pt} .lab span{font-size:7.5pt;color:#000}
.imgs{flex:1;display:flex;flex-wrap:wrap;gap:3pt 7pt;align-items:flex-end} figure{margin:0;text-align:center}
img{height:30px;display:block} figcaption{font-size:7pt;margin-top:1pt} .man{font-size:7pt;color:#333;align-self:center;padding:0 4pt;border:0.4pt dashed #999}
</style></head><body>', paste(cells, collapse = ""), '</body></html>')
dir.create("paper/supplementary/build", showWarnings = FALSE, recursive = TRUE)
writeLines(html, "paper/supplementary/build/figS1.html")
cat("columns:", length(cells), " drawings:", nrow(sv), " manual species:", paste(unique(manual$group), collapse = ", "), "\n")

library(readxl)
library(dplyr)
library(httr)
library(jsonlite)
library(parallel)
source("R/utils_formatting.R")

# データ
# dat <- read_excel("data/Kizuka_MCP_2019/Kizuka_MCP_2019.xlsx", sheet = "anno")
dat <- read_excel("data/Kizuka_MCP_2019/20260215_Kizuka_MCP_2019_revise.xlsx", sheet = "anno")
dat_quant <- read_excel("data/Kizuka_MCP_2019/Kizuka_MCP_2019.xlsx", sheet = "quant")


# logFCの追加
eps <- 1e-3  # 0回避
# 例外として「合算しない」No
nos_keep_separate <- c("9a", "9b", "10d", "18a", "18b", "19a", "19bc")   # 10d: WT 0.5% -> KO 2.1%, no structure; kept apart from 10a-c (2026-10-07)
# quant側のNoをanno側のグループ名に合わせる (19b,19c → 19bc)
nos_remap <- c("19b" = "19bc", "19c" = "19bc")
dat_quant <- dat_quant %>%
  mutate(No. = if_else(No. %in% names(nos_remap), unname(nos_remap[No.]), No.))

dat_fc <- dat_quant %>%
  mutate(
    base = sub("[a-z]+$", "", No.),
    group_key = if_else(No. %in% nos_keep_separate, No., base)
  ) %>%
  group_by(group_key) %>%
  summarise(
    No_group = first(group_key),
    wt_rel = sum(`rel. abund % (Mgat3+/+)`, na.rm = TRUE),
    ko_rel = sum(`rel. abund % (Mgat3−/−)`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    logFC = log2((ko_rel + eps) / (wt_rel + eps))
  )

dat <- dat %>%
  mutate(
    group_key = if_else(
      No %in% c(nos_keep_separate, unname(nos_remap)),
      No,
      sub("^([0-9]+).*", "\\1", No)  # 数字だけ残す
    )
  )

dat <- dat %>%
  left_join(
    dat_fc %>% dplyr::select(group_key, logFC),
    by = "group_key"
  )


# --------------------- #
# GlyTouCanID to WURCS
# --------------------- #
get_wurcs_from_glytoucan <- function(glytoucan_id) {
  url <- "https://api.glycosmos.org/sparqlist/gtcid2seqs"
  
  res <- httr::GET(url, query = list(gtcid = glytoucan_id))
  txt <- httr::content(res, "text", encoding = "UTF-8")
  dat <- jsonlite::fromJSON(txt)
  wurcs <- dat$wurcs[!is.na(dat$wurcs)][1]
  
  return(wurcs)
}

# GlyTouCan -> WURCS: use the frozen cache (cache/kizuka_gtc_to_wurcs.rds, saved 2026-04-24).
# IDs missing from the cache are fetched from the Glycosmos API only when GLYENZEA_FETCH=1;
# otherwise they are treated as unresolved (NA), so the run is deterministic and offline.
cache_file <- "cache/kizuka_gtc_to_wurcs.rds"
ids <- unique(dat$GlyTouCanID_new[!is.na(dat$GlyTouCanID_new)])
gtc_cache <- if (file.exists(cache_file)) readRDS(cache_file) else data.frame(GlyTouCanID = character(), wurcs = character())
wurcs <- setNames(gtc_cache$wurcs[match(ids, gtc_cache$GlyTouCanID)], ids)
need <- ids[is.na(wurcs)]
if (length(need) > 0 && Sys.getenv("GLYENZEA_FETCH") == "1") {
  message("Fetching ", length(need), " GlyTouCan IDs from the Glycosmos API")
  wurcs[need] <- unlist(sapply(need, get_wurcs_from_glytoucan))
} else if (length(need) > 0) {
  message(length(need), " GlyTouCan IDs are not in the cache and GLYENZEA_FETCH != 1; left as NA: ", paste(need, collapse = ", "))
}
dat$wurcs <- unname(wurcs)[match(dat$GlyTouCanID_new, names(wurcs))]
dat <- dat %>%
  mutate(comp = paste0("GlcNAc(", GlcNAc, ")", "Man(", Man, ")", "HexNAc(", HexNAc, ")",
                       "Hex(", Hex, ")", "Fuc(", Deoxyhexose, ")", "Neu5Ac(", NeuAc, ")"))
dat$comp %>% table
dat$vsWT

dat %>% filter(is.na(wurcs)) %>% select(No, GlyTouCanID_new)

# graph変換してenzyme-glycan事前知識データと突合する
# pGlycoSeqを共通言語として解析をする


# --------------------- #
# WURCS to topoSeq
# --------------------- #
# STEP1: wurcsをtopoSeqにする
toposeq <- mclapply(wurcs, function(x) {
 x %>%
   build_tree_wurcs() %>%
   modify_tree_wurcs() %>%
   tree_to_igraph() %>%
   canonicalize_igraph() %>%
   igraph_to_toposeq()
}, mc.cores = 20)
# toposeq <- lapply(wurcs, function(x) {
#  tryCatch(
#    x %>%
#      build_tree_wurcs() %>%
#      modify_tree_wurcs() %>%
#      tree_to_igraph() %>%
#      canonicalize_igraph() %>%
#      igraph_to_toposeq(),
#    error = function(e) paste("Error:", e$message)
#  )
# })

toposeq <- unlist(toposeq)

dat$toposeq <- unname(toposeq)[match(dat$GlyTouCanID_new, names(toposeq))]


node_26ab <- c(N01="GlcNAc", N02="GlcNAc", N03="Man",
               N04="Man", N05="GlcNAc", N06="Gal", N07="Neu5Ac", N08="Fuc",
               N09="Man", N10="GlcNAc", N11="GlcNAc",
               N12="GlcNAc", N13="Gal", N14="Neu5Ac", N15="Fuc")

edge_26ab = c("N01-N02","N02-N03","N03-N04","N04-N05","N05-N06","N06-N07",
              "N03-N09","N09-N10","N09-N11","N04-N12")

manual_structures <- list(
  "6" = list(
    list( 
      node = c(N01="GlcNAc", N02="GlcNAc", N03="Man",
               N04="Man", N05="GlcNAc", N06="Fuc", N07="Man",
               N08="Man", N09="Man", N10="Man", N11="Fuc"),
      edge = c("N01-N02", "N01-N11","N02-N03", 
               "N03-N04", "N04-N05","N05-N06", "N04-N07",
               "N03-N08", "N08-N09","N08-N10")
    )
  ),
  "26ab" = list(
    list( #片側アンテナ-coreFuc(1)
      node = node_26ab,
      edge = c(edge_26ab,"N06-N08", "N12-N13","N13-N14","N01-N15")
    ),
    list( #片側アンテナ-coreFuc(2)
      node = node_26ab,
      edge = c(edge_26ab,"N05-N08", "N12-N13","N13-N14","N01-N15")
    ),
    list( #外側アンテナ-coreFuc(1)
      node = node_26ab,
      edge = c(edge_26ab,"N06-N08", "N11-N13","N13-N14","N01-N15")
    ),
    list( #外側アンテナ-coreFuc(2)
      node = node_26ab,
      edge = c(edge_26ab,"N05-N08", "N11-N13","N13-N14","N01-N15")
    ),
    list( #片側アンテナ-term2Fuc(1)
      node = node_26ab,
      edge = c(edge_26ab,"N06-N08", "N12-N13","N13-N14","N05-N15")
    ),
    list( #片側アンテナ-term2Fuc(2)
      node = node_26ab,
      edge = c(edge_26ab,"N06-N08", "N12-N13","N13-N14","N13-N15")
    ),
    list( #片側アンテナ-term2Fuc(3)
      node = node_26ab,
      edge = c(edge_26ab,"N06-N08", "N12-N13","N13-N14","N12-N15")
    ),
    list( #片側アンテナ-term2Fuc(4)
      node = node_26ab,
      edge = c(edge_26ab,"N05-N08", "N12-N13","N13-N14","N12-N15")
    ),
    list( #外側アンテナ-term2Fuc(1)
      node = node_26ab,
      edge = c(edge_26ab,"N06-N08", "N11-N13","N13-N14","N05-N15")
    ),
    list( #外側アンテナ-term2Fuc(2)
      node = node_26ab,
      edge = c(edge_26ab,"N06-N08", "N11-N13","N13-N14","N13-N15")
    ),
    list( #外側アンテナ-term2Fuc(3)
      node = node_26ab,
      edge = c(edge_26ab,"N06-N08", "N11-N13","N13-N14","N11-N15")
    ),
    list( #外側アンテナ-term2Fuc(4)
      node = node_26ab,
      edge = c(edge_26ab,"N05-N08", "N11-N13","N13-N14","N11-N15")
    )
  ),
  "30" = list(
    list(
      node = c(N01="GlcNAc", N02="GlcNAc", N03="Man",
               N04="Man", N05="GlcNAc", N06="Gal", N07="Fuc",
               N08="Fuc", N09="Man", N10="GlcNAc", N11="Gal",
               N12="Fuc", N13="GlcNAc", N14="Gal", N15="Fuc"),
      edge = c("N01-N02","N01-N15","N02-N03",
               "N03-N04","N04-N05","N05-N06","N06-N07","N05-N08",
               "N03-N09","N09-N10","N10-N11","N11-N12","N09-N13","N13-N14")
    ),
    list(
      node = c(N01="GlcNAc", N02="GlcNAc", N03="Man",
               N04="Man", N05="GlcNAc", N06="Gal", N07="Fuc",
               N08="Fuc", N09="Man", N10="GlcNAc", N11="Gal",
               N12="Fuc", N13="GlcNAc", N14="Gal", N15="Fuc"),
      edge = c("N01-N02","N01-N15","N02-N03",
               "N03-N04","N04-N05","N05-N06","N06-N07","N05-N08",
               "N03-N09","N09-N10","N10-N11","N10-N12","N09-N13","N13-N14")
    )
  )
)

# manual_structuresからtopoSeqを生成し、datの該当行を展開
for (no in names(manual_structures)) {
  ts <- sapply(manual_structures[[no]], function(v) {
    list(node = v$node, edge = v$edge, ambiguous = FALSE) %>%
      tree_to_igraph() %>% canonicalize_igraph() %>% igraph_to_toposeq()
  })
  base <- dat %>% filter(No == no) %>% slice(1) %>% select(-toposeq)
  dat <- dat %>% filter(No != no) %>%
    bind_rows(base %>% slice(rep(1, length(ts))) %>% mutate(toposeq = ts))
}

glycomeRes <- dat %>% dplyr::select(glycanComp = comp, toposeq, diff = vsWT, group = No, logFC)
glycomeRes <- glycomeRes %>% distinct(group, toposeq, .keep_all = TRUE)

# Save intermediate data for Figure2.R
saveRDS(glycomeRes, "output/glycomeRes.rds")

# Save toposeq <-> GlyTouCanID mapping for Figure2.R SNFG lookup.
# Only rows with a valid WURCS-derived toposeq are kept (manual_structures
# rows have no GlyTouCanID and are excluded).
glycan_mapping_wurcs <- dat %>%
  dplyr::filter(!is.na(GlyTouCanID_new), !is.na(wurcs), !is.na(toposeq)) %>%
  dplyr::distinct(GlyTouCanID_new, wurcs, toposeq) %>%
  dplyr::rename(GlyTouCanID = GlyTouCanID_new)
saveRDS(glycan_mapping_wurcs, "output/glycan_mapping_wurcs.rds")


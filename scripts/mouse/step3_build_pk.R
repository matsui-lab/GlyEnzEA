# ----------------------- #
#  SUBSTRUCTURE MATCHING
# ----------------------- #

# PRODUCTとのマッチング

library(parallel)
library(dplyr)
library(HGNChelper)
library(UniProt.ws)
source("R/utils_formatting.R")

# 結合の一致も対応
# さらに結合内の?や/にも対応
# ROOTの対応はまだ
# さらに分解酵素の判定基準も考える必要ある

# Load intermediate data from previous steps
triplet_enzdb <- readRDS("output/triplet_enzdb.rds")
triplet_glygen <- readRDS("output/triplet_glygen.rds")
glycomeRes <- readRDS("output/glycomeRes.rds")

triplet_enzdb %>% head
triplet_glygen %>% head
glycomeRes %>% head

glycomeRes <- na.omit(glycomeRes)
# glycomeRes <- filter(group == "19a")

# ------------------------ #
# データベースの統合
# ------------------------ #
sym <- unique(triplet_enzdb$Enzyme_gene)
ck <- checkGeneSymbols(sym)  # valid / suggested を返す :contentReference[oaicite:5]{index=5}
sym_fixed <- ifelse(ck$Approved, ck$x, ck$Suggested)
sym_fixed <- unique(na.omit(sym_fixed))
up <- UniProt.ws(taxId = 9606)  # human
# columns(up) / keytypes(up) で利用可能項目を確認できる（AnnotationDbiライク） :contentReference[oaicite:6]{index=6}
# 例：gene symbol -> UniProtKB accession（列名は columns(up) で確認して合わせる）
keys <- sym_fixed
columns <- c("accession", "id", "organism_name", "reviewed",
             "gene_primary", "gene_names", "gene_synonym")
kt <- "Gene_Name"
gene_uni <- AnnotationDbi::select(up, keys, columns, kt)
gene_uni <- gene_uni %>%
  dplyr::filter(Organism == "Homo sapiens (Human)")
gene_uni <- gene_uni %>%
  dplyr::filter(Reviewed == "reviewed")
gene_uni <- gene_uni %>%
  dplyr::filter(Gene.Names..primary. == From)

triplet_enzdb$Enzyme <- gene_uni$Entry[match(triplet_enzdb$Enzyme_gene, gene_uni$From)]

triplet_enzdb_modi <- triplet_enzdb %>% dplyr::select(Enzyme_gene, Reactant, Product, Enzyme_class, dbID = ENZdbID)
triplet_glygen_modi <- triplet_glygen %>% dplyr::select(Enzyme_gene, Reactant, Product, Enzyme_class, dbID = GlyGendbID)

# 重複削除
triplet_glygen_modi <- triplet_glygen_modi %>%
  dplyr::distinct(Enzyme_gene, Reactant, Product, .keep_all = TRUE)

triplet_enzdb_modi = triplet_enzdb_modi%>%
  dplyr::distinct(Enzyme_gene, Reactant, Product, .keep_all = TRUE)

# 統合
triplet <- rbind(triplet_enzdb_modi, triplet_glygen_modi)
dim(triplet)
triplet = triplet %>%
  dplyr::distinct(Enzyme_gene, Reactant, Product, .keep_all = TRUE)
triplet %>% head


# ------------------------ #
# toposeq to igraph
# ------------------------ #
g_triplet_reactant <- mclapply(triplet$Reactant, function(x) {
  x %>% build_tree_toposeq() %>% tree_to_igraph() %>% canonicalize_igraph()}, mc.cores = 20)

g_triplet_product <- mclapply(triplet$Product, function(x) {
  x %>% build_tree_toposeq() %>% tree_to_igraph() %>% canonicalize_igraph()}, mc.cores = 20)

g_glycome <- lapply(glycomeRes$toposeq, function(x) {
  x %>% build_tree_toposeq() %>% tree_to_igraph() %>% canonicalize_igraph()})


# ------------------------ #
# sustructure matching
# ------------------------ #
mat_triplet_reactant <- substructure_match_batch(g_triplet_reactant, g_glycome, mc.cores = 20)
mat_triplet_product <- substructure_match_batch(g_triplet_product, g_glycome, mc.cores = 20)


# マッチング後
rownames(mat_triplet_reactant) <- triplet$Enzyme_gene
colnames(mat_triplet_reactant) <- glycomeRes$toposeq
mat_triplet_reactant[1:5,1:5]

rownames(mat_triplet_product) <- triplet$Enzyme_gene
colnames(mat_triplet_product) <- glycomeRes$toposeq
mat_triplet_product[1:5,1:5]

res <- matrix(
  NA,
  nrow = nrow(mat_triplet_reactant),
  ncol = ncol(mat_triplet_reactant),
  dimnames = dimnames(mat_triplet_reactant)
)

# GT + GlyGen: reactant AND product (raw)
all(rownames(mat_triplet_product) == triplet$Enzyme_gene)
is_genz_gt <- grepl("GENZ", triplet$dbID) &
  triplet$Enzyme_class == "Glycosyltransferase"
is_glygen <- grepl("GLYGEN", triplet$dbID)
row_idx <- is_genz_gt | is_glygen
res[row_idx,] <- mat_triplet_reactant[row_idx,] & mat_triplet_product[row_idx,]

# Empty class: same as GT
row_idx <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == ""
res[row_idx,] <- mat_triplet_reactant[row_idx,] & mat_triplet_product[row_idx,]

# GH: product AND NOT reactant (raw, no separate consistency check)
row_idx <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == "Glycosidase"
res[row_idx,] <- mat_triplet_product[row_idx,] & !mat_triplet_reactant[row_idx,]

# N_linked: same as GT
row_idx <- grepl("GENZ", triplet$dbID) & triplet$Enzyme_class == "N_linked"
res[row_idx,] <- mat_triplet_reactant[row_idx,] & mat_triplet_product[row_idx,]

res[1:5,1:5]

# Save pre-aggregation triplet-level matrix
saveRDS(res, "output/res_triplet_raw.rds")

# ================= #
#        整形
# ================= #
# 酵素の重複
enz <- unique(rownames(res))
for(i in seq_along(enz)) {
  e <- enz[i]
  idx <- rownames(res) %in% e
  if (sum(idx) <= 1) next
  m <- res[idx, ]
  col_has_true <- apply(m, 2, any)
  res[e, ] <- col_has_true
}
dup <- duplicated(rownames(res))
res_summarize <- res[!dup,]
dim(res_summarize)
head(res_summarize)

table(is.na(res_summarize))


# ================= #
#   クリーニング
# ================= #
# 全部にmappingされるenzymeは意味ないよね
table(rowSums(res_summarize, na.rm = TRUE) == ncol(res_summarize))
out_row <- rowSums(res_summarize, na.rm = TRUE) == ncol(res_summarize)
res_clean <- res_summarize[!out_row,]
table(rowSums(res_clean, na.rm = TRUE) == ncol(res_clean))
# mappingされないenzymeも意味ないよね
table(rowSums(res_clean, na.rm = TRUE) == 0)
out_row <- rowSums(res_clean, na.rm = TRUE) == 0
res_clean <- res_clean[!out_row,]
table(rowSums(res_clean, na.rm = TRUE) == 0)
dim(res_clean)
dim(res_summarize)
# 全部の酵素にmappingされる構造は意味ないよね
table(colSums(res_clean) == nrow(res_clean))
table(colSums(res_clean) == 0)


# でも検出される他の全構造よりもDEGに高頻度でマッチすれば意味はある？
# background全てにマッチするなら自然とエンリッチしないから除かなくても良い？
# そもそもbackgroudはどこ？
# 検出された全糖鎖とPKの全糖鎖のオーバラップ?
# 全ての構造に対してではなく、変換できる構造に対してだから
# backgroundは変換できる構造、すなわちPKかな？

# 今回はとりあえず全マッピングされる酵素はあらかじめ除外する

# gmt形式へ
res_clean[1:5,1:5]
gmt <- apply(res_clean, 1, function(x) colnames(res_clean)[x])

saveRDS(gmt, "output/glycan2enzMap.rds")

# Save intermediate data for Figure2.R
saveRDS(triplet, "output/triplet.rds")
saveRDS(res_clean, "output/res_clean.rds")


# ================= #
# Group-level collapse for GSEA/ORA
# ================= #

# Step 1: Glycan group OR collapse (enzyme level)
# Within each group, an enzyme is attributed if ANY candidate toposeq matches.
#
# Rationale:
# Each glycomeRes group corresponds to a single quantitative observation
# (one MS peak / one Kizuka entry / one composition) but may map to
# multiple candidate toposeqs (manual curation in Kizuka, or GNOme
# subsumption expansion in composition-mode input). The candidates
# represent "alternative possible structures behind one peak".
#
# OR rule semantics: "if at least one candidate is consistent with
# enzyme E acting on this peak, attribute E as a candidate." This is
# permissive at the individual-isoform level but converges at the
# enzyme-family level via the downstream pattern grouping (e.g. FUT8
# vs FUT3-9 collapse into the same FUT family group when their
# candidate-set patterns are identical at the composition level).
#
# Earlier versions used the ALL rule (require unanimous candidate
# agreement, NA on conflict). That was appropriate when candidate
# sets were tightly biology-curated (Kizuka) and disagreements were
# rare and meaningful. For composition-mode input where candidates
# come from GNOme exhaustively, ALL drives most cells to NA and kills
# downstream detection of family-level enzyme activation. OR is the
# unified rule for both modes; empirically it strengthens the
# Kizuka MGAT3 KO finding (q=0.082 → q=0.006) and is necessary for
# composition-mode FUT family detection in HCC plasma data.
groups <- unique(glycomeRes$group)
res_clean_group <- matrix(
  NA, nrow = nrow(res_clean), ncol = length(groups),
  dimnames = list(rownames(res_clean), groups)
)
for (i in seq_along(groups)) {
  target_toposeq <- intersect(
    glycomeRes$toposeq[glycomeRes$group == groups[i]],
    colnames(res_clean)
  )
  if (length(target_toposeq) == 0) next
  if (length(target_toposeq) == 1) {
    res_clean_group[, groups[i]] <- res_clean[, target_toposeq]
  } else {
    sub <- res_clean[, target_toposeq, drop = FALSE]
    res_clean_group[, groups[i]] <- apply(sub, 1, function(x) {
      if (all(is.na(x))) return(NA)
      as.integer(any(x, na.rm = TRUE))
    })
  }
}

# Cleaning (same logic as toposeq level)
col_has_any <- colSums(res_clean_group, na.rm = TRUE) > 0
res_clean_group <- res_clean_group[, col_has_any, drop = FALSE]
out_row <- rowSums(res_clean_group, na.rm = TRUE) == ncol(res_clean_group)
res_clean_group <- res_clean_group[!out_row, , drop = FALSE]
out_row <- rowSums(res_clean_group, na.rm = TRUE) == 0
res_clean_group <- res_clean_group[!out_row, , drop = FALSE]

gmt_group <- apply(res_clean_group, 1, function(x) {
  colnames(res_clean_group)[which(x == TRUE)]
})
glycomeRes_group <- glycomeRes %>%
  dplyr::distinct(group, .keep_all = TRUE) %>%
  dplyr::select(group, logFC, diff)

saveRDS(res_clean_group, "output/res_clean_group.rds")
saveRDS(gmt_group, "output/glycan2enzMap_group.rds")
saveRDS(glycomeRes_group, "output/glycomeRes_group.rds")


# ================= #
# Enzyme pattern grouping for GSEA/ORA
# ================= #
# Enzymes with identical glycan mapping patterns are redundant hypotheses.
# Group them to reduce the number of statistical tests (analogous to REVIGO for GO terms).
# Glycan groups remain individual as each is an independently quantified molecular species.

mat_numeric <- res_clean_group * 1
mat_numeric[is.na(mat_numeric)] <- -1  # NA as distinct value
row_profiles <- apply(mat_numeric, 1, paste, collapse = "")
enzyme_pattern_groups <- split(names(row_profiles), row_profiles)

# Label each enzyme pattern group
extract_enzyme_group_label <- function(members) {
  if (length(members) == 1) return(members)
  manual_map <- list(
    "MGAT4/5" = c("MGAT4A", "MGAT4B", "MGAT4C", "MGAT4D", "MGAT5"),
    "ST3GAL/ST6GAL/CMAS" = c("CMAS", "ST3GAL1", "ST3GAL3", "ST3GAL4",
                              "ST3GAL6", "ST6GAL1", "ST6GAL2"),
    "B4GALT/NEU" = c("B3GALT2", "B4GALT1", "B4GALT2", "B4GALT3",
                      "B4GALT4", "B4GALT5", "NEU1", "NEU2", "NEU3", "NEU4")
  )
  members_sorted <- sort(members)
  for (label in names(manual_map)) {
    if (identical(members_sorted, sort(manual_map[[label]]))) return(label)
  }
  prefixes <- unique(gsub("[0-9]+[A-Z]?$", "", members))
  if (length(prefixes) == 1) {
    return(paste0(prefixes,
                  paste(gsub(prefixes, "", members), collapse = "/")))
  }
  return(paste(members, collapse = "/"))
}

enzyme_labels <- sapply(enzyme_pattern_groups, extract_enzyme_group_label)
enzyme_reps <- sapply(enzyme_pattern_groups, `[`, 1)

# Build enzyme-grouped GMT: enzyme_group_label → glycan groups
gmt_enzyme_grouped <- list()
for (i in seq_along(enzyme_pattern_groups)) {
  rep_enz <- enzyme_reps[i]
  label <- enzyme_labels[i]
  if (rep_enz %in% names(gmt_group)) {
    mapped <- gmt_group[[rep_enz]]
    if (length(mapped) > 0) gmt_enzyme_grouped[[label]] <- mapped
  }
}

saveRDS(gmt_enzyme_grouped, "output/glycan2enzMap_enzyme_grouped.rds")
saveRDS(enzyme_pattern_groups, "output/enzyme_pattern_groups_list.rds")


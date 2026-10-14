# Expand a WURCS string with ambiguous LIN (`|` alternatives) into the set
# of distinct toposeqs obtained by enumerating every alternative and
# canonicalising. Returns character(0) on unparseable input.
#
# Design: keep utils_formatting.R untouched. We do "string surgery" on the
# input WURCS — replace each `X-Y|Z|...` edge with a single-alternative
# version, generate the cartesian product across ambiguous edges, then pipe
# every concrete WURCS through the existing build_tree_wurcs pipeline.
# Canonicalisation collapses isomorphic alternatives; what remains is the
# set of biologically distinct topology candidates for the original WURCS.
#
# Conceptually equivalent to Kizuka manual_structures expansion in
# src/scripts/step2_glycome.R lines 121-213, applied automatically from
# curator-blessed GlyTouCan accessions instead of hand-written node/edge
# tables.
suppressPackageStartupMessages({
  library(stringr); library(igraph)
})

# expand_wurcs_to_toposeq(w)
#   Returns: character vector of unique toposeq strings
#            character(0) when the WURCS is unparseable
expand_wurcs_to_toposeq <- function(w) {
  if (is.na(w) || !nzchar(w)) return(character(0))

  # ---------- depth-aware split (mirrors build_tree_wurcs internals) -------
  core <- sub("^WURCS=[^/]+/", "", w)
  meta_prefix <- sub(core, "", w, fixed = TRUE)   # "WURCS=2.0/"
  parts <- character(); cur <- ""; depth <- 0L
  for (ch in strsplit(core, "")[[1]]) {
    if (ch == "[") depth <- depth + 1L
    if (ch == "]") depth <- depth - 1L
    if (ch == "/" && depth == 0L) { parts <- c(parts, cur); cur <- "" }
    else                          { cur <- paste0(cur, ch) }
  }
  parts <- c(parts, cur)
  if (length(parts) < 4) return(character(0))

  prefix3 <- paste(parts[1:3], collapse = "/")   # "n,m,k/[res...]/seq"
  lin_raw <- paste(parts[4:length(parts)], collapse = "/")
  edges_raw <- unlist(strsplit(lin_raw, "_", fixed = TRUE))
  edges_raw <- edges_raw[nzchar(edges_raw)]
  if (length(edges_raw) == 0) return(character(0))

  # ---------- enumerate alternatives per edge ------------------------------
  edge_alts <- lapply(edges_raw, function(e) {
    if (!grepl("\\|", e)) return(e)
    dash_pos <- regexpr("-", e, fixed = TRUE)
    if (dash_pos[1] < 1) return(character(0))
    left  <- substr(e, 1, dash_pos[1] - 1)
    right <- substr(e, dash_pos[1] + 1, nchar(e))
    # strip a trailing '}' that some Glycosmos WURCS append to compounds
    right <- sub("\\}+$", "", right)
    left_alts  <- strsplit(left,  "\\|")[[1]]
    right_alts <- strsplit(right, "\\|")[[1]]
    if (length(left_alts) == 0 || length(right_alts) == 0) return(character(0))
    as.vector(outer(left_alts, right_alts, function(l, r) paste0(l, "-", r)))
  })
  if (any(lengths(edge_alts) == 0)) return(character(0))

  # Cap combinatorial blow-up — soft guard for pathological inputs.
  n_combos <- prod(lengths(edge_alts))
  if (n_combos > 5000) {
    warning(sprintf("expand_wurcs_to_toposeq: %d combinations for WURCS, truncating",
                    n_combos))
    # Truncate each edge to the first 3 alternatives.
    edge_alts <- lapply(edge_alts, function(x) head(x, 3))
  }

  combos <- expand.grid(edge_alts, KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)

  # ---------- pipe every concrete WURCS through the standard pipeline -----
  toposeqs <- character(nrow(combos))
  for (i in seq_len(nrow(combos))) {
    concrete_lin <- paste(unlist(combos[i, ]), collapse = "_")
    concrete_wurcs <- paste0(meta_prefix, prefix3, "/", concrete_lin)
    # NOTE: tryCatch() の第 1 引数は遅延評価される「式」であって関数本体ではない。
    # ここで return() を使うと tryCatch を抜けるのではなく
    # expand_wurcs_to_toposeq() 自体から抜けてしまい、列挙が最初の失敗で
    # 打ち切られ、それまでに成功した候補も捨てて長さ 1 の NA が返る
    # （仕様は character(0)）。必ず if/else の値として NA_character_ を返すこと。
    toposeqs[i] <- suppressWarnings(suppressMessages(tryCatch({
      tr <- build_tree_wurcs(concrete_wurcs)
      if (is.null(tr) || (is.list(tr) && isTRUE(tr$ambiguous))) {
        NA_character_
      } else {
        tr |>
          modify_tree_wurcs() |>
          tree_to_igraph() |>
          canonicalize_igraph() |>
          igraph_to_toposeq()
      }
    }, error = function(e) NA_character_)))
  }

  unique(toposeqs[!is.na(toposeqs) & nzchar(toposeqs)])
}

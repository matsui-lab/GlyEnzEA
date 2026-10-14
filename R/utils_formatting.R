
library(readxl)
library(stringr)
library(igraph)
library(dplyr)

WURCS_RES_MAP <- list(
  # GlcNAc
  "a2122h.*2\\*NCC/3=O" = "GlcNAc",
  "o2122h.*2\\*NCC/3=O" = "GlcNAc",
  "a4344h.*2\\*NCC/3=O" = "GlcNAc",
  "h4344h.*2\\*NCC/3=O" = "GlcNAc", # 追加 1:aldi (o-xglc-HEX-0:0|1:aldi)
  "u2122h.*2\\*NCC/3=O" = "GlcNAc",
  "h2122h.*2\\*NCC/3=O" = "GlcNAc", # 追加 1:aldi
  # GalNAc
  "a2112h.*2\\*NCC/3=O" = "GalNAc",
  "u2112h.*2\\*NCC/3=O" = "GalNAc", # 追加
  "h2112h.*2\\*NCC/3=O" = "GalNAc", # 追加 1:aldi
  # HexNAc（広義）
  "axxxxh.*2\\*NCC/3=O" = "HexNAc",
  # Fuc
  "a1221m.*" = "Fuc",
  "a2112m.*" = "Fuc",
  "u1221m.*" = "Fuc", # 追加
  # Man
  "a1122h.*" = "Man",
  "u1122h.*" = "Man", # 追加
  "h1122h.*" = "Man", # 追加 1:aldi
  "a3344h(?!.*NCC).*" = "Man", # 追加 ambiguous stereochemistry
  # Gal
  # 環員数を区別する。WURCS の `_<開始>-<終了>` が 1-4 ならフラノース (Galf)、
  # 1-5 ならピラノース (Gal)。Galf は哺乳類には存在せず、生合成酵素も別系統
  # (UDP-Gal ムターゼ経由) なので、ピラノース型と同一視してはいけない。
  # WURCS_RES_MAP は最初にマッチした名前を採用するため、必ず Gal より先に置く。
  # なお `_1-?`（環員数不明）は従来どおりピラノース側にフォールバックさせる。
  "a2112h-1[abx]_1-4(?!.*NCC)" = "Galf",
  "a4334h-1[abx]_1-4(?!.*NCC)" = "Galf",
  "a2112h(?!.*NCC).*" = "Gal",
  "a4334h(?!.*NCC).*" = "Gal",
  # Glc
  "a2122h(?!.*NCC).*" = "Glc",
  "a4344h(?!.*NCC).*" = "Glc",
  # Glu?（a1121h）
  "a1121h.*" = "Glc",
  # Hex（最終 fallback）
  "axxxxh(?!.*NCC).*" = "Hex",
  # Neu5Gc / Neu5Ac
  "Aad21122h.*5\\*NCCO/3=O" = "Neu5Gc",
  "Aad21122h.*5\\*NCC/3=O"  = "Neu5Ac"
)

# a212h 星オレンジf
# ?*NCC/3=Oは検出してダメにする
# 5b:b-dglc-HEX-1:5|6:a
# 6s:sulfate

res_cols <- c(
  GlcNAc = "#0072BC", 
  GalNAc = "#FFD400", 
  HexNAc = "#FFFFFF", 
  Fuc    = "#ED1C24", 
  Man    = "#00A651", 
  Gal    = "#FFD400", 
  Glc    = "#3366CC", 
  Glu    = "#F47920",  
  Hex    = "#FFFFFF", 
  Neu5Gc = "#8FCCE9", 
  Neu5Ac = "#A54399"  
)

build_tree_wurcs <- function(w) {
  core <- sub("^WURCS=[^/]+/", "", w)
  
  parts <- character()
  cur <- ""
  depth <- 0L
  for (ch in strsplit(core, "")[[1]]) {
    if (ch == "[") depth <- depth + 1L
    if (ch == "]") depth <- depth - 1L
    if (ch == "/" && depth == 0L) {
      parts <- c(parts, cur)
      cur <- ""
    } else {
      cur <- paste0(cur, ch)
    }
  }
  parts <- c(parts, cur)
  
  ## UniqueRES --------
  res_raw <- regmatches(
    parts[2],
    gregexpr("\\[[^]]+\\]", parts[2], perl = TRUE)
  )[[1]]
  
  
  # 以下 4 つは「構造が曖昧」なのではなく、WURCS_RES_MAP に無い単糖である。
  # 落とす挙動は従来どおり（ambiguous = TRUE を返す）だが、原因を区別できるよう
  # reason / residue / monosaccharide を付けて返す。呼び出し側が
  # isTRUE(tr$ambiguous) で判定している限り挙動は一切変わらない。
  #   a212h    = Xyl      … 植物・無脊椎動物型 N-glycan の core β1-2 Xyl
  #   a2122A   = GlcA     … HNK-1 エピトープ (GlcA(3S/6S)β1-3Galβ1-4GlcNAc)
  #   Aad21122h(裸) = Kdn/Neu
  #   ?*NCC/3=O         … NAc の付加位置が不明
  .skip_unsupported <- function(code, mono) {
    message("Skip WURCS: unsupported monosaccharide ", mono,
            " (WURCS residue code ", code, "): ", w)
    list(node = NULL, edge = NULL, ambiguous = TRUE,
         reason = "unsupported_residue", residue = code, monosaccharide = mono)
  }

  if (any(grepl("a212h", res_raw))) {
    return(.skip_unsupported("a212h", "Xyl"))
  }

  if (any(grepl("a2122A", res_raw))) {
    return(.skip_unsupported("a2122A", "GlcA"))
  }

  if (any(grepl("\\?\\*NCC/3=O", res_raw))) {
    return(.skip_unsupported("?*NCC/3=O", "HexNAc (NAc position unknown)"))
  }

  if (any( (grepl("Aad21122h", res_raw, fixed = TRUE) &
            !grepl("\\*NCC/3=O|\\*NCCO/3=O", res_raw)) )) {
    return(.skip_unsupported("Aad21122h", "Kdn/Neu"))
  }
  
  
  res_raw <- sapply(res_raw, function(x)
    WURCS_RES_MAP[[ names(WURCS_RES_MAP)[ str_detect(x, names(WURCS_RES_MAP)) ][1] ]]
  )
  
  ## RES-Sequence -------
  res_idx <- as.integer(strsplit(parts[3], "-", fixed = TRUE)[[1]])
  node <- res_raw[res_idx]
  
  ## LIN --------
  lin_raw <- paste(parts[4:length(parts)], collapse = "/")
  edge <- unlist(strsplit(lin_raw, "_", fixed = TRUE))
  if (length(edge) != 0) { 
    
    edge <- unlist(strsplit(edge, "_", fixed = TRUE))
    
    check_ambiguous_topo <- sapply(edge, function(x) {
      # Strip the left side (residue letter + position digits/? + dash) so
      # that `right` is the right-hand specification only. The original
      # regex `^[A-Za-z?]*-` did not allow digits in the position part, so
      # e.g. `f1-e3|e6` never matched and the un-stripped string was
      # incorrectly classified as having multiple parent residues. Allow
      # `[0-9?]+` in the position component.
      right <- sub("^[A-Za-z?][0-9?]+-", "", x)
      parts <- unlist(strsplit(right, "\\|"))
      let <- sub("^([A-Za-z]).*", "\\1", parts)
      length(unique(let)) > 1
    })
    
    if (any(check_ambiguous_topo)) {
      warning(sprintf(
        "Ambiguous LIN topology (edges %s) in WURCS: %s",
        paste(edge[check_ambiguous_topo], collapse = ", "),
        w
      ))
      # 本当は複数のノード候補を;で繋げて返してambigous = TRUEで返すのが理想
      # これは未対応単糖ではなく、LIN のトポロジーが真に一意に決まらないケース。
      # expand_wurcs_to_toposeq() が全候補を列挙して扱う。
      return(list(node = NULL, edge = NULL, ambiguous = TRUE,
                  reason = "ambiguous_lin_topology"))

    } else {
      
      edge <- unlist(strsplit(edge, "\\|"))
      edge <- grep("^[A-Za-z][0-9?]+-[A-Za-z][0-9?]+", edge, value = TRUE)
      edge <- sub("^([A-Za-z])[0-9?]+-([A-Za-z])[0-9?]+.*$", "\\1-\\2", edge, perl = TRUE)
      # edge <- sub("^([A-Za-z])[0-9?]+-([A-Za-z])[0-9?]+", "\\1-\\2", edge, perl = TRUE)
      
    }
    
  }
  
  list(node = node, edge = edge, ambiguous = FALSE)
}


# ------------ IUPACのから必要な情報を抽出 ----------- #
build_tree_iupac <- function(iupac) {
  
  iupac <- gsub("\\([^()]*-[^()]*\\)", ";", iupac)
  pat <- "[A-Z][A-Za-z0-9]*|\\[|\\]"
  
  tokens <- str_extract_all(iupac, pat)[[1]]
  tokens <- rev(tokens)
  tokens <- chartr("[]", "][", tokens)
  
  node_ids   <- character()
  node_label <- character()
  edge_chr   <- character()
  
  stack     <- character()       # 分岐用スタック（親ノードIDを積む）
  last_node <- NA_character_     # 直前ノードID
  node_idx  <- 0L
  
  for (tok in tokens) {
    
    if (tok == "[") {
      # push: 合流点（いまの last_node）を記憶するだけ
      stack <- c(stack, last_node)
      next
    }
    
    if (tok == "]") {
      # pop: 分岐を閉じたら last_node を合流点に巻き戻す（LIFO）
      if (length(stack) == 0L) stop("Unbalanced brackets: ']' without matching '['")
      last_node <- stack[length(stack)]
      stack <- stack[-length(stack)]
      next
    }
    
    # ---- 糖ノード ----
    node_idx <- node_idx + 1L
    nid <- sprintf("N%02d", node_idx)
    
    node_ids   <- c(node_ids, nid)
    node_label <- c(node_label, tok)
    
    # 親は常に直前（last_node）
    parent <- last_node
    
    if (!is.na(parent)) {
      edge_chr <- c(edge_chr, sprintf("%s-%s", parent, nid))
    }
    
    last_node <- nid
  }
  
  # node: N1, N2, ... を name に持つ named vector
  node_vec <- node_label
  names(node_vec) <- node_ids
  
  list(
    node      = node_vec,   # ex. N1 N2 N3 ... "GlcNAc" "GlcNAc" "Man" ...
    edge      = edge_chr,   # ex. "N1-N2" "N2-N3" ...
    ambiguous = FALSE
  )
}


modify_tree_wurcs <- function(tree) {
  
  if (isTRUE(tree$ambiguous)) {
    return(tree)
  }
  
  names(tree$node) <- sprintf("N%02d", seq_along(tree$node))
  
  if (length(tree$edge) == 0) {
    return(tree)
  }
  
  p <- str_split(tree$edge, "-", simplify = T)[,1]
  c <- str_split(tree$edge, "-", simplify = T)[,2]
  p <- names(tree$node)[match(p, c(letters, LETTERS))]
  c <- names(tree$node)[match(c, c(letters, LETTERS))]
  tree$edge <- paste(p, c, sep = "-")
  return(tree)
}



tree_to_igraph <- function(tree) {
  
  if (isTRUE(tree$ambiguous)) {
    return(NULL)
  }
  
  if (all(is.na(tree$node))) {
    return(igraph::make_empty_graph(n = 0, directed = TRUE))
  }
  
  
  if (length(tree$edge) == 0) {
    g <- make_empty_graph(n = length(tree$node), directed = TRUE)
    V(g)$name <- names(tree$node)
    V(g)$label <- tree$node
    V(g)$col <- res_cols[V(g)$label]
    return(g)
  } else {
    
    edge_df <- do.call(rbind, strsplit(tree$edge, "-", fixed = TRUE)) %>% as.data.frame()
    colnames(edge_df) <- c("parent", "child")
    if (any(edge_df$parent > edge_df$child)) {
      edge_df <- edge_df %>%
        rowwise() %>%
        mutate(
          # parent > child のときは入れ替える
          p = min(parent, child),
          c = max(parent, child)
        ) %>%
        ungroup() %>%
        select(parent = p, child = c) %>%
        arrange(child)
    }
    
    vertice_df <- data.frame(name = names(tree$node), label = tree$node)
    
    # 親が子にないかチェック -> 捨てる
    if (any (edge_df$parent != "N01" & !(edge_df$parent %in% edge_df$child) )) {
      warning("Invalid glycan topology: child nodes exist whose parent cannot be determined.")
      return(NULL)
    }
    
    g <- graph_from_data_frame(edge_df, vertices = vertice_df, directed = TRUE)
    V(g)$col <- res_cols[V(g)$label]
    
    return(g)
  }
}



canonicalize_igraph <- function(g) {
  
  if (is.null(g)) {
    return(NULL)
  } 
  
  # empty graph（ノードが0）はそのまま返す
  if (igraph::vcount(g) == 0L) {
    return(g)
  }
  
  ahu_label <- function(g, v) {
    children <- neighbors(g, v, "out")
    node_label <- V(g)$label[v]
    if (length(children) == 0) return(paste0(node_label, "()"))
    
    # 子部分木の canonical form を再帰的に取得
    child_labs <- vapply(children, function(ch) ahu_label(g, ch), character(1))
    
    # AHU: 子の canonical を辞書順にソート
    child_labs <- sort(child_labs)
    
    # 自分のラベル + "(" + 子ラベル列 + ")"
    paste0(node_label, "(", paste(child_labs, collapse = ""), ")")
  }
  
  compute_ahu_labels <- function(g) {
    vapply(V(g), function(v) ahu_label(g, v), character(1))
  }
  
  traverse_ahu <- function(g, root="N01") {
    
    labels <- compute_ahu_labels(g)   # ここで一括計算しておく
    
    root_v <- V(g)[name == root]
    
    dfs <- function(v) {
      children <- neighbors(g, v, "out")
      
      # 葉なら自分だけ
      if (length(children) == 0) return(V(g)$name[v])
      
      extract_tokens <- function(x) {
        stringr::str_extract_all(x, "[A-Za-z]+")[[1]]
      }

      lab_child  <- labels[children]        # AHU ラベル（括弧を含む完全形）
      len_child  <- vapply(lab_child, function(x) length(extract_tokens(x)),
                           integer(1))      # AHU ラベル長 = 部分木のノード数

      # 長さで降順（大きい部分木を優先）→ 同じ長さなら AHU ラベル全体で昇順。
      # 括弧を剥がした文字列で比較してはいけない。直鎖の腕と分岐の腕が同じ
      # トークン列に潰れて tie が破れず、order() の安定性で入力順が漏れる。
      # その結果、同一分子でも WURCS の残基ラベル割当順によって別の toposeq が
      # 出力される（かつ一方が別分子の toposeq と衝突しうる）。
      ord <- order(-len_child, lab_child)
      ch_sorted <- children[ord]
      
      c(
        V(g)$name[v],
        unlist(lapply(ch_sorted, dfs), use.names = FALSE)
      )
    }
    
    dfs(root_v)
  }
  
  vertex_ord <- traverse_ahu(g, "N01")
  
  g <- permute(g, match(V(g)$name, vertex_ord))
  V(g)$name <- sprintf("N%02d", seq_along(V(g)$name))
  g <- permute(g, sub("N", "", V(g)$name))
  
  return(g)
}


## ---------- ✦ igraph to pGlycoSeq ✦ ----------

igraph_to_toposeq <- function(g) {
  
  if (is.null(g)) {
    return(NULL)
  } 
  
  if (igraph::vcount(g) == 0L) {
    return(NA_character_)
  }
  
  seq <- character(0)
  root <- "N01"  
  depth <- distances(g, v = root, to = V(g), mode = "out")[1, ]

  for (k in seq_along(V(g)$name)) {
    cur_node  <- V(g)$name[k]
    cur_label <- V(g)$label[k]
    parent_node <- neighbors(g, cur_node, mode = "in") %>% as_ids()
    prev_node <- V(g)$name[k - 1]
    parent_depth <- depth[parent_node]
    cur_depth <- depth[prev_node]
    
    if (k == 1) {
      out <- paste0("(", cur_label)
      
    } else if (prev_node == parent_node) {
      out <- paste0("(", cur_label)
      
    } else {
      n_closeBracket <- cur_depth - parent_depth
      out <- paste0(strrep(")", n_closeBracket), "(", cur_label)
    }
    
    seq <- paste0(seq, out)
    
  }
  final_depth <- depth[cur_node]
  final_depth <- ifelse(is.na(final_depth), 0, final_depth)
  seq <- paste0(seq, strrep(")", final_depth + 1) )

  return(seq)
}


build_tree_toposeq <- function(toposeq) {
  
  if (is.na(toposeq)) {
    return(list(node = NA_character_, edge = NA_character_, ambiguous = FALSE))
  }

  tokens <- str_split(toposeq, "(?=[()])|(?<=[()])")[[1]]
  tokens <- tokens[tokens != ""]  
  
  node_ids   <- character()
  node_label <- character()
  edge_chr   <- character()
  
  stack     <- character()       # 分岐用スタック（親ノードIDを積む）
  last_node <- NA_character_     # 直前ノードID
  node_idx  <- 0L
  
  for(tok in tokens) {
    
    if(tok == "(") {
      stack <- c(stack, last_node)
      next
    }
    
    if (tok == ")") {
      last_node <- stack[length(stack)]
      stack <- stack[-length(stack)]
      next
    }
    
    node_idx <- node_idx + 1L
    nid <- sprintf("N%02d", node_idx)
    node_ids   <- c(node_ids, nid)
    node_label <- c(node_label, tok)
    parent <- last_node
    if (!is.na(parent)) {
      edge_chr <- c(edge_chr, sprintf("%s-%s", parent, nid))
    }
    last_node <- nid
  }
  
  # node: N1, N2, ... を name に持つ named vector
  node_vec <- node_label
  names(node_vec) <- node_ids

  list(
    node      = node_vec,   # ex. N1 N2 N3 ... "GlcNAc" "GlcNAc" "Man" ...
    edge      = edge_chr,   # ex. "N1-N2" "N2-N3" ...
    ambiguous = FALSE
  )
}


# ------------ GPSeeker g-Linkage parser (linkage-prefix + single-letter monosaccharides) ----------- #
# GPSeeker (Tian group) writes glycan structure as a parenthetical tree
# where each token is <linkage_prefix><letter>:
#   Y=GlcNAc(HexNAc), M=Man, L=Gal, G=Glc, F=Fuc, S=NeuAc, V=GalNAc
# The 2-digit linkage prefix encodes anomeric position + linkage position
# (e.g. "01Y" root, "41Y" beta-1,4-Y, "61M" alpha-1,6-M, "32S" alpha-2,3-S).
# For glyEnzMap topology-only matching against existing GlyEnzDB/GlyGen
# triplets, the linkage prefix is stripped — only the single-letter
# token remains as the node label, mirroring how Kizuka toposeq is
# linkage-free.
build_tree_gpseeker <- function(gpseeker_str) {

  if (length(gpseeker_str) != 1L || is.na(gpseeker_str) || !nzchar(gpseeker_str)) {
    return(list(node = character(0), edge = character(0), ambiguous = FALSE))
  }

  # Tokenize: parenthetical splits keep '(' and ')' as their own tokens;
  # everything between is a run that contains one or more
  # <2-digit><letter> chunks back-to-back (e.g. "01Y41Y41M" = 3 tokens).
  raw <- str_split(gpseeker_str, "(?=[()])|(?<=[()])")[[1]]
  raw <- raw[raw != ""]

  # Expand non-bracket runs into individual <2-digit><letter> tokens.
  tokens <- character()
  for (chunk in raw) {
    if (chunk == "(" || chunk == ")") {
      tokens <- c(tokens, chunk)
    } else {
      # Each glyco token is two digits + one letter.
      pieces <- str_extract_all(chunk, "[0-9]{2}[A-Za-z]")[[1]]
      if (length(pieces) == 0L ||
          paste(pieces, collapse = "") != chunk) {
        warning("Cannot fully tokenize chunk '", chunk,
                "' in g-Linkage: ", gpseeker_str)
        return(list(node = character(0), edge = character(0), ambiguous = TRUE))
      }
      tokens <- c(tokens, pieces)
    }
  }

  # Strip linkage prefix to keep just the letter.
  is_glyco <- !(tokens %in% c("(", ")"))
  letters_v <- substr(tokens[is_glyco], 3L, 3L)

  # Validate letter alphabet.
  valid_letters <- c("Y", "M", "L", "G", "F", "S", "V")
  bad_let <- setdiff(letters_v, valid_letters)
  if (length(bad_let) > 0L) {
    warning("Unknown GPSeeker monosaccharide letter(s): ",
            paste(unique(bad_let), collapse = ", "),
            " in: ", gpseeker_str)
    return(list(node = character(0), edge = character(0), ambiguous = TRUE))
  }

  # Replace the glyco tokens with their stripped letter for the parser loop.
  tokens[is_glyco] <- letters_v

  node_ids   <- character()
  node_label <- character()
  edge_chr   <- character()

  stack     <- character()
  last_node <- NA_character_
  node_idx  <- 0L

  for (tok in tokens) {
    if (tok == "(") {
      stack <- c(stack, last_node)
      next
    }
    if (tok == ")") {
      if (length(stack) == 0L) {
        warning("Unbalanced ')' in g-Linkage: ", gpseeker_str)
        return(list(node = character(0), edge = character(0), ambiguous = TRUE))
      }
      last_node <- stack[length(stack)]
      stack <- stack[-length(stack)]
      next
    }
    node_idx <- node_idx + 1L
    nid <- sprintf("N%02d", node_idx)
    node_ids   <- c(node_ids, nid)
    node_label <- c(node_label, tok)
    parent <- last_node
    if (!is.na(parent)) {
      edge_chr <- c(edge_chr, sprintf("%s-%s", parent, nid))
    }
    last_node <- nid
  }

  node_vec <- node_label
  names(node_vec) <- node_ids

  list(
    node      = node_vec,
    edge      = edge_chr,
    ambiguous = FALSE
  )
}


# ------------ GPSeeker single-letter -> full-name label expansion ----------- #
# Maps GPSeeker single-letter labels to glyEnzMap-canonical labels.
# No positional logic needed: GPSeeker uses distinct letters for Glc (G)
# and Gal (L), and for GlcNAc (Y) and GalNAc (V), unlike pGlyco's H
# which is ambiguous. Biology is encoded in the letter choice itself.
assign_gpseeker_labels <- function(g) {
  if (is.null(g) || igraph::vcount(g) == 0L) return(g)

  letter_map <- c(Y = "GlcNAc", M = "Man",   L = "Gal",
                  G = "Glc",    F = "Fuc",   S = "NeuAc",
                  V = "GalNAc")

  labs    <- V(g)$label
  new_lab <- labs
  hit     <- labs %in% names(letter_map)
  new_lab[hit] <- letter_map[labs[hit]]
  if (!all(hit)) {
    unknowns <- unique(labs[!hit])
    warning("assign_gpseeker_labels: unmapped node label(s): ",
            paste(unknowns, collapse = ", "))
  }

  V(g)$label <- new_lab
  V(g)$col   <- res_cols[V(g)$label]
  g
}


substructure_match_one <- function(g_small, g_large) {
  
  label_union <- union(igraph::V(g_small)$label, igraph::V(g_large)$label)
  color_tbl <- setNames(seq_along(label_union), label_union)
  
  igraph::V(g_small)$color <- unname(color_tbl[igraph::V(g_small)$label])
  igraph::V(g_large)$color <- unname(color_tbl[igraph::V(g_large)$label])
  
  igraph::subgraph_isomorphic(g_small, g_large, method = "vf2")
}

substructure_match_batch <- function(g_small_list, g_large_list, mc.cores = 4) {
  
  res_list <- mclapply(
    g_small_list,
    function(gs) {
      vapply(
        g_large_list,
        function(gl) substructure_match_one(gs, gl),
        logical(1)
      )
    }, mc.cores = mc.cores
  )
  
  do.call(rbind, res_list)
}


## ---------- ✦ Enzyme Group Labeling ✦ ----------

#' Label an enzyme pattern group with a human-readable name
#'
#' Given a character vector of enzyme gene names that share the same
#' glycan-interaction profile, produce a concise label.
#'
#' Rules (applied in order):
#'   1. Single-member group → return the member name directly.
#'   2. Manual mapping for known multi-family groups (MGAT4/5, etc.).
#'   3. If all members share a single prefix (letters before trailing
#'      digits + optional letter), collapse the numeric suffixes with "/".
#'   4. Fallback: paste all member names with "/".
#'
#' @param members Character vector of gene names (e.g., c("FUT3", "FUT5", "FUT6"))
#' @return A single character string label
extract_enzyme_group_label <- function(members) {
  if (length(members) == 1) return(members)

  # Manual mapping for known complex groups

  manual_map <- list(
    "MGAT4/5"              = c("MGAT4A", "MGAT4B", "MGAT4C", "MGAT4D", "MGAT5"),
    "ST3GAL/ST6GAL/CMAS"   = c("CMAS", "ST3GAL1", "ST3GAL3", "ST3GAL4",
                                "ST3GAL6", "ST6GAL1", "ST6GAL2"),
    "B4GALT/NEU"           = c("B3GALT2", "B4GALT1", "B4GALT2", "B4GALT3",
                                "B4GALT4", "B4GALT5", "NEU1", "NEU2", "NEU3", "NEU4")
  )

  members_sorted <- sort(members)
  for (label in names(manual_map)) {
    if (identical(members_sorted, sort(manual_map[[label]]))) return(label)
  }

  # Auto-extract: strip trailing digits (+ optional letter) to get family prefix
  prefixes <- unique(gsub("[0-9]+[A-Z]?$", "", members))

  if (length(prefixes) == 1) {
    # All share the same prefix — collapse numeric suffixes with "/"
    return(paste0(prefixes,
                  paste(gsub(prefixes, "", members), collapse = "/")))
  }

  # Multiple prefixes with no manual mapping — list all members
  return(paste(members, collapse = "/"))
}


## ---------- ✦ Plot ✦ ----------

plot_glycan_tree <- function(g) {
  plot(
    g,
    layout = igraph::layout_as_tree(g, root = "N01"),
    vertex.label = igraph::V(g)$label,
    vertex.color = igraph::V(g)$col,
    vertex.size = 30,
    edge.arrow.size = 1,
    main = "Glycan Topology"
  )
}

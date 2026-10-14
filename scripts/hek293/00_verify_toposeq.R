# Verification: can every Huang 2021 Structure ID be converted to toposeq?
#
# Pipeline: GlyTouCan ID -> WURCS via gtcid2seqs -> toposeq via utils_formatting.R
# Output: pass/fail report per ID, cached WURCS table.

# Run from the repository root.

suppressMessages({
  library(readxl); library(dplyr); library(httr); library(jsonlite); library(parallel); library(tibble)
})
source("R/utils_formatting.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

XLSX <- "data/Huang_2021/NIHMS1683936-supplement-5.xlsx"
CACHE_DIR <- "cache/hek293"
OUT_DIR   <- "output/hek293"
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)

# ---------- 1. read Table_S4_1 N-glycan rows only ----------
# Sheet layout: rows 1-2 are multi-level header, data rows 3..N.
# Table_S4_1 contains 73 rows total: rows 1..61 (data) are N-glycans,
# rows 62..73 are GSL (Hex2Cer-anchored). We keep only N-glycans.
# Table_S4_2 has a different schema (MS2 signature ions) — skip.
raw <- suppressMessages(read_excel(XLSX, sheet = "Table_S4_1",
                                   skip = 2, col_names = FALSE))
struct_name <- as.character(raw[[2]])   # "Structural name"
struct_id   <- as.character(raw[[5]])   # "Structure ID"

is_gsl <- grepl("Cer", struct_name, fixed = TRUE)
n_total <- nrow(raw)
n_ngly  <- sum(!is_gsl)

cat(sprintf("Total Table_S4_1 rows: %d (N-glycan: %d, GSL: %d)\n",
            n_total, n_ngly, sum(is_gsl)))

ng_ids <- struct_id[!is_gsl]
all_ids <- unlist(strsplit(ng_ids, ";"))
all_ids <- trimws(all_ids)
all_ids <- all_ids[!is.na(all_ids) & nzchar(all_ids) & all_ids != "NA"]

# Sanity: every token should match GlyTouCan accession format G[A-Z0-9]{7}
bad_format <- !grepl("^G[A-Z0-9]{7}$", all_ids)
if (any(bad_format)) {
  cat(sprintf("WARNING: %d tokens not in GlyTouCan format:\n",
              sum(bad_format)))
  print(unique(all_ids[bad_format]))
}
all_ids <- all_ids[!bad_format]
unique_ids <- unique(all_ids)

cat(sprintf("Total N-glycan ID tokens: %d\n", length(all_ids)))
cat(sprintf("Unique GlyTouCan IDs: %d\n", length(unique_ids)))

# ---------- 2. fetch WURCS via gtcid2seqs ----------
get_wurcs <- function(gtcid) {
  url <- "https://api.glycosmos.org/sparqlist/gtcid2seqs"
  res <- tryCatch(httr::GET(url, query = list(gtcid = gtcid), httr::timeout(20)),
                  error = function(e) NULL)
  if (is.null(res) || httr::status_code(res) != 200) return(NA_character_)
  txt <- httr::content(res, "text", encoding = "UTF-8")
  dat <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
  if (is.null(dat) || is.null(dat$wurcs)) return(NA_character_)
  w <- dat$wurcs[!is.na(dat$wurcs)]
  if (length(w) == 0) return(NA_character_)
  w[1]
}

cache_file <- file.path(CACHE_DIR, "gtc_to_wurcs.rds")
if (file.exists(cache_file)) {
  cache <- readRDS(cache_file)
} else {
  cache <- list()
}

new_ids <- setdiff(unique_ids, names(cache))
cat(sprintf("Need to fetch %d new IDs (cached: %d)\n",
            length(new_ids), length(cache)))

if (length(new_ids) > 0) {
  # serial: httr + macOS objc fork interaction breaks mclapply
  fetched <- vector("list", length(new_ids))
  names(fetched) <- new_ids
  for (i in seq_along(new_ids)) {
    fetched[[i]] <- get_wurcs(new_ids[i])
    Sys.sleep(0.05)
    if (i %% 20 == 0) cat(sprintf("  fetched %d / %d\n", i, length(new_ids)))
  }
  cache <- c(cache, fetched)
  saveRDS(cache, cache_file)
}

wurcs_vec <- unlist(cache[unique_ids])
n_wurcs_ok <- sum(!is.na(wurcs_vec))
n_wurcs_fail <- sum(is.na(wurcs_vec))
cat(sprintf("WURCS lookup: %d ok / %d fail\n", n_wurcs_ok, n_wurcs_fail))

failed_ids_wurcs <- unique_ids[is.na(wurcs_vec)]
if (length(failed_ids_wurcs) > 0) {
  cat("WURCS-failed IDs:\n")
  print(failed_ids_wurcs)
}

# ---------- 3. WURCS -> toposeq via existing pipeline ----------
ok_ids   <- unique_ids[!is.na(wurcs_vec)]
ok_wurcs <- wurcs_vec[!is.na(wurcs_vec)]

convert_one <- function(w) {
  stage <- "start"
  out <- suppressWarnings(suppressMessages(tryCatch({
    stage <- "build_tree_wurcs"; tr <- build_tree_wurcs(w)
    if (is.list(tr) && isTRUE(tr$ambiguous)) {
      return(list(toposeq = NA_character_, status = "ambiguous_skip",
                  error = NA_character_, stage_fail = NA_character_))
    }
    stage <- "modify_tree_wurcs"; tr <- modify_tree_wurcs(tr)
    stage <- "tree_to_igraph"; g  <- tree_to_igraph(tr)
    stage <- "canonicalize_igraph"; g  <- canonicalize_igraph(g)
    stage <- "igraph_to_toposeq"; ts <- igraph_to_toposeq(g)
    list(toposeq = ts, status = "ok",
         error = NA_character_, stage_fail = NA_character_)
  }, error = function(e) {
    list(toposeq = NA_character_, status = "error",
         error = conditionMessage(e),
         stage_fail = stage)
  })))
  out
}

# serial too (canonicalize uses igraph internals that may not be fork-safe)
results <- lapply(ok_wurcs, convert_one)

tab <- tibble(
  glytoucan  = ok_ids,
  wurcs      = ok_wurcs,
  toposeq    = vapply(results, function(r) r$toposeq %||% NA_character_, character(1)),
  status     = vapply(results, function(r) r$status %||% NA_character_, character(1)),
  error      = vapply(results, function(r) r$error   %||% NA_character_, character(1)),
  stage_fail = vapply(results, function(r) r$stage_fail %||% NA_character_, character(1))
)

cat("\n===== toposeq conversion result =====\n")
print(table(tab$status, useNA = "ifany"))
cat(sprintf("End-to-end ok: %d / %d (%.1f%%)\n",
            sum(tab$status == "ok"), length(unique_ids),
            100 * sum(tab$status == "ok") / length(unique_ids)))

if (any(tab$status == "ambiguous_skip")) {
  cat(sprintf("\n%d structures classified as ambiguous_skip (truly multi-parent LIN, biologically uncertain assignment).\n",
              sum(tab$status == "ambiguous_skip")))
  cat("First 5:\n")
  print(tab |> dplyr::filter(status == "ambiguous_skip") |>
          dplyr::select(glytoucan, wurcs) |> head(5))
}

if (any(tab$status == "error")) {
  cat("\nHard errors (unexpected):\n")
  print(tab |> dplyr::filter(status == "error") |>
          dplyr::count(stage_fail, error, sort = TRUE))
}

# ---------- 4. save report ----------
saveRDS(tab,
        file.path(OUT_DIR, "verify_toposeq_table.rds"))
write.csv(tab,
          file.path(OUT_DIR, "verify_toposeq_table.csv"),
          row.names = FALSE)


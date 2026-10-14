# Regenerate Supplementary Table S1 (enzyme pattern groups) from the CURRENT
# 16-group analysis. The previous output/TableS1_enzyme_pattern_groups.csv used a
# stale grouping that merged FUT3/4/5/6/7/9 into one group, contradicting the
# current GMT/resGSEA (FUT3/5/6 and FUT4/7/9 are separate) and the manuscript.
# Source of truth: output/enzyme_pattern_groups_list.rds (saved by step3_build_pk.R).
# Output: output/TableS1_enzyme_pattern_groups.csv
suppressMessages(library(data.table))

epg     <- readRDS("output/enzyme_pattern_groups_list.rds")     # profile -> member enzymes (16 groups)
triplet <- as.data.frame(readRDS("output/triplet.rds"))
rcg     <- readRDS("output/res_clean_group.rds")                 # enzyme x glycan-group binary matrix

## ---- group label (verbatim from src/scripts/step3_build_pk.R) --------------
extract_enzyme_group_label <- function(members) {
  if (length(members) == 1) return(members)
  manual_map <- list(
    "MGAT4/5" = c("MGAT4A", "MGAT4B", "MGAT4C", "MGAT4D", "MGAT5"),
    "ST3GAL/ST6GAL/CMAS" = c("CMAS", "ST3GAL1", "ST3GAL3", "ST3GAL4",
                              "ST3GAL6", "ST6GAL1", "ST6GAL2"),
    "B4GALT/NEU" = c("B3GALT2", "B4GALT1", "B4GALT2", "B4GALT3",
                      "B4GALT4", "B4GALT5", "NEU1", "NEU2", "NEU3", "NEU4")
  )
  ms <- sort(members)
  for (label in names(manual_map)) if (identical(ms, sort(manual_map[[label]]))) return(label)
  prefixes <- unique(gsub("[0-9]+[A-Z]?$", "", members))
  if (length(prefixes) == 1) return(paste0(prefixes, paste(gsub(prefixes, "", members), collapse = "/")))
  paste(members, collapse = "/")
}
classify_enzyme <- function(ec) {
  if (is.na(ec) || ec == "") return("Other")
  if (grepl("transferase", ec, ignore.case = TRUE)) return("GT")
  if (grepl("Glycosidase|hydrolase", ec, ignore.case = TRUE)) return("GH")
  "Other"
}
cls_map <- setNames(sapply(triplet$Enzyme_class, classify_enzyme), triplet$Enzyme_gene)

rows <- rbindlist(lapply(epg, function(members) {
  label <- extract_enzyme_group_label(members)
  rep_e <- members[members %in% rownames(rcg)][1]
  nglyc <- if (!is.na(rep_e)) sum(rcg[rep_e, ] > 0) else NA_integer_
  classes <- paste(sort(unique(unname(cls_map[members]))), collapse = "/")
  data.table(Pattern_group = label,
             N_enzymes     = length(members),
             N_glycan_groups = nglyc,
             Enzyme_class  = classes,
             Member_enzymes = paste(sort(members), collapse = ", "))
}))

## order: multi-enzyme groups first (by size desc), then singletons; alpha within
setorder(rows, -N_enzymes, Pattern_group)
fwrite(rows, "output/TableS1_enzyme_pattern_groups.csv")
cat("WROTE output/TableS1_enzyme_pattern_groups.csv  (", nrow(rows), "groups )\n")
print(rows)

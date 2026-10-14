library(readr)
library(httr)
library(jsonlite)
library(dplyr)
library(parallel)
source("R/utils_formatting.R")

eg <- read_csv("data/GlyGen/glycan_pathway_glycotree.csv")
eg <- eg %>% filter(enzyme_tax_name == "Homo sapiens")


enz <- read_csv("data/GlyGen/glycan_enzyme.csv.gz")
eg$gene_name <- enz$gene_name[match(eg$enzyme_uniprotkb_ac, enz$uniprotkb_ac)]
eg$enzyme_class <- enz$recommended_name_full[match(eg$enzyme_uniprotkb_ac, enz$uniprotkb_ac)]

# Enzymeのフィルタリング -------------------------------------------------
enz <- as.data.frame(enz)
enz_filt <- enz %>% filter(species == "Homo sapiens")  
disc <- unique(enz_filt$recommended_name_full)
kw <- c( "Polypeptide N-acetylgalactosaminyltransferase",
         "N-acetylgalactosaminyltransferase 7",
         "N-acetylgalactosaminide",
         "N-acetylgalactosamine",
         "N-acetyl-galactosaminide",
         "fringe",
         "C1GALT1-specific chaperone 1", 
         "O-glycosyl", "O-glucosyl", "O-mannosyl", "O-linked", "O-fucosyl", 
         "xylosyl", "Lactosylceramide" )
hit <- Reduce(`|`, lapply(kw, function(k) grepl(k, disc, ignore.case = TRUE)))
out_disc <- disc[hit]
enz_filt <- enz_filt %>% filter(!(recommended_name_full %in% out_disc))
use_enzyme <- unique(enz_filt$uniprotkb_ac)
# ------------------------------------------------------------------------

eg <- eg %>% filter(eg$enzyme_uniprotkb_ac %in% use_enzyme)

eg$GlyGendbID <- sprintf("GLYGEN%05d", seq_len(nrow(eg)))

# STEP1: glytoucanidをwurksに変換する
glytoucan <- unique(c(eg$glytoucan_ac, eg$source, eg$target)) 

get_wurcs_from_glytoucan <- function(glytoucan_id) {
  url <- "https://api.glycosmos.org/sparqlist/gtcid2seqs"
  
  res <- httr::GET(url, query = list(gtcid = glytoucan_id))
  txt <- httr::content(res, "text", encoding = "UTF-8")
  dat <- jsonlite::fromJSON(txt)
  wurcs <- dat$wurcs[!is.na(dat$wurcs)][1]

  return(wurcs)
}

# Load from cache if exists, otherwise fetch from API
if (file.exists("cache/glygen_enzglycan_wurcs.rds")) {
  wurcs <- readRDS("cache/glygen_enzglycan_wurcs.rds")
} else {
  message("Fetching WURCS from GlyCosmos API (this may take a while)...")
  wurcs <- sapply(glytoucan, get_wurcs_from_glytoucan)
  saveRDS(wurcs, "cache/glygen_enzglycan_wurcs.rds")
}

# STEP2: igraphに変換する
toposeq <- mclapply(wurcs, function(x) {
  x %>%
    build_tree_wurcs() %>%
    modify_tree_wurcs() %>%
    tree_to_igraph() %>%
    canonicalize_igraph() %>%
    igraph_to_toposeq()
}, mc.cores = 20)
table(sapply(toposeq, is.null))
toposeq_filt <- toposeq[!sapply(toposeq, is.null)]

keep_idx <- eg$source %in% names(toposeq_filt) & eg$target %in% names(toposeq_filt)
eg <- eg[keep_idx,]

toposeq_filt <- unlist(toposeq_filt)

eg <- eg %>%
  mutate(source_toposeq = toposeq_filt[source],
         target_toposeq = toposeq_filt[target])


triplet_glygen <- data.frame(Enzyme = eg$enzyme_uniprotkb_ac, Reactant = eg$source_toposeq, Product = eg$target_toposeq, Enzyme_gene = eg$gene_name, Enzyme_class = eg$enzyme_class, GlyGendbID = eg$GlyGendbID)

# Save intermediate data for Figure2.R
saveRDS(triplet_glygen, "output/triplet_glygen.rds")


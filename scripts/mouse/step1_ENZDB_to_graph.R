library(data.table)
library(stringr)
library(tidyverse)
library(igraph)
source("R/utils_formatting.R")

# -------------- #
#   glycoEnzDB
# -------------- #
enzdb <- fread("data/GlyEnzDB/gene_other_information.csv", data.table = F)




# ------------ ENZDBの整理 ----------- #
table(enzdb$GlyClass)
enzdb_N <- enzdb[grep("N_linked", enzdb$GlyClass), ]

product1 <- enzdb_N$Product
product1 <- str_split(product1, ";| or ")
names(product1) <- enzdb_N$GeneSymbol
product1 <- enframe(product1, name = "enzyme", value = "product") %>% unnest(product)

reactant1 <- enzdb_N$Reactant
reactant1 <- sub(";$", "", reactant1)
reactant1 <- str_split(reactant1, ";| or ")
names(reactant1) <- enzdb_N$GeneSymbol
reactant1 <- enframe(reactant1, name = "enzyme", value = "reactant") %>% unnest(reactant)
pair1 <- cbind(product1, dplyr::select(reactant1, -enzyme))


product2 <- enzdb_N$Product_1
product2 <- str_split(product2, " or ")
names(product2) <- enzdb_N$GeneSymbol
product2 <- enframe(product2, name = "enzyme", value = "product") %>% unnest(product)

reactant2 <- enzdb_N$Substrate
reactant2 <- str_split(reactant2, " or ")
names(reactant2) <- enzdb_N$GeneSymbol
reactant2 <- enframe(reactant2, name = "enzyme", value = "reactant") %>% unnest(reactant)
add_row   <- reactant2 %>% filter(enzyme == "ABO")   # ① 行を取り出す
reactant2 <- reactant2 %>%
  add_row(.after   = 1, !!! reactant2[1, ]) 
pair2 <- cbind(product2, dplyr::select(reactant2, -enzyme))

enzdb_tb <- rbind(pair1, pair2)
enzdb_tb$enzclass <- enzdb_N$EnzClass[match(enzdb_tb$enzyme, enzdb_N$GeneSymbol)]




# ------------ ENZDBのフィルター ----------- #
# gdbにSulfateはないので除く [OK]
enzdb_filt <- enzdb_tb %>% filter(enzclass != "Sulfate synthase")
# ""を除く [OK]
noreact <- c("none","None","no reaction","No reaction","no clearly defined reaction available","no reaction (accessory protein)")
enzdb_filt <- enzdb_filt %>% filter(! (product == "" & reactant %in% noreact))




# ------------ 使用しない反応をフィルター ----------- #
exclude_patterns <- c(
  "-R",
  "Sia",
  "P-Man",
  "Man6P",
  "\\.\\.\\.",   # "..."
  "…",           # 全角…
  "GlcNAc6S",
  "Gal6S",
  "6SGlcNAc",
  "GlcA3S",
  "\\[6S\\]",
  "\\[3S\\]",
  "Ser",
  "Cer",
  "GlcA", # GlcN, Fuc1Pは保持
  ".+GlcNAc6P.+", 
  ".+GlcNAc\\(6P\\).+"
)
pattern_regex <- paste(exclude_patterns, collapse = "|")
enzdb_filt <- enzdb_filt %>% 
  filter(!str_detect(product, pattern_regex)) %>%
  filter(!str_detect(reactant, pattern_regex))



# ------------ IUPACの修正 ----------- #
repl_rules <- c(
  "α"                             = "a",          # α → a
  "β"                             = "b",          # β → b
  "Fuc\\(1P\\)"                    = "Fuc1P",      # Fuc(1P) → Fuc1P
  "GlcNAc\\(6P\\)"                 = "GlcNAc6P",   # GlcNAc(6P) → GlcNAc6P
  "GlcN\\(6P\\)"                   = "GlcN6P",     # GlcN(6P) → GlcN6P
  "Dol-PP|Dol-P-P|-PP-Dol|PP-Dol"  = "",           # Dol-PP 系削除
  "Asn|-Asn|\\)Asn"                = "",           # Asn 削除
  "Dol-P|-P-Dol"                   = "",           # Dol-P 系削除
  "GDP-|CMP-"                      = "",           # GDP- / CMP- prefix 削除
  "NeuAc"                          = "Neu5Ac"      # NeuAc → Neu5Ac
)
repl_rules_6P <- c(
  "6P" = ""   # GlcNAc6P / GlcN6P → GlcNAc / GlcN など
)

enzdb_filt <- enzdb_filt %>%
  mutate(
    product  = str_replace_all(product,  repl_rules),
    reactant = str_replace_all(reactant, repl_rules)
  ) %>%
  mutate(
    product  = str_replace_all(product,  repl_rules_6P),
    reactant = str_replace_all(reactant, repl_rules_6P)
  ) %>%
  mutate(
    product  = trimws(product),
    reactant = trimws(reactant)
  )


# ============================== #
#           　完成
enzdb_filt
enzdb_filt$ENZdbID <- sprintf("GENZ%05d", seq_len(nrow(enzdb_filt)))
# ============================== #



enzdb_reactant_toposeq <- sapply(enzdb_filt$reactant, function(x) {
  x %>%
    build_tree_iupac() %>%
    tree_to_igraph() %>%
    canonicalize_igraph() %>%
    igraph_to_toposeq()
})

enzdb_product_toposeq <- sapply(enzdb_filt$product, function(x) {
  x %>%
    build_tree_iupac() %>%
    tree_to_igraph() %>%
    canonicalize_igraph() %>%
    igraph_to_toposeq()
})


triplet_enzdb <- data.frame(Enzyme_gene = enzdb_filt$enzyme, Reactant = enzdb_reactant_toposeq, Product = enzdb_product_toposeq, Enzyme_class = enzdb_filt$enzclass, ENZdbID = enzdb_filt$ENZdbID)

# Save intermediate data for Figure2.R
saveRDS(triplet_enzdb, "output/triplet_enzdb.rds")

# GlcNAc(6P)、GlcN(6P)は、GlcNAc6P、GlcN6Pに変換し、
# 単体ではなく糖鎖の中に組み込まれたGlcNAc6Pはgdbにないので除く [要確認]
# 単体ではなく糖鎖の中に組み込まれたGlcNAc6PやGlcNAc(6P)はgdbにないので除く [要確認]
# Fuc(1P)はFuc1Pに変換する
# PP-DolとAsnはgdbにはない。切り取ってもマッチする時は根本でマッチするはずだからトリミング
# 供与体も一応保持。-P-Dolをトリミングしても良いと思う。
# transferaseでproductが単糖１つという条件
# nucleotide sugar("GDP-Fuc", "CMP-Neu5Ac")の処理
# transferaseとして扱いproductが単糖１つという条件
# GlcNAc6Pの処理
# gdbにはGlcNAc6Pはない。UDP-GlcNAcすなわちGlcNAcの供給を促進すると仮定する。
# transferaseとして扱いproductが単糖１つという条件
# NeuAcをNeu5Acに置換


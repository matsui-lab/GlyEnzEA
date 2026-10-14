#!/usr/bin/env Rscript
# Emit the data behind the Huang-KO Supplementary Table(s):
#   (1) full GSEA results for all 7 KOs (long), with the a-priori target flagged
#   (2) WT N-glycome composition by antennary class (exclusion rationale)
#   (3) high-mannose (M5..M9 + Glc1M9) total rel % per condition (substrate side)
# Read-only inputs; writes TSVs under .../output/ for the docx builder.

suppressMessages({library(data.table)})
od <- "output/hek293"
r  <- readRDS(file.path(od, "resGSEA_by_ko.rds"))
rm <- as.data.frame(readRDS(file.path(od, "rowmeta.rds")))

ko_order <- c("MGAT1_KO","MGAT2_KO","T_KO_MANA1_A2_B1_","B4GALNT3_4_KO",
              "SLC35C1_KO","MGAT4A_4B_KO","MGAT5_KO")
ko_lab <- c(MGAT1_KO="MGAT1-KO", MGAT2_KO="MGAT2-KO",
            T_KO_MANA1_A2_B1_="T-KO (MAN1A1&A2&B1)", B4GALNT3_4_KO="B4GALNT3&4-KO",
            SLC35C1_KO="SLC35C1-KO", MGAT4A_4B_KO="MGAT4A&4B-KO", MGAT5_KO="MGAT5-KO")
target <- c(MGAT1_KO="MGAT1/2/5B", MGAT2_KO="MGAT1/2/5B",
            T_KO_MANA1_A2_B1_="MAN1A1/MAN1A2/MAN1C1",
            B4GALNT3_4_KO="B3GALNT2/B4GALNT3/B4GALNT4", SLC35C1_KO="FUT8",
            MGAT4A_4B_KO="MGAT4/5", MGAT5_KO="MGAT4/5")
exp_dir <- c(MGAT1_KO="DOWN", MGAT2_KO="DOWN", T_KO_MANA1_A2_B1_="UP",
             B4GALNT3_4_KO="DOWN", SLC35C1_KO="DOWN(FUT)",
             MGAT4A_4B_KO="DOWN", MGAT5_KO="DOWN")

## (1) full GSEA long table
g <- rbindlist(lapply(ko_order, function(k){
  x <- as.data.frame(r[[k]])
  data.table(KO=ko_lab[k], pathway=x$pathway, NES=round(x$NES,2),
             p=signif(x$pval,3), q=signif(x$qvalue,3), size=x$size,
             is_target=(x$pathway==target[k]))
}))
g <- g[order(factor(KO, levels=ko_lab[ko_order]), q)]
fwrite(g, file.path(od,"suppT_gsea_all.tsv"), sep="\t")

## (2) WT composition by antennary class
nm  <- rm[!is.na(rm$name),]
ant <- sapply(nm$name, function(x){ m<-regmatches(x,regexpr("N\\d",x))
  if(length(m)==0) "other" else c(N1="N1 (mono)",N2="N2 (bi)",N3="N3 (tri)",
    N4="N4 (tetra)")[m] })
comp <- aggregate(nm$avg_WT, list(class=ant), sum)
comp <- comp[order(-comp$x),]; comp$x <- round(comp$x,1)
names(comp)[2] <- "WT_pct"
fwrite(comp, file.path(od,"suppT_composition.tsv"), sep="\t")

## (3) high-mannose total rel % per condition
hm_names <- c("M5N2","M6N2","M7N2","M8N2","M9N2","Glc1M9N2")
avgcols  <- grep("^avg_", colnames(rm), value=TRUE)
hm <- data.table(condition=sub("^avg_","",avgcols),
                 highmannose_pct=round(colSums(rm[rm$name %in% hm_names, avgcols],
                                               na.rm=TRUE),1))
fwrite(hm, file.path(od,"suppT_highmannose.tsv"), sep="\t")

## small per-KO target verdict summary (for the caption / quick reference)
verdict <- rbindlist(lapply(ko_order, function(k){
  x <- as.data.frame(r[[k]]); row <- x[x$pathway==target[k],]
  if(nrow(row)==0) return(NULL)
  q <- row$qvalue[1]; p <- row$pval[1]
  v <- if(q<0.05) "hit" else if(q<0.10 || p<0.05) "partial" else "miss"
  data.table(KO=ko_lab[k], target=target[k], expected=exp_dir[k],
             NES=round(row$NES[1],2), q=signif(q,3), verdict=v)
}))
fwrite(verdict, file.path(od,"suppT_target_verdict.tsv"), sep="\t")

cat("wrote suppT_gsea_all / composition / highmannose / target_verdict TSVs\n")
print(verdict); print(comp); print(hm)

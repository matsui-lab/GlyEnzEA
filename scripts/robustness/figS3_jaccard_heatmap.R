## Supplementary Figure S3 — Jaccard similarity heatmap of the 16 enzyme pattern groups.
## Source: output/glycan2enzMap_enzyme_grouped.rds (GMT: group -> glycan member set).
suppressMessages({library(ggplot2);library(data.table)})
gmt <- readRDS("output/glycan2enzMap_enzyme_grouped.rds")
G <- names(gmt)
J <- outer(seq_along(G), seq_along(G), Vectorize(function(i,j){
  a<-gmt[[i]]; b<-gmt[[j]]; length(intersect(a,b))/length(union(a,b)) }))
dimnames(J) <- list(G,G)

## order rows/cols by hierarchical clustering on (1 - Jaccard)
ord <- hclust(as.dist(1-J), method="average")$order
Go <- G[ord]
df <- as.data.table(as.table(J)); setnames(df, c("g1","g2","jac"))
df[, g1 := factor(g1, levels=Go)]; df[, g2 := factor(g2, levels=Go)]

p <- ggplot(df, aes(g1, g2, fill=jac)) +
  geom_tile(color="white", linewidth=0.3) +
  geom_text(aes(label=ifelse(jac>=0.005, sprintf("%.2f", jac), "")), size=2.5,
            color=ifelse(df$jac>0.6,"white","grey20")) +
  scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B",
                       midpoint=0.5, limits=c(0,1), name="Jaccard") +
  coord_equal() +
  labs(x=NULL, y=NULL,
       title="Pairwise Jaccard similarity of enzyme pattern groups") +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        axis.text.y=element_text(size=7), panel.grid=element_blank())
ggsave("output/FigureS3_jaccard_heatmap.pdf", p, width=8.2, height=7)
cat("WROTE output/FigureS3_jaccard_heatmap.pdf\n")

## report ranges relevant to the main text claims
late <- c("ST3GAL/ST6GAL/CMAS","FUT3/5/6","FUT4/7/9","B4GALT/NEU","B3GALT5")
offdiag <- J[upper.tri(J)]
cat(sprintf("\nglobal off-diagonal Jaccard range: %.2f - %.2f\n", min(offdiag), max(offdiag)))
sub <- J[late,late]; cat("\nlate-pathway terminal groups submatrix:\n"); print(round(sub,3))
ld <- sub[upper.tri(sub)]
cat(sprintf("\nlate-pathway off-diagonal range: %.2f - %.2f\n", min(ld), max(ld)))
fwrite(data.table(as.table(round(J,4))), "output/TableSx_jaccard_matrix.tsv", sep="\t")

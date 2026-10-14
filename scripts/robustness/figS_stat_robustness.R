suppressMessages({library(fgsea);library(data.table);library(dplyr);library(ggplot2)})
lab3 <- c("MGAT3","FUCA1","ST3GAL/ST6GAL/CMAS")

## ---- Panel A: size-stratified (curveball degree-preserving) null vs fgsea ----
nn <- as.data.table(read.delim("output/gsea_network_null.tsv", check.names=FALSE))
cat("network-null cols:", paste(colnames(nn),collapse=" | "),"\n")
nn[, dir := ifelse(fgsea_NES>0,"increased (red)","lost (blue)")]
nn[, lab := ifelse(pathway %in% lab3, pathway, NA_character_)]
pA <- ggplot(nn, aes(-log10(fgsea_pval), -log10(netnull_pval))) +
  geom_abline(slope=1,intercept=0,linetype=2,color="grey60") +
  geom_hline(yintercept=-log10(0.05),linetype=3,color="grey70")+
  geom_vline(xintercept=-log10(0.05),linetype=3,color="grey70")+
  geom_point(aes(size=size,color=dir),alpha=0.8) +
  geom_text(aes(label=lab),size=3,vjust=-0.9,na.rm=TRUE)+
  scale_color_manual(values=c("increased (red)"="#C0392B","lost (blue)"="#2E86C1"))+
  scale_size_continuous(range=c(2,7))+
  labs(x=expression(-log[10]~"p (fgsea, score permutation)"),
       y=expression(-log[10]~"p (curveball degree/size-preserving null)"),
       size="set size", color="direction",
       title="A  Size-stratified null (curveball) vs fgsea")+
  theme_bw(base_size=11)+theme(legend.position="right")

## ---- Panel B: permutation-number sensitivity ----
gmt <- readRDS("output/glycan2enzMap_enzyme_grouped.rds")
grp <- as.data.table(readRDS("output/glycomeRes_group.rds"))
lfc <- setNames(grp$logFC,grp$group); lfc<-sort(lfc)
nperms<-c(1000,5000,10000,50000,100000); seeds<-1:5
R <- rbindlist(lapply(nperms,function(np) rbindlist(lapply(seeds,function(s){
  set.seed(s); r<-suppressWarnings(fgsea(gmt,lfc,minSize=1,maxSize=500,nPermSimple=np))
  data.table(pathway=r$pathway,pval=r$pval,nperm=np,seed=s)}))))
focus<-c("MGAT3","FUCA1","ST3GAL/ST6GAL/CMAS","FUT4/7/9","FUT3/5/6","B3GALT5")
S <- R[pathway %in% focus, .(mean_p=mean(pval),sd_p=sd(pval)),by=.(pathway,nperm)]
fwrite(S,"output/FigureS_nperm_sensitivity.tsv",sep="\t")
pB <- ggplot(S, aes(nperm, mean_p, color=pathway))+
  geom_hline(yintercept=0.05,linetype=2,color="grey50")+
  geom_line()+geom_point(size=1.6)+
  geom_errorbar(aes(ymin=pmax(mean_p-sd_p,1e-5),ymax=mean_p+sd_p),width=0.06)+
  scale_x_log10(breaks=nperms)+scale_y_log10()+
  labs(x="nPermSimple (log scale)",y="p-value (mean ± sd, 5 seeds)",
       color="enzyme group", title="B  Permutation-number sensitivity")+
  theme_bw(base_size=11)+theme(legend.position="right")

## ---- combine ----
ok_pw <- requireNamespace("patchwork", quietly=TRUE)
ok_ge <- requireNamespace("gridExtra", quietly=TRUE)
pdf("output/FigureS_stat_robustness.pdf", width=13, height=5)
if(ok_pw){ library(patchwork); print(pA + pB + patchwork::plot_layout(widths=c(1,1))) } else if(ok_ge){ gridExtra::grid.arrange(pA,pB,ncol=2) } else { print(pA); print(pB) }
dev.off()
cat("WROTE output/FigureS_stat_robustness.pdf\n")

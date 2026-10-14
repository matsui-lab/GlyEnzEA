suppressMessages({library(fgsea);library(data.table);library(ggplot2)}); set.seed(1)
gmt <- readRDS("output/glycan2enzMap_enzyme_grouped.rds")
grp <- as.data.table(readRDS("output/glycomeRes_group.rds"))
lfc <- setNames(grp$logFC, grp$group); lfc <- sort(lfc, decreasing=TRUE)  # fgsea convention
G <- names(lfc); enz <- names(gmt)
M <- matrix(0L,length(enz),length(G),dimnames=list(enz,G))
for(e in enz) M[e, intersect(gmt[[e]],G)] <- 1L
es_one <- function(members){ idx<-which(G %in% members); if(length(idx)==0) return(NA_real_); fgsea::calcGseaStat(stats=lfc, selectedStats=idx) }
obs <- es_one(gmt[["MGAT3"]]); k <- length(gmt[["MGAT3"]])

## Null A: score permutation = random size-k glycan sets
nullA <- replicate(20000, es_one(sample(G,k)))
## Null B (reviewer 1): enzyme-LABEL permutation (edges fixed; MGAT3 label -> a random existing enzyme set)
nullB <- sapply(gmt, es_one)
## Null C (reviewer 2): curveball degree-preserving rewiring
curveball <- function(M, nt=5L*sum(M)){
  R<-lapply(seq_len(nrow(M)),function(i) which(M[i,]==1L))
  for(t in seq_len(nt)){ ij<-sample(length(R),2); A<-R[[ij[1]]]; Bp<-R[[ij[2]]]
    cm<-intersect(A,Bp); oa<-setdiff(A,cm); ob<-setdiff(Bp,cm); pool<-c(oa,ob)
    if(!length(pool)) next; na<-sample(pool,length(oa)); R[[ij[1]]]<-c(cm,na); R[[ij[2]]]<-c(cm,setdiff(pool,na)) }
  Mb<-matrix(0L,nrow(M),ncol(M),dimnames=dimnames(M)); for(i in seq_along(R)) Mb[i,R[[i]]]<-1L; Mb }
mi<-which(enz=="MGAT3"); nullC<-numeric(1000)
for(b in 1:1000){ Mb<-curveball(M); nullC[b]<-es_one(G[Mb[mi,]==1L]) }

pA<-(sum(nullA<=obs,na.rm=T)+1)/(length(nullA)+1)
pB<-sum(nullB<=obs,na.rm=T)/length(nullB)
pC<-(sum(nullC<=obs,na.rm=T)+1)/(length(nullC)+1)
cat(sprintf("\nobserved MGAT3 ES = %.3f (set size %d)\n",obs,k))
cat(sprintf("Null A  score-permutation (random %d-glycan sets) : p = %.4f\n",k,pA))
cat(sprintf("Null B  enzyme-label permutation [reviewer (1)]   : p = %.4f  (MGAT3 rank %d / %d enzyme groups)\n",pB,sum(nullB<=obs),length(nullB)))
cat(sprintf("Null C  curveball degree-preserving [reviewer (2)]: p = %.4f\n",pC))

df <- rbind(data.table(ES=nullA,null="A: score-perm (random sets)"),
            data.table(ES=nullC,null="C: curveball (degree-preserving)"))
pdf("output/FigureS_three_null_MGAT3.pdf", width=8, height=5)
print(ggplot(df,aes(ES,fill=null))+geom_density(alpha=0.45,color=NA)+
  geom_vline(xintercept=obs,linetype=2,linewidth=0.8)+
  geom_rug(data=data.table(ES=as.numeric(nullB),null="B: enzyme-label perm (real sets)"),aes(color=null),inherit.aes=FALSE,linewidth=0.8)+
  annotate("text",x=obs,y=0,label=sprintf("observed MGAT3 ES=%.2f",obs),hjust=-0.05,vjust=-0.5,size=3)+
  scale_fill_manual(values=c("A: score-perm (random sets)"="#7FB3D5","C: curveball (degree-preserving)"="#F5B041"))+
  scale_color_manual(values=c("B: enzyme-label perm (real sets)"="#C0392B"))+
  labs(x="MGAT3 enrichment score (ES)",y="null density",fill="null model",color="",
       title="Three null models for MGAT3 (lost-direction)",
       subtitle=sprintf("p: score-perm=%.3f | label-perm=%.3f | curveball=%.3f",pA,pB,pC))+
  theme_bw(base_size=11))
dev.off()
fwrite(data.table(null=c("A_score_perm","B_label_perm","C_curveball"),p=c(pA,pB,pC),obs_ES=obs),
       "output/FigureS_three_null_MGAT3.tsv",sep="\t")
cat("\nWROTE output/FigureS_three_null_MGAT3.pdf / .tsv\n")

library(seqinr)
library(phyclust)
library("PerformanceAnalytics")
set.seed(123)

# Generate a tree
tree.ms <- read.nexus("/Users/chewyiting/Desktop/Mycobacterium/output/simtophylo_out/shortening/YT01_p10_s3.newick")
write.tree(tree.ms,file="/Users/chewyiting/Desktop/Mycobacterium/output/simtophylo_out/shortening/YT01_p10_s3_actual.newick")
# Generate nucleotide sequences
# (anc.HKY <- rep(0:3, 3))
# paste(nid2code(anc.HKY, lower.case = FALSE), collapse = "")

genome_length <- 4411532
seq_length <- 1000

myMutRate <- function(orgMutRate, lenAln, lenGenome, branchScale){
  clock_SNP <- orgMutRate * (lenGenome/lenAln)
  clock_SNP_days <- branchScale * clock_SNP
  return(clock_SNP_days)
}

mut_rate <- myMutRate(2.5E-7,seq_length, genome_length, 1/365)
bases <- c("A","C","G","T")

reference1 <- sample(bases, seq_length, prob =  c(0.1734,0.3268,0.3263,0.1735), replace = TRUE)
write(reference1,"myrefgenome.txt",sep="",ncolumns=1000)
anc.HKY <- code2nid(reference1)

pi.HKY <- c(0.1734,0.3263,0.3268,0.1735) # A,G,C,T

fconstCalc <- function(genLen, alnLen, piFreqs){
  nConst <- genLen - alnLen
  pi.IQTREE <- c(pi.HKY[1],pi.HKY[3],pi.HKY[2],pi.HKY[4]) # A,C,G,T
  fConst <- round(nConst*pi.IQTREE)
  return(fConst)
}

fconst <- fconstCalc(genome_length,seq_length,pi.HKY)

kappa <- 4.217

HKY.1 <- gen.seq.HKY(tree.ms, pi.HKY, kappa, seq_length, rate.scale = mut_rate,
                      anc.seq = anc.HKY)

# if rate scale = 1, 1 sub/site/


ret <- read.seqgen(HKY.1)

write.fasta(ret$org, "myseqgendummy.fasta", classid = NULL, seqname = NULL,
            width.line = 60, lower.case = FALSE, code.type = .code.type[1],
            sep = "")


#############################################
##  Comparing ttree to ptree               ##
#############################################
# Do branch lengths of ptree and ttree correlate?
# If yes, do they correlate with lifetime of lineages wrt sampling time ?
load("/Users/chewyiting/Desktop/Mycobacterium/output/simtophylo_out/shortening/YT_01_p01_s1.Rdata")
rm(mystateDF,endemicStage,endStageTime,IndexIn,obsHeight,point_lagtime,point_prevspec,
   RootNode,stree_start,stree_stop,tree3_edge,lifeHistory,rewriteEdgeMatrix,tipSampling,treeSampling)

myHistory1 <- myHistory %>% filter(tips %in% sample_tips)
myHistory2 <- myHistory1 %>% mutate(life=tsample-births) %>% select(-c(births,deaths))

ptree <- read.tree("/Users/chewyiting/Desktop/Mycobacterium/output/simtophylo_out/phylotrees/YT01_p07_s1.fasta.treefile")
ttree <- read.tree("/Users/chewyiting/Desktop/Mycobacterium/output/simtophylo_out/shortening/YT01_p07_s1_actual.newick")


sampled_tips <- sample(ptree$tip.label,200,replace=FALSE)

treeSampling <- function(sampledTips, wholeTree){
  tips <- wholeTree$tip.label
  tips <- tips[! tips %in% sampledTips]
  s_tree <- drop.tip(wholeTree, tips)
  return(s_tree)
}

ptree_sampled <- treeSampling(sampled_tips,ptree)
ttree_sampled <- treeSampling(sampled_tips,ttree)


#creation of the association matrix:
association <- cbind(ttree_sampled$tip.label, ptree_sampled$tip.label)

cophyloplot(ttree_sampled, ptree_sampled, assoc = association,
            length.line = 100, space = 150, gap = 10,show.tip.label = F)

cophyloplot(ptree01, ptree10, assoc = NULL,
            show.tip.label = F)

#plot with rotations
## Not run: 
# cophyloplot(tree1, tree2, assoc=association, length.line=4, space=28, gap=3, rotate=TRUE)

ptreeBrlens<- setNames(ptree$edge.length[sapply(1:length(ptree$tip.label),
                                  function(x,y) which (y==x),y=ptree$edge[,2])],ptree$tip.label)
ttreeBrlens<- setNames(ttree$edge.length[sapply(1:length(ttree$tip.label),
                                                function(x,y) which (y==x),y=ttree$edge[,2])],ttree$tip.label)

ptreeBlens <- cbind(as.integer(names(ptreeBrlens)),unname(ptreeBrlens))
ttreeBlens <- cbind(as.integer(names(ttreeBrlens)),unname(ttreeBrlens))

reorderpTreeLens <- ptreeBlens[match(ttreeBlens[,1],ptreeBlens[,1]),2]
reorderHistory <- myHistory2[match(ttreeBlens[,1],myHistory2$tips),2]
compareBrlens <- data.frame(id=ttreeBlens[,1],
                            History=reorderHistory,
                            ttreeLen=ttreeBlens[,2],
                            ptreeLen=reorderpTreeLens)
rm(ptreeBrlens,ttreeBrlens,ptreeBlens,ttreeBlens)

corcoff <- cor(compareBrlens$ttreeLen,compareBrlens$ptreeLen)

chart.Correlation(compareBrlens[,2:4], histogram=TRUE, pch=19)

q <- ggplot(compareBrlens) +
  geom_point(aes(x=id,y=ttreeLen))
p <- ggplot(compareBrlens) +
  geom_point(aes(x=id,y=ptreeLen))
r <- ggplot(compareBrlens) +
  geom_point(aes(x=id,y=History))
library(gridExtra)
grid.arrange(p,q,r)


stree_stat <- read.csv("/Users/chewyiting/Desktop/Mycobacterium/output/stree_asymmetree.csv")

qC <- ggplot(data=stree_stat,aes(x=Rhc,y=Norm_cherries)) +
  geom_line(aes(colour=factor(prevalence))) +
  geom_hline(yintercept=0.33333333333) +
  theme_bw()+
  theme(legend.position="none")
qK <- ggplot(data=stree_stat,aes(x=Rhc,y=Norm_sackin)) +
  geom_line(aes(colour=factor(prevalence))) +
  geom_hline(yintercept=0) +
  theme_bw()+
  theme(legend.position="none")
grid.arrange(qC,qK,ncol=2,nrow=1)

#########################################
## 
#########################################
tree.ms <- read.tree("/Users/chewyiting/Desktop/Mycobacterium/output/checkScaleTree.newick")
plot(tree.ms,show.tip.label = F)
axisPhylo(backward=F)
genome_length <- 4411532
seq_length <- 1000

myMutRate <- function(orgMutRate, lenAln, lenGenome, branchScale){
  clock_SNP <- orgMutRate * (lenGenome/lenAln)
  clock_SNP_days <- branchScale * clock_SNP
  return(clock_SNP_days)
}

mut_rate <- myMutRate(1.3E-3, seq_length, genome_length, 1/365)
bases <- c("A","C","G","T")

reference1 <- sample(bases, seq_length, prob =  c(0.1734,0.3268,0.3263,0.1735), replace = TRUE)
write(reference1,"myrefgenome.txt",sep="",ncolumns=1000)

##########################################
# Check scaling factor accounts for days -> years conversion of branch lengths and for genome --> aln conversion of fasta
##########################################
dummytree<- rtree(200,rooted=TRUE)
plot(dummytree)
axisPhylo(backward=F)
dummytree2 <- dummytree
dummytree2$edge.length <- dummytree$edge.length*365
plot(dummytree2)
axisPhylo(backward=F)
genome_length <- 100000
seq_length <- 1000

myMutRate <- function(orgMutRate, lenAln, lenGenome, branchScale){
  clock_SNP <- orgMutRate * (lenGenome/lenAln)
  clock_SNP_days <- branchScale * clock_SNP
  return(clock_SNP_days)
}

mut_rate <- myMutRate(1.3E-3,seq_length, genome_length, 1/365)
bases <- c("A","C","G","T")
reference1 <- sample(bases, seq_length, prob =  c(0.1734,0.3268,0.3263,0.1735), replace = TRUE)
write(reference1,"dummyscale.txt",sep="",ncolumns=seq_length)
write.tree(dummytree2,'dummyscale.newick')
dummytree2 <- read.tree('dummyscale.newick')
dates <- cbind(dummytree2$tip.label,1/365*dist.nodes(dummytree2)[201,1:200]) # sampling time given in years and with respect to age of MRCA of sampled lineages.
write(t(dates),file="dummyscale_dates.txt",ncolumns=2,sep = "\t")
dates <- read_tsv("dummyscale_dates.txt",col_names = F)


# tree1 <- rtree(40)
# tree2 <- rtree(20)
# 
# #creation of the association matrix:
# association <- cbind(tree2$tip.label, tree2$tip.label)
# 
# cophyloplot(tree1, tree2, assoc = association,
#             length.line = 4, space = 28, gap = 3)

################################
##  Load libraries            ##
################################
library("apTreeshape") # Must be installed for treeImbalance
library("treeImbalance")
library("adephylo")
library("phylobase")
library("ape")
library("cluster")
library("phytools")
library("dplyr")

load("simtophylo_p02.Rdata")
tree <- myphylo_out
rm(myphylo_out)

# removing variables we don't need
rm(myedges_out,myepi_out,myreactionDF,myreactions,myrxns_cleaned_out,myrxns_out,reactions,rname_all)
rm(branches,brlengths,edges,fetcher,fetcher2,imputeRxns,myrxns_cleaned,newbranches,paringDown,phylotonewick,shortenSim,reorder_rxns)

# finding timepoint when a specific prevalence was reached
point_prevspec <- 0.4
point_lagtime <- 10*365 # Adding 10 years
endemicStage <- min(which(mystateDF$Ic>point_prevspec*(mystateDF$Ic+mystateDF$Sc)))
endStageTime <- mystateDF[endemicStage,]$time + point_lagtime

#######################################
##  Getting life history for tips    ##
#######################################
lifeHistory <- function(branches){
  n_tips <- max(sort(unique(branches$infector)))
  birthtimes <- branches %>% 
    filter(reaction != "D.Ic" & reaction != "D.Ih" & reaction !="endSim") %>%
    arrange(infectee)
  deathtimes <- branches %>% 
    filter(reaction == "D.Ic" | reaction == "D.Ih") %>% 
    arrange(infector)
  
  # Assign birth and death times
  life_df <- data.frame(tips=seq(1,n_tips,1), births=rep(0,n_tips), deaths=rep(max(mybranches_out$time),n_tips))
  life_df$births[2:nrow(life_df)] <- birthtimes$time
  life_df$deaths[match(deathtimes$infector,life_df$tips)] <- deathtimes$time
  return(life_df)
}

myHistory <- lifeHistory(mybranches_out)

#######################################
##  Sample transmission tree         ##
#######################################
# Include tips from whole transmission tree only if sample is alive at time of sampling (tsample)

tsample <- endStageTime

tipSampling <- function(tsample,life_df,fsampled){
  all_tipDF <- life_df %>% filter(births < tsample & deaths > tsample)
  all_tips <- all_tipDF$tips
  nsampled <- round(fsampled*length(all_tips),digits=0)
  s_tips <- sort(sample(all_tips,nsampled,replace=F))
  return(s_tips)
}

sample_tips <- tipSampling(tsample,myHistory,1)
stree_start <- min(myHistory %>% filter(tips %in% sample_tips) %>% select(births)) # Earliest birth
stree_stop <- max(myHistory %>% filter(tips %in% sample_tips) %>% select(deaths)) # Latest death

# Will sampled tree include index? 
IndexIn <- 1 %in% sample_tips
orgRoot <- tree$root.edge
tree$root.edge <- 0 # Remove root edge from whole transmission tree to be used as input for treeSampling

treeSampling <- function(sampledTips, wholeTree){
  tips <- c(1:length(wholeTree$tip.label))
  tips <- tips[! tips %in% sampledTips]
  s_tree <- drop.tip(wholeTree, tips)
  return(s_tree)
}

tree2 <- treeSampling(sample_tips,tree)
rm(tree)

RootNode <- length(sample_tips) +1 
obsHeight <- max(dist.nodes(tree2)[RootNode,])
banana <- stree_stop-orgRoot-obsHeight

#######################################
##  Shorten sampled tree             ##
#######################################
# rewrite terminal branches to reflect sampling time
# because branch lengths represent amount of time a lineage can accumulate substitutions for, 
# when this shortened sampled transmission tree is used as input for simulating sequence evolution

rewriteEdgeMatrix <- function(my_tree,tsample){
  tiplabels <- as.numeric(my_tree$tip.label)
  tiplabels_old <- seq(1,length(tiplabels),1)
  new_edges <- my_tree$edge.length
   
  # find full node path for each tip and its whole length
  
  for(j in 1:length(tiplabels)){
    # Full node path and length, from root to tip
    fullpath <- nodepath(my_tree,RootNode,tiplabels_old[j])
    fullLength <- dist.nodes(my_tree)[RootNode,tiplabels_old[j]]
    deduct <- fullLength-(tsample-orgRoot-banana)
      
    # Finding terminal branch
    node_indices <- tail(fullpath,2)
    tip_term <- which(my_tree$edge[,1]==node_indices[1] & my_tree$edge[,2]==node_indices[2])
    old.term.length <- my_tree$edge.length[tip_term] # terminal branch length 
    new.term.length <- old.term.length - deduct 
    
    # Update terminal branch length
    new_edges[tip_term] <- new.term.length
    print(paste(c("tip=",j,"of",length(tiplabels))))
  }
  return(new_edges)
}

tree3_edge <- rewriteEdgeMatrix(tree2,tsample)
tree3 <- tree2 
tree3$edge.length <- tree3_edge

ape::write.nexus(tree3,file="YT01_p02_s3.newick")
save.image(file="YT_01_p02_s3.Rdata")

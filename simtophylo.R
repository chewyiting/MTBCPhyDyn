#########################################################
##     Assigning inheritance relationships             ##
#########################################################
load.image("epi_history.Rdata")

# 01. Defining removal of infectors
# deathp assigns probabilities of death (D.Ic or D.Ih), to be linearly proportional to time of birth. 

deathp <- function(times){
  if(length(times)>1){
    p.death <- times
    p.death = 0.5*(1 - (times/sum(times)))
    
  } else {
    p.death <- 1
    
  }
  return(p.death)
}
# timevec <- c(1,39,95)
# plot(timevec,deathp(timevec))

# 02. Assigning infector-infectee relationships
# Each transmission event is assigned an infector-infectee pair
# Each removal event is assigned an infector, no infectee

branches <- function(reactionsDF,reactions){
  # Initialize state
  outDF <- reactionsDF %>% mutate(infector=NA,infectee=NA)
  Ic_vec <- seq(1,(reactions$nXcc + reactions$nXhc + 1),1)
  Ih.start <- length(Ic_vec)+1
  Ih.stop <- length(Ic_vec) + reactions$nXhh + reactions$nXch
  Ih_vec <- c(Ih.start:Ih.stop)
  Ic_avail <- 1
  Ih_avail <- head(Ih_vec,1)
  Ic_born <- 0
  Ih_born <- 0
  i <- 2
  j <- 1
  r <- 1
  nrxns <- nrow(reactionsDF)
  # Check reaction
  for(r in 1:nrxns) {
    reaction <- reactionsDF[r,]$reaction
    reaction.id <- which(row.names(rname_all)==reaction)
    reaction.time <- reactionsDF[r,]$time
    
    if(reaction.id ==3){ # Xcc (cattle-cattle transmission)
      #Choose infector & infectee
      chosen.one <- sample(length(Ic_avail),1)
      Ic.infector <- Ic_avail[chosen.one]
      Ic.infectee <- Ic_vec[i]
      #Update 
      outDF[r,]$infector <- Ic.infector
      outDF[r,]$infectee <- Ic.infectee
      Ic_avail <- c(Ic_avail,Ic.infectee)
      
      newIcind <- which(Ic_avail==Ic.infectee)
      Ic_born[newIcind] <- reaction.time
      
      i <- i + 1
    } else if(reaction.id==4) { # Xch (cattle->human transmission)
      # Choose infector & infectee
      chosen.one <- sample(length(Ic_avail),1)
      Ic.infector <- Ic_avail[chosen.one]
      Ih.infectee <- Ih_vec[j]
      # Update
      outDF[r,]$infector <- Ic.infector
      outDF[r,]$infectee <- Ih.infectee
      Ih_avail <- unique(c(Ih_avail,Ih.infectee))
      
      newIhind <- which(Ih_avail==Ih.infectee)
      Ih_born[newIhind] <- reaction.time
      print(paste(c(Ih_avail,"/",Ih_born)))
      j <- j + 1
    } else if(reaction.id==5){ #Xhc (human->cattle transmission)
      # Choose infector & infectee
      chosen.one <- sample(length(Ih_avail),1)
      Ih.infector <- Ih_avail[chosen.one]
      Ic.infectee <- Ic_vec[i]
      # Update
      outDF[r,]$infector <- Ih.infector
      outDF[r,]$infectee <- Ic.infectee
      Ic_avail <- c(Ic_avail,Ic.infectee)
      
      newIcind <- which(Ic_avail==Ic.infectee)
      Ic_born[newIcind] <- reaction.time
      
      i <- i + 1
    } else if(reaction.id==6){ #Xhh (human-human transmission)
      # Choose infector & infectee
      chosen.one <- sample(length(Ih_avail),1)
      Ih.infector <- Ih_avail[chosen.one]
      Ih.infectee <- Ih_vec[j]
      # Update
      outDF[r,]$infector <- Ih.infector
      outDF[r,]$infectee <- Ih.infectee
      Ih_avail <- c(Ih_avail,Ih.infectee)
      
      newIhind <- which(Ih_avail==Ih.infectee)
      Ih_born[newIhind] <- reaction.time
      print(paste(c(Ih_avail,"/",Ih_born)))
      j <- j + 1
    } else if(reaction.id==1){ # D.Ic (removal of infected cow)
      # Choose infector
      Ic_deathprob <- deathp(Ic_born)
      chosen.one <- sample(length(Ic_avail),1,prob=Ic_deathprob)
      Ic.infector <- Ic_avail[chosen.one]
      # Update
      outDF[r,]$infector <- Ic.infector
      
      newindices <- !Ic_avail %in% Ic.infector
      Ic_avail <- Ic_avail[newindices]
      Ic_born <- Ic_born[newindices]
      
    } else if(reaction.id==2){ # D.Ih (removal of infected human)
      # Choose
      Ih_deathprob <- deathp(Ih_born)
      chosen.one <- sample(length(Ih_avail),1,prob=Ih_deathprob)
      Ih.infector <- Ih_avail[chosen.one]
      # Update
      outDF[r,]$infector <- Ih.infector
      newindices <- !Ih_avail %in% Ih.infector
      Ih_avail <- Ih_avail[newindices]
      Ih_born <- Ih_born[newindices]
      print(paste(c("reaction ",r,Ih_avail,"/",Ih_born)))
    } 
    r <- r + 1
  }
  return(outDF)
}

myrxns_out <- branches(myreactionDF,myreactions) 

# 03. Identifying all extant tips
# This is important for data clean up before nodes are assigned
# This is because extant tips (infectious individuals which have not been removed when simulation has terminated)
# need to be assigned a terminal branch leading from its last transmission / its inception
# this terminal branch is labeled with its corresponding reaction being "endSim" 

myrxns_cleaned <- function(mybranches){
  branchout <- mybranches
  x <- unique(c(mybranches$infectee,mybranches$infector))
  tips <- sort(x[!is.na(x)])
  
  death_vec <- mybranches %>% filter(reaction=="D.Ic" | reaction=="D.Ih")
  annoyingtips <- tips[!tips %in% death_vec$infector]
  slice <- data.frame(reaction="endSim",rownumber=NA,
                      time=max(branchout$time),infector=annoyingtips,infectee=NA)
  branchout <- rbind(branchout,slice)
  return(branchout)
}

myrxns_cleaned_out <- myrxns_cleaned(myrxns_out)

# 04. Assigning parental and child nodes defining each branch 
# outputDF stores these nodes in the 'from' and 'to' columns respectively.

edges <- function(mybranches){
  # Initialize 
  edgelist <- mybranches %>% mutate(from=NA,to=NA)
  x <- unique(c(mybranches$infectee,mybranches$infector))
  tips <- sort(x[!is.na(x)])
  y1 <- tail(tips,1) + 1
  y2 <- tail(tips,1) + length(tips) 
  node_vec <- c(y1:y2)
  rm(x,y1,y2)
  
  nodepaths <- vector("list",length(tips))
  # Finding infectors, 
  infectorsDF <- mybranches %>% filter(reaction!="D.Ic" & reaction!="D.Ih" & reaction!="endSim")
  infector_vec <- unique(infectorsDF$infector)
  infectee_vec <- tips[!tips %in% infector_vec]
  for(t in 1:length(infector_vec)){
    t <- infector_vec[t]
    print(paste(c("infector=",t)))
    
    # Find nodepath for each infector
    if(is.null(nodepaths[[t]])){ # this is only true for the first infected individual.
      nrxns <- nrow(mybranches %>% filter(infector==t))
      nodepath <- c(node_vec[1:nrxns],t)
      
    } else { 
      nodestart <- nodepaths[[t]][1]
      nrxns <- nrow(mybranches %>% filter(infector==t))
      
      if(nrxns==1){
        nodepath <- c(nodestart,t)
        
      } else {
        nextnode <- max(unlist(nodepaths,`[[`,1)) + 1 # first unassigned node
        j <- which(node_vec==nextnode)
        k <- j + nrxns -2
        nodepath <- c(nodestart,node_vec[j:k],t)
        
      }
    }
    # Update nodepaths of tip
    nodepaths[[t]] <- nodepath
    
    # Update edge list
    indices <- which(edgelist$infector==t)
    edgelist[indices,]$from <- head(nodepath,-1)
    edgelist[indices,]$to <- tail(nodepath,-1)
    
    # Update nodepaths of infectees
    infectees <- sort(na.exclude(edgelist[indices,]$infectee))
    if(length(infectees) != 0){
      for(i in 1:length(infectees)){
        in.fectee <- infectees[i]
        slice <- edgelist %>% filter(infectee==in.fectee)
        nodepaths[[in.fectee]] <- slice$to
      }
    }
    print(nodepath)
  }
  
  
  for(l in 1:length(infectee_vec)){
    theinfectee <- infectee_vec[l]
    # print(nodepaths[[theinfectee]])
    nodestart <- nodepaths[[theinfectee]][1]
    nodepath <- c(nodestart,theinfectee)
    nodepaths[[theinfectee]] <- nodepath
    # Update edge list
    indices <- which(edgelist$infector==theinfectee)
    edgelist[indices,]$from <- head(nodepath,-1)
    edgelist[indices,]$to <- tail(nodepath,-1)
  }
  return(edgelist)
}

myedges_out <- edges(myrxns_cleaned_out)

# 05. Calculating branch lengths 
brlengths <- function(myedges){
  mybranches <- myedges %>%mutate(edge.list=NA)
  
  # Preparing root node's branch length
  rootnode <- min(myedges$from)
  rootind <- which(mybranches$from==rootnode)
  mybranches[rootind,]$edge.list <- mybranches[rootind,]$time
  allinds <- seq(1,nrow(mybranches),1)
  otherind <- allinds[!allinds %in% rootind]
  nbranches <- length(otherind) 
  
  # For all other branch lengths..
  for(b in 1:nbranches){
    branchind <- otherind[b]
    branch.pNode <- myedges[branchind,]$from
    branch.pTime <- myedges$time[which(myedges$to==branch.pNode)]
    branch.cTime <- myedges[branchind,]$time
    mybranches[branchind,]$edge.list <- branch.cTime - branch.pTime
  }
  return(mybranches)
}

mybranches_out <- brlengths(myedges_out)

# 06. Reordering edge list suitable for phylo 

newbranches <- function(mybrlengths){ # mybranches_out = mybrlengths
  tips <- sort(unique(mybrlengths$infector))
  Nnode <-  length(tips) -1
  node_vec <- (seq(1,Nnode,1))+max(tips)
  root1 <- node_vec[1]
  root2 <- node_vec[2]
  unroot <- mybrlengths %>% filter(mybrlengths$from!=root1 & mybrlengths$to!=root2)
  newfrom <- unroot$from - 1
  newto_ind <- which(unroot$to > max(tips)) 
  newto <- unroot$to
  newto[newto_ind] <- unroot$to[newto_ind] - 1 
  unroot$from <- newfrom
  unroot$to <- newto
  return(unroot)
}

myepi_out <- newbranches(mybranches_out)

# 07. Obtaining phylo object
phylotonewick <- function(epiDF,mybranches){ 
  x <- matrix(c(epiDF$from,epiDF$to),ncol=2)
  myroot <- mybranches[1,]$edge.list 
  tips <- sort(unique(epiDF$infector))
  numnodes <- length(tips) -1 
  tr <- list(edge=x, tip.label=tips,
              Nnode= numnodes, edge.length=epiDF$edge.list,
              root.edge=myroot)
  class(tr) <- "phylo"
  return(tr)
}

myphylo_out <- phylotonewick(myepi_out, mybranches_out)
plot(myphylo_out)

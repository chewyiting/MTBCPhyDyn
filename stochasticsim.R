#########################################################
##      Load libraries                                 ##
#########################################################
library(dplyr)
library(rootSolve)
library(adaptivetau) 
library(compiler)

#########################################################
##      Parameter space                                ##
#########################################################

RhcFinder <- function(Ic,Rr,Ib,Rcc){
  Rhc <- (Ic*(1+Rcc*(Ic-1)))/(Ib*Rr*(1-Ic))
  return(Rhc)
}

xRcc <- seq(0.01,5/3,length.out=10)
xRhc <- RhcFinder(0.4,10/3,3E-3,xRcc)

params.setR <- cbind(
  Rhh = 0.5,
  Rcc = xRcc,
  Rhc = xRhc,
  Rch = 0.01257523
)

xmuC <- 0.1/365 #(10 year lifespan)
xmuH <- 0.3333333/365 #1/3

params.setB <- as.data.frame(params.setR) %>%
  mutate(muC=xmuC,muH=xmuH)%>%
  mutate(Bhh=(Rhh*xmuH),
         Bcc=(Rcc*xmuC),
         Bhc=(Rhc*xmuH),
         Bch=(Rch*xmuC)) %>%
  select(muC,muH,Bhh,Bcc,Bhc,Bch)

#########################################################
##      Deterministic behavior                         ##
#########################################################
# Run this section to check system's deterministic behavior

# SIS_dydt_days <- function(time,y,parms){
  Sc <- y[1]
  Sh <- y[2]
  Ic <- y[3]
  Ih <- y[4]

  muC <- parms[1]
  muH <- parms[2]
  Bhh <- parms[3]
  Bcc <- parms[4]
  Bhc <- parms[5]
  Bch <- parms[6]

  #Differential Equations
  dSc <- muC*201 - Bcc*Sc*(Ic/201) - Bhc*Sc*(Ih/201)-muC*Sc
  dSh <- muH*201 - Bhh*Sh*(Ih/201) - Bch*Sh*(Ic/201)-muH*Sh
  dIc <- Bcc*Sc*(Ic/201) + Bhc*Sc*(Ih/201) - muC*Ic
  dIh <- Bhh*Sh*(Ih/201) + Bch*Sh*(Ic/201) - muH*Ih
  return(list(c(dSc,dSh,dIc,dIh)))
}

# dummy <- runsteady(c(Sc=200,Sh=201, Ic=1,Ih=0),
                      func = SIS_dydt_days, parms = unlist(params.setB[1,]), times = c(0,100*365))
                      
                      
#########################################################
##     Stochastic simulation                           ##
#########################################################
# Define state transitions and demivariables recording each transition

transitions = list(c(Sc=+1,B.Sc=+1), #birth.Sc
                   c(Sc=-1,D.Sc=+1), #death.Sc
                   c(Ic=-1,D.Ic=+1), #death.Ic
                   c(Sh=+1,B.Sh=+1), #birth.Sh
                   c(Sh=-1,D.Sh=+1), #death.Sh
                   c(Ih=-1,D.Ih=+1), #death.Ih
                   c(Sc=-1,Ic=+1,Xcc=+1), # Cattle-cattle transmission events recorded by Xcc
                   c(Sc=-1,Ic=+1,Xhc=+1), # Human->cattle transmission events recorded by Xhc
                   c(Sh=-1,Ih=+1,Xch=+1), # Cattle->human transmission events recorded by Xch
                   c(Sh=-1,Ih=+1,Xhh=+1) # Human-human transmission events recorded by Xhh
                   )

# Define probabilities for each transition
probabilities <- function(x,p,t){
  if(x["Ic"]!=0){ 
    return(c(
      (sum(x["Sc"],x["Ic"]))*(p$muC),                 #birth.Sc
      p$muC*x["Sc"],                                  #death.Sc
      p$muC*x["Ic"],                                  #death.Ic
      (sum(x["Sh"],x["Ih"]))*(p$muH),                 #birth.Sh
      p$muH*x["Sh"],                                  #death.Sh
      p$muH*x["Ih"],                                  #death.Ih
      p$Bcc*x["Sc"]*(x["Ic"]/sum(x["Sc"],x["Ic"])),   #cattle-cattle
      p$Bhc*x["Sc"]*(x["Ih"]/sum(x["Sh"],x["Ih"])),   #human-cattle
      p$Bch*x["Sh"]*(x["Ic"]/sum(x["Sc"],x["Ic"])),   #cattle-human
      p$Bhh*x["Sh"]*(x["Ih"]/sum(x["Sh"],x["Ih"]))    #human-human
    ))} else { #If Ic = 0, halt simulation by making all rates = 0
      return(c(0,0,0,0,0,0,0,0,0,0))
    } 
}

# Running stochastic simulations
# Each parameter combination is run for ntrials number of simulations

ntrials <- 3 
tmax <- 100*365 # 100 years

twoSImodel <- function(var1,var2){
  n.cond1 <- length(var1)
  l <- vector("list",n.cond1*ntrials)
  dim(l) <- c(n.cond1,ntrials)
  for (Rhc_ind in 1:n.cond1){
    init.values=c(Sc=200,Ic=1,Sh=201,Ih=0,
                  B.Sc=0,D.Sc=0,D.Ic=0,
                  B.Sh=0,D.Sh=0,D.Ih=0,
                  Xcc=0,Xch=0,Xhc=0,Xhh=0)
    print(paste(c("Bhc=",var1[Rhc_ind])))
          params=list(muC=params.setB$muC[1],muH=params.setB$muH[1],
                      Bcc=var2[Rhc_ind], Bhc=var1[Rhc_ind],
                      Bch=params.setB$Bch[1],Bhh=params.setB$Bhh[1])
          for(i in 1:ntrials){
            l[Rhc_ind,i] <- list(ssa.adaptivetau(
              init.values,transitions,probabilities,params,
              tf=tmax,tl.params=list(maxtau=7))) # Max tau limited to 1 week, to capture small changes.
          }
  }
  l
}
enableJIT(1) 
system.time(SI_ccs.1 <- twoSImodel(params.setB$Bhc,params.setB$Bcc)) 
saveRDS(SI_ccs.1,"SI_01_Rhc.rds")

#########################################################
##     Obtaining transmission history                  ##
#########################################################
# 01. Read in transmission simulation
SI_ccs.1<- readRDS("SI_01_Rhc.rds")
dummy <- SI_ccs.1[7,2][[1]] # Choosing the 2nd trial of the 7th parameter combination as a demo
dummy2 <- as.data.frame(dummy)

reactions <- data.frame(nD.Ic=max(dummy2$D.Ic),nD.Ih=max(dummy2$D.Ih),
               nXcc=max(dummy2$Xcc), nXch=max(dummy2$Xch),
               nXhc=max(dummy2$Xhc), nXhh=max(dummy2$Xhh))

# Reaction DF
rxnDF=data.frame(reaction=c(rep("D.Ic",reactions$nD.Ic),
                                rep("D.Ih",reactions$nD.Ih),
                                rep("Xcc",reactions$nXcc),
                                rep("Xhc",reactions$nXhc),
                                rep("Xch",reactions$nXch),
                                rep("Xhh",reactions$nXhh)
                                ),
                     number=c(seq(1,reactions$nD.Ic,1),
                              seq(1,reactions$nD.Ih,1),
                              seq(1,reactions$nXcc,1),
                              seq(1,reactions$nXhc,1),
                              seq(1,reactions$nXch,1),
                              seq(1,reactions$nXhh,1)),
                     rownumber=NA,
                 time=NA,
                 index=c(1:sum(reactions)))

 

# 02. Selecting only reactions we are interested in, i.e. those affecting dIc/dt and dIh/dt
# First, a dataframe specifying these reactions
rname_all <- data.frame(row.names=c("D.Ic","D.Ih",
                                    "Xcc","Xch",
                                    "Xhc","Xhh"),
                        val=c(8,11,12,13,14,15))

# Secondly, a function that fetches row corresponding to a specific occurrence of a specific reaction
# rxnInd = column numbers in dummy2 corresponding to each reaction
# rxnInd = 8 for d.Ic, 11 for D.Ih, 12 for Xcc, 13 for Xch, 14 for Xhc, 15 for Xhh
fetcher <- function(dataset, reaction, rxnInd){
  rowInd <- head(which(dataset[,reaction]==rxnInd),1)
  return(rowInd)
}

# Lastly, selecting from original rxnDF based on rname_all
paringDown <- function(dataset){
  for(r in c(1:6)){
  name.id <- rname_all[r,]
  name <- row.names(rname_all)[which(rname_all[,1]==name.id)]
  reactionNumber <- as.numeric(reactions[r])
  for(i in c(1:reactionNumber)){
    output <- fetcher(dataset,name.id,i)
    slice <- rxnDF %>% filter(reaction==name & number == i)
    rxnDF[slice$index,]$rownumber <- output
    rxnDF[slice$index,]$time <- dataset[output,]$time
  }
  }
  outDF <- rxnDF %>% select(-c(number,index)) %>% arrange(time)
return(outDF)
}

reactionsDF<- paringDown(dummy2)

# Select from original states DF corresponding to states we're interested in
# i.e., those which correspond to when interesting reactions occur

myreactionDF <- reactionsDF
mystateDF3 <- dummy2[myreactionDF$rownumber,]
mystateDF2 <- rbind(dummy2[1,],mystateDF3)
mystateDF <- mystateDF2 %>% 
  select(-c("B.Sc","B.Sh","D.Sc","D.Sh"))
myreactions <- data.frame(nD.Ic=max(mystateDF$D.Ic),nD.Ih=max(mystateDF$D.Ih),
                          nXcc=max(mystateDF$Xcc), nXch=max(mystateDF$Xch),
                          nXhc=max(mystateDF$Xhc), nXhh=max(mystateDF$Xhh))
rm(mystateDF3,mystateDF2)
rm(dummy,dummy2,reactionsDF,rxnDF,SI_ccs.1)
save.image(file="epi_history.Rdata")

#############################################
##  Load libraries                         ##
#############################################
library(simecol)
library(rootSolve)
library(dplyr)
library(compiler)

#############################################
##  Parameter space                        ##
#############################################

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

xmuC <- 0.1
xmuH <- 0.3333333 #1/3

params.setB <- as.data.frame(params.setR) %>%
  mutate(muC=xmuC,muH=xmuH)%>%
  mutate(Bhh=Rhh*xmuH,
         Bcc=Rcc*xmuC,
         Bhc=Rhc*xmuH,
         Bch=Rch*xmuC) %>%
  select(muC,muH,Bhh,Bcc,Bhc,Bch)

#############################################
##  Deterministic population dynamics      ##
#############################################

# Check state variables at different times  

# SIS_dydt <- function(time,y,parms){
#   Sc <- y[1]
#   Sh <- y[2]
#   Ic <- y[3]
#   Ih <- y[4]
#   Nc <- Sc + Ic 
#   Nh <- Sh + Ih 
# 
#   muC <- parms[1]
#   muH <- parms[2]
#   Bhh <- parms[3]
#   Bcc <- parms[4]
#   Bhc <- parms[5]
#   Bch <- parms[6]
# 
#   #Differential Equations
#   dSc <- muC*Nc - Bcc*Sc*(Ic/Nc) - Bhc*Sc*(Ih/Nh) - muC*Sc
#   dSh <- muH*Nh - Bhh*Sh*(Ih/Nh) - Bch*Sh*(Ic/Nc) - muH*Sh
#   dIc <- Bcc*Sc*(Ic/Nc) + Bhc*Sc*(Ih/Nh) - muC*Ic
#   dIh <- Bhh*Sh*(Ih/Nh) + Bch*Sh*(Ic/Nc) - muH*Ih
#   return(list(c(dSc,dSh,dIc,dIh)))
# }
# 
# dummy <- list()
# times <- c(50,75,100,200,300)
# 
# for (i in 1:5){
#   params <- unlist(params.setB[i,])
#   tlast <- times[i]
#   dummy[[i]] <- runsteady(c(Sc=200,Sb=201, Ic=1,Ib=0),
#               func = SIS_dydt, parms = params, times = c(0,tlast))
# }

#######################################################
##  Deterministic approximation of tree metrics      ##
#######################################################

reversetime <- function(d){
  newtime <- max(d$time)-d$time
  newd <- d
  newd$time <- newtime
  o <- order(newd$time)
  newd <- newd[o,]
  newd
}


tbmodel.dr <- new("odeModel",
                  main = function(time, init, parms, ...){
                    with(as.list(c(init,parms)),{
                      # ODEs
                      Nc <- Sc+Ic
                      Nh <- Sh+Ih
                      dSc <- muC*Nc - Bcc*Sc*(Ic/Nc) - Bhc*Sc*(Ih/Nh) - muC*Sc
                      dIc <- Bcc*Sc*(Ic/Nc) + Bhc*Sc*(Ih/Nh) - muC*Ic
                      dSh <- muH*Nh - Bhh*Sh*(Ih/Nh) - Bch*Sh*(Ic/Nc) - muH*Sh
                      dIh <- Bhh*Sh*(Ih/Nh) + Bch*Sh*(Ic/Nc) - muH*Ih
                      list(c(dSc,dIc,dSh,dIh))
                    })},
                  parms = list(muC=params.setB[1,]$muC, muH=params.setB[1,]$muH, 
                               Bhh=params.setB[1,]$Bhh, Bcc=params.setB[1,]$Bcc, 
                               Bch=params.setB[1,]$Bch, Bhc=params.setB[1,]$Bhc),
                  times = c(from=0,to=200,by=1/52),
                  init = c(Sc=200,Ic=1,Sh=201,Ih=0),
                  solver = "lsoda"
)
tbmodel.dr <- sim(tbmodel.dr)

# Check that population dynamics match that of deterministic model above
# plot(tbmodel.dr)

tbmodel.dr.makemat <- function(s,direction="B"){
  if(direction=="B"){
    inp <- reversetime(out(s))
  }
  else{
    inp <- out(s)
  }
  ntimes <- dim(inp)[1]
  # parameters
  Bhh <- parms(s)[["Bhh"]]
  Bcc <- parms(s)[["Bcc"]]
  Bhc <- parms(s)[["Bhc"]]
  Bch <- parms(s)[["Bch"]]
  muC <- parms(s)[["muC"]]
  muH <- parms(s)[["muH"]]
  # initialize arrays
  F.arr <- array(0.0,c(2,2,ntimes))
  G.arr <- array(0.0,c(2,2,ntimes))
  eta.arr <- array(0.0,c(1,2,ntimes))
  mu.arr <- array(0.0,c(2,1,ntimes))
  # Cycle through times
  for(i in 1:ntimes){
    Sc <- inp[i,]$Sc
    Sh <- inp[i,]$Sh
    Ic <- inp[i,]$Ic
    Ih <- inp[i,]$Ih
    Nc <- Sc+Ic
    Nh <- Sh+Ih
    F.arr[1,1,i] <- Bcc*Sc*(Ic/Nc)
    F.arr[2,1,i] <- Bhc*Sc*(Ih/Nh)
    F.arr[1,2,i] <- Bch*Sh*(Ic/Nc)
    F.arr[2,2,i] <- Bhh*Sh*(Ih/Nh)
    # G.arr is 0
    # eta.arr is 0
    mu.arr[1,1,i] <- (muC)*Ic
    mu.arr[2,1,i] <- (muH)*Ih
  }
  list(F.arr=F.arr,G.arr=G.arr,eta.arr=eta.arr,mu.arr=mu.arr,taxis=inp$time)
}
tbmodel.dr.fg <- tbmodel.dr.makemat(tbmodel.dr)

anc <- new("odeModel",
           main=function(time,init,parms,...){
             input <- approxTime(inputs,time)
             with(as.list(c(init,parms,input)),{
               dA1 <- -f11*(A1/Ic)*(A1/Ic) + # coalescence between two type 1
                 - f21*(A1/Ic)*(A2/Ih) + # coalescence with type 2 no. 1, coalescence with type 2 no. 2 doesn't give rise to net change in A1
                 + g12*(A2/Ih) + # type 2 becomes type 1 (backwards)
                 - g21*(A1/Ic) + # type 1 become type 2 (backwards)
                 - f21*(A1/Ic)*((Ih-A2)/Ih)  + # invisible transmission
                 + f12*(A2/Ih)*((Ic-A1)/Ic) # invisible transmission
               dA2 <- -f22*(A2/Ih)*(A2/Ih) + # coalescence between two type 2
                 - f12*(A2/Ih)*(A1/Ic) + # coalescence with type 1 no. 1; coalescence with type 1 no. 2 doesn't give rise to net change in A2
                 + g21*(A1/Ic) + # type 1 becomes type 2 (backwards)
                 - g12*(A2/Ih) + # type 2 become type 1 (backwards)
                 - f12*(A2/Ih)*((Ic-A1)/Ic) + # invisible transmission
                 + f21*(A1/Ic)*((Ih-A2)/Ih) # invisible transmission
               # Lineages initially in type 1 who are currently in type 1
               dL11 <- #-(L11/Ic)*(f11*L11/Ic+f12*L21/Ih+f11*L12/Ic+f12*L22/Ih) +
                 #-(L11/Ic)*(f11*L11/Ic+f21*L21/Ih+f11*L12/Ic+f21*L22/Ih) +
                 #-(L11/Ic)*(f11*B1/Ic+f12*B2/Ih)+
                 #-(L11/Ic)*(f11*B1/Ic+f21*B2/Ih)+
                 -(L11/Ic)*(2*f11*A1/Ic+(f21+f12)*A2/Ih)+
                 # Now invisible transmissions away from this state
                 -(L11/Ic)*((f21)*(Ih-A2)/Ih)+
                 # Now invisible transmissions towards this state
                 +(L21/Ih)*((f12)*(Ic-A1)/Ic)+
                 # Now migrations away from this state
                 - (L11/Ic)*g21 +
                 # Migrations towards this state
                 + (L21/Ih)*g12
               dL21 <- -(L21/Ih)*(2*f22*A2/Ih+(f12+f21)*A1/Ic)+
                 -(L21/Ih)*((f12)*(Ic-A1)/Ic)+
                 +(L11/Ic)*((f21)*(Ih-A2)/Ih)+
                 -(L21/Ih)*g12+
                 +(L11/Ic)*g21
               dL12 <- -(L12/Ic)*(2*f11*A1/Ic+(f21+f12)*A2/Ih)+
                 -(L12/Ic)*((f21)*(Ih-A2)/Ih)+
                 +(L22/Ih)*((f12)*(Ic-A1)/Ic)+
                 -(L12/Ic)*g21+
                 +(L22/Ih)*g12
               dL22 <- -(L22/Ih)*(2*f22*A2/Ih+(f12+f21)*A1/Ic)+
                 -(L22/Ih)*((f12)*(Ic-A1)/Ic)+
                 +(L12/Ic)*((f21)*(Ih-A2)/Ih)+
                 -(L22/Ih)*g12+
                 +(L12/Ic)*g21
               dC11 <- f11*(L11/Ic)^2+(f12+f21)*(L11/Ic)*(L21/Ih)+f22*(L21/Ih)^2
               dC12 <- 2*f11*(L11/Ic)*(L12/Ic)+(f12+f21)*(L11/Ic)*(L22/Ih)+(f12+f21)*(L12/Ic)*(L21/Ih)+2*f22*(L21/Ih)*(L22/Ih)
               dC22 <- f22*(L22/Ih)^2+(f12+f21)*(L12/Ic)*(L22/Ih)+f11*(L12/Ic)^2
               dX1 <- -g21*(A1/Ic)*(X1/A1)+g12*(A2/Ih)*(X2/A2)+f12*(A2/Ih)*(X2/A2)-f21*(A1/Ic)*(X1/A1)
               dX2 <-  -g12*(A2/Ih)*(X2/A2)+g21*(A1/Ic)*(X1/A1)+f21*(A1/Ic)*(X1/A1)-f12*(A2/Ih)*(X2/A2)
               dK <- 2*f11*((A1/Ic)^2)*(X1/A1)+2*f22*((A2/Ih)^2)*(X2/A2)+f12*(A1/Ic)*(A2/Ih)*(X1/A1+X2/A2)+f21*(A1/Ic)*(A2/Ih)*(X1/A1+X2/A2)
               list(c(dA1,dA2,dL11,dL21,dL12,dL22,dC11,dC12,dC22,dX1,dX2,dK))
               #list(c(dA1,dA2,dL11,dL21,dL12,dL22,dC11,dC12,dC22))
             })},
           times = c(from=0,to=200,by=1/52),
           #parms = c(X1=0,X2=0),
           parms = c(x=0),
           init = c(A1=0,A2=0,L11=0,L21=0,L12=0,L22=0,C11=0,C12=0,C22=0,X1=0,X2=0,K=0),
           #init = c(A1=0,A2=0,L11=0,L21=0,L12=0,L22=0,C11=0,C12=0,C22=0),
           inputs=NULL,
           solver = "lsoda"
)

calcclus <- function(C11,C12,C22){
  M <- matrix(c(C11,C12/2,C12/2,C22),ncol=2,nrow=2,byrow=TRUE)#try this, which splits C12 into two
  #M <- matrix(c(C11,C12/3,2*C12/3,C22),ncol=2,nrow=2,byrow=TRUE)#try this, which splits C12 into two
  mi <- rowSums(M)
  mj <- colSums(M)
  Mhat <- (matrix(rep(mi,2),ncol=2,byrow=FALSE)*matrix(rep(mj,2),ncol=2,byrow=TRUE))/sum(M)
  M-Mhat
  E <- M/sum(M)
  r <- (sum(diag(E))-sum(E%*%E))/(1-sum(E%*%E)) # check
  #ai <- colSums(E)
  #bj <- rowSums(E)
  r
}
calcpclus <- function(pvec){
  M <- matrix(c(pvec[1],pvec[2],pvec[3],pvec[4]),ncol=2,nrow=2,byrow=TRUE)
  mi <- rowSums(M)
  mj <- colSums(M)
  Mhat <- (matrix(rep(mi,2),ncol=2,byrow=FALSE)*matrix(rep(mj,2),ncol=2,byrow=TRUE))/sum(M)
  #M-Mhat
  E <- M/sum(M)
  r <- (sum(diag(E))-sum(E%*%E))/(1-sum(E%*%E)) # check
  #ai <- colSums(E)
  #bj <- rowSums(E)
  r
}

coalsim <- function(fsim,fsiminits,params,fsimtime,tstep,matfun,bsim,frac1,frac2,numseq=NULL){
  result <- list()
  fs <- fsim # simecol models forwards simulation
  parms(fs) <- params # parameter values
  init(fs) <- fsiminits # initial conditions
  times(fs) <- seq(0,fsimtime,by=tstep)
  maxindex <- length(times(fs))
  fs <- sim(fs)
  fsout <- out(fs)
  fsoutrev <- reversetime(fsout)
  mat <- matfun(fs)
  F.arr <- mat[[1]]
  G.arr <- mat[[2]]
  if(is.null(numseq)){
    # A1.init <- frac1*fsout$Ic[maxindex]
    # A2.init <- frac2*fsout$Ih[maxindex]
    A1.init <- fsout$Ic[maxindex]
    A2.init <- fsout$Ih[maxindex]
  }else{
    if(length(numseq)==1){
      Atotal <- fsout$Ic[maxindex]+fsout$Ih[maxindex]
      A1.init <- numseq*fsout$Ic[maxindex]/Atotal
      A2.init <- numseq*fsout$Ih[maxindex]/Atotal
    }else{
      A1.init <- numseq[1]
      A2.init <- numseq[2]
    }
  }
  A.init <- c(A1=A1.init,A2=A2.init,L11=A1.init,L21=0,L12=0,L22=A2.init,C11=0,C12=0,C22=0,X1=A1.init,X2=A2.init,K=0)
  #A.init <- c(A1=A1.init,A2=A2.init,L11=A1.init,L21=0,L12=0,L22=A2.init,C11=0,C12=0,C22=0)
  ancinp <- data.frame(time=seq(0,fsimtime,by=tstep),Ic=fsoutrev$Ic,Ih=fsoutrev$Ih,f11=F.arr[1,1,],f12=F.arr[1,2,],f21=F.arr[2,1,],f22=F.arr[2,2,],g11=G.arr[1,1,],g12=G.arr[1,2,],g21=G.arr[2,1,],g22=G.arr[2,2,])
  bsim <- anc
  init(bsim) <- A.init
  inputs(bsim) <- ancinp
  bstimes <- seq(0,fsimtime,by=tstep)
  #bstimes[maxindex] <- bstimes[maxindex]
  times(bsim) <- bstimes
  #parms(bsim) <- c(X1=A1.init,X2=A2.init)
  bsim <- sim(bsim)
  bsout <- out(bsim)
  ancfn <- splinefun(bsout$time,bsout$A1+bsout$A2)
  ancfNh <- function(x){
    ancfn(x)-1
  }
  tmrca <- uniroot.all(ancfNh,lower=0,upper=max(bstimes))
  kfn <- splinefun(bsout$time,bsout$K)
  sack <- as.double(kfn(as.double(tmrca)))
  result[["parms"]] <- parms
  result[["fsout"]] <- fsout
  result[["bsout"]] <- bsout
  result[["A1init"]] <- A1.init
  result[["A2init"]] <- A2.init
  result[["C11"]] <- bsout$C11[maxindex-1]
  result[["C12"]] <- bsout$C12[maxindex-1]
  result[["C22"]] <- bsout$C22[maxindex-1]
  result[["C"]] <- bsout$C11[maxindex-1]+bsout$C12[maxindex-1]+bsout$C22[maxindex-1]
  result[["Cnorm"]] <- result[["C"]]/(result[["A1init"]]+result[["A2init"]])
  result[["r"]] <- calcclus(result[["C11"]],result[["C12"]],result[["C22"]])
  result[["tmrca"]] <- tmrca
  result[["K"]] <- sack
  result[["nlineages"]] <- fsout$Ic[maxindex-1] + fsout$Ih[maxindex-1]
  return(result)
}

r1 <- xRhc 
numr1 <- length(r1)

dr.col.prop.sim.list <- list()

# Simulating individually, from index cow, for simtime years, with weekly timestamps. 
simulate <- function(i,paras,simtime){
  dr.col.prop.par <- as.list(c(muC=0.1, muH=1/3, 
                               Bhh=paras[i,]$Bhh, Bcc=paras[i,]$Bcc, 
                               Bch=paras[i,]$Bch, Bhc=paras[i,]$Bhc))
  dr.col.init <- c(Sc=200,Ic=1,Sh=201,Ih=0)
  dr.col.sim <- coalsim(fsim=tbmodel.dr,
                        fsiminits=dr.col.init,
                        params=dr.col.prop.par,
                        fsimtime=simtime,
                        tstep=1/52,
                        matfun=tbmodel.dr.makemat,
                        bsim=anc,
                        frac1=NULL,
                        frac2=NULL,
                        numseq=NULL)
  return(dr.col.sim)
}
enableJIT(1)

# 300 years only
system.time(dr.col.simOUT <- simulate(1,params.setB,300))
dr.col.prop.sim.list[[1]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim1.rds")
system.time(dr.col.simOUT <- simulate(2,params.setB,300))
dr.col.prop.sim.list[[2]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim2.rds")
system.time(dr.col.simOUT <- simulate(3,params.setB,300))
dr.col.prop.sim.list[[3]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim3.rds")
system.time(dr.col.simOUT <- simulate(4,params.setB,300))
dr.col.prop.sim.list[[4]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim4.rds")
system.time(dr.col.simOUT <- simulate(5,params.setB,300))
dr.col.prop.sim.list[[5]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim5.rds")
system.time(dr.col.simOUT <- simulate(6,params.setB,300))
dr.col.prop.sim.list[[6]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim6.rds")
system.time(dr.col.simOUT <- simulate(7,params.setB,300))
dr.col.prop.sim.list[[7]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim7.rds")
system.time(dr.col.simOUT <- simulate(8,params.setB,300))
dr.col.prop.sim.list[[8]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim8.rds")
system.time(dr.col.simOUT <- simulate(9,params.setB,300))
dr.col.prop.sim.list[[9]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim9.rds")
system.time(dr.col.simOUT <- simulate(10,params.setB,300))
dr.col.prop.sim.list[[10]] <- dr.col.simOUT
saveRDS(dr.col.simOUT,"Rhc15_sim10.rds")

#######################################################
##  Visualizing tree metrics                         ##
#######################################################
dr.col.prop.sim.list <- list()
dr.col.prop.sim.list[[1]] <- readRDS("Rhc15_sim1.rds")
dr.col.prop.sim.list[[2]] <- readRDS("Rhc15_sim2.rds")
dr.col.prop.sim.list[[3]] <- readRDS("Rhc15_sim3.rds")
dr.col.prop.sim.list[[4]] <- readRDS("Rhc15_sim4.rds")
dr.col.prop.sim.list[[5]] <- readRDS("Rhc15_sim5.rds")
dr.col.prop.sim.list[[6]] <- readRDS("Rhc15_sim6.rds")
dr.col.prop.sim.list[[7]] <- readRDS("Rhc15_sim7.rds")
dr.col.prop.sim.list[[8]] <- readRDS("Rhc15_sim8.rds")
dr.col.prop.sim.list[[9]] <- readRDS("Rhc15_sim9.rds")
dr.col.prop.sim.list[[10]] <- readRDS("Rhc15_sim10.rds")

dr.col.prop.sim.Cnorm <- unlist(lapply(dr.col.prop.sim.list,"[","Cnorm"))
dr.col.prop.sim.r <- unlist(lapply(dr.col.prop.sim.list,"[","r"))
dr.col.prop.sim.K <- unlist(lapply(dr.col.prop.sim.list,"[","K"))
nseqs <- unlist(lapply(dr.col.prop.sim.list,"[","nlineages"))
dr.col.prop.sim.Knorm <-(dr.col.prop.sim.K-2*nseqs*log(nseqs))/(2*nseqs*log(nseqs))

# For visualizing set of simulations for given sampling time with different Rhc values 
# # pdf("newfig7.pdf",width=6,height=4)
# par(mfrow=c(1,3),pty="s",font.lab=2,las=1)
# plot(dr.col.prop.sim.Cnorm~r1,ylim=c(0.27,0.37),main="(a) Normalised cherries",
#      ylab=expression(bold(C[norm])),xlab=expression(paste(bold(R)[HC])),
#      col="black",type="n",lwd=2,xlim=c(0,66.3))
# lines(dr.col.prop.sim.Cnorm~r1,lwd=2)
# abline(h=0.3333333333333,col="black",lty=2)
# plot(dr.col.prop.sim.r~r1,ylim=c(-0.2,0.2),main="(b) Assortativity",
#      ylab="r",xlab=expression(paste(bold(R)[HC])),
#      col="black",type="n",lwd=2,xlim=c(0,66.3))
# abline(h=0,col="black",lty=2)
# lines(dr.col.prop.sim.r~r1,lwd=2)
# plot(dr.col.prop.sim.Knorm~r1,ylim=c(-0.1,0.2),main="(c) Normalised Sackin",
#      ylab=expression(bold(bar(K))),xlab=expression(paste(bold(R)[HC])),
#      col="black",type="n",lwd=2,xlim=c(0,66.3))
# lines(dr.col.prop.sim.Knorm~r1,lwd=2)
# abline(h=0,col="black",lty=2)
# dev.off()


# visualising metrics at different sampling times

RhcTsamples <- data.frame(tsample=rep(c("50","75","100","200","300"),each=10),
                          Rhc=rep(xRhc,5),
                          Cnorm=NA,
                          r=NA,
                          K=NA,
                          Knorm=NA)

# Load Rhc11*.rds to populate first 10 rows
# RhcTsamples$Cnorm[1:10] <- dr.col.prop.sim.Cnorm
# RhcTsamples$r[1:10] <- dr.col.prop.sim.r
# RhcTsamples$K[1:10] <- dr.col.prop.sim.K
# RhcTsamples$Knorm[1:10] <- dr.col.prop.sim.Knorm
# rm(dr.col.prop.sim.list,dr.col.prop.sim.Cnorm,dr.col.prop.sim.r,dr.col.prop.sim.Knorm,dr.col.prop.sim.K)

# Do the same for Rhc12*.rds
# RhcTsamples$Cnorm[11:20] <- dr.col.prop.sim.Cnorm
# RhcTsamples$r[11:20] <- dr.col.prop.sim.r
# RhcTsamples$K[11:20] <- dr.col.prop.sim.K
# RhcTsamples$Knorm[11:20] <- dr.col.prop.sim.Knorm
# rm(dr.col.prop.sim.list,dr.col.prop.sim.Cnorm,dr.col.prop.sim.r,dr.col.prop.sim.Knorm,dr.col.prop.sim.K)

# Do the same for Rhc13*.rds
# RhcTsamples$Cnorm[21:30] <- dr.col.prop.sim.Cnorm
# RhcTsamples$r[21:30] <- dr.col.prop.sim.r
# RhcTsamples$K[21:30] <- dr.col.prop.sim.K
# RhcTsamples$Knorm[21:30] <- dr.col.prop.sim.Knorm
# rm(dr.col.prop.sim.list,dr.col.prop.sim.Cnorm,dr.col.prop.sim.r,dr.col.prop.sim.Knorm,dr.col.prop.sim.K)

# Do the same for Rhc14*.rds
# RhcTsamples$Cnorm[31:40] <- dr.col.prop.sim.Cnorm
# RhcTsamples$r[31:40] <- dr.col.prop.sim.r
# RhcTsamples$K[31:40] <- dr.col.prop.sim.K
# RhcTsamples$Knorm[31:40] <- dr.col.prop.sim.Knorm
# rm(dr.col.prop.sim.list,dr.col.prop.sim.Cnorm,dr.col.prop.sim.r,dr.col.prop.sim.Knorm,dr.col.prop.sim.K)

# Do the same for Rhc15*.rds
# RhcTsamples$Cnorm[41:50] <- dr.col.prop.sim.Cnorm
# RhcTsamples$r[41:50] <- dr.col.prop.sim.r
# RhcTsamples$K[41:50] <- dr.col.prop.sim.K
# RhcTsamples$Knorm[41:50] <- dr.col.prop.sim.Knorm

RhcTsamples$time <- rep(c(50,75,100,200,300),each=10)
# ggplot. seeing metrics change with Rhc and sampling time.
library(ggplot2)
library(gridExtra)

qC <- ggplot(data=RhcTsamples,aes(x=Rhc,y=Cnorm)) +
  geom_line(aes(colour=factor(time))) +
  geom_hline(yintercept=0.33333333333) +
  theme_bw()+
  theme(legend.position="none")
qr <- ggplot(data=RhcTsamples,aes(x=Rhc,y=r)) +
  geom_line(aes(colour=factor(time))) +
  geom_hline(yintercept=0) +
  theme_bw()+
  theme(legend.position="none")
qK <- ggplot(data=RhcTsamples,aes(x=Rhc,y=Knorm)) +
  geom_line(aes(colour=factor(time))) +
  geom_hline(yintercept=0) +
  theme_bw()+
  theme(legend.position="none")
grid.arrange(qC,qr,qK,ncol=3,nrow=1)

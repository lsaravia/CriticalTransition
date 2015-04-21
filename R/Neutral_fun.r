
# Generate inp parameters file for simulation of neutral model
#
# fname: file name
# side: side of simulation square grid 
# mprob: metacomunity vector of species frecuencies
# birth: birth rate
# mortality: 
# dispersal: parameter of the dispersal kernel (power/exponential)
# colonization: migration rate from metacommunity
# replace: replacement rate
#
genNeutralParms <- function(fname,side,mprob,birth,mortality,dispersal,colonization,replace=0){
  S <- length(mprob)
  parm <- data.frame( n=c(rep(0,3),1:S) ,prob=c(rep(0,3),prob=mprob))
  parm$text <- "" 
  parm$z <- 0
  parm$text[1] <- paste(side,side,sep="\t")
  parm$text[2] <- S
  parm$text[3] <- paste(0,birth,mortality,dispersal,colonization,replace,sep="\t")
  nff <- S+3
  parm$text[4:nff] <-   with(parm[4:nff,], paste(n,z,z,z,format(sort(prob),scientific=F),sep="\t"))                    
  if(!grepl(".inp",fname)) fname <-paste0(fname,".inp")

  write.table(parm$text,fname,sep="\t",row.names=F,col.names=F,quote=F)
}


# Generates pomac.lin paramenter file for multiple simulations
#
genPomacParms <- function(fname,GrowthR,MortR,DispD,ColonR,ReplaceR,numRep=1)
  {
  pom <- expand.grid(GrowthRate=GrowthR,MortalityRate=MortR, DispersalDistance=DispD, 
    ColonizationRate=ColonR,ReplacementRate=ReplaceR)
  pom <- pom[rep(seq_len(nrow(pom)), numRep), ]
  if(!grepl(".lin",fname)) fname <-paste0(fname,".lin")
  write.table(pom,fname,sep="\t",row.names=F,col.names=T,quote=F)
}


# Calculate Ranks for every parameter combination, for ggplot2  
#
calcRankSAD  <- function(den)
{
  require(plyr)
  hh <- function(x) { 
  x1 <- x[x$value>0,]
  x1$parms <- paste(unique(x[,1:4]),collapse="_")
  x1$Rank <- nrow(x1) - rank(x1$value) +1
  return(x1)
  }
  ddply(den, .(MortalityRate,DispersalDistance,ColonizationRate,ReplacementRate),hh )
}

# Calculate Ranks for every parameter combination in columns cols, for variable vv, using data.frame den 
# to use ggplot2  
#
calcRankSAD_by  <- function(den,vv,cols)
{
  require(plyr)
  hh <- function(x) { 
  x1 <- x[x[[vv]]>0,]
  x1$parms <- paste(unique(x[,cols]),collapse="_")
  x1$Rank <- nrow(x1) - rank(x1[[vv]],ties.method="random") +1
  return(x1)
  }
  ddply(den, cols,hh )
}



# pairwise Anderson-Darling K test for all combinations of parameters in variable parms 
# vv: name of the variable where density or proportion is
# denl: data.frame with all data
#
pairwiseAD_SAD <- function(denl,vv,parms){
  require(kSamples)

  parms <- unique(parms)
  combo <- combn(nrow(parms),2)
  
  require(plyr)
  mks <-adply(combo,2, function(x) {
    p1 <- parms[x[1],]
    p2 <- parms[x[2],]
    d1 <- merge(denl,p1)
    d2 <- merge(denl,p2)
    ks <- ad.test(list(d1[,vv],d2[,vv]),method="simulated",nsim=1000)
    out <-data.frame(p1,p2,stat=ks$ad[2,2],p.value=ks$ad[2,4],stringsAsFactors=F)
    ln <-length(names(p1))*2
    names(out)[1:(ln)]<-c(paste0(abbreviate(names(p1)),1),paste0(abbreviate(names(p1)),2))
    return(out)      
  })
  mks$p.adjust <- p.adjust(mks$p.value, method="hommel")
  return(mks)
}


# pairwise Anderson-Darling K test for all combinations of 
# parameters in variable parms with repetitions, BEWARE last variable in parms must be 
# the repetition
# 
# vv: name of the variable where density or proportion is
# denl: data.frame with all data
#
pairwiseAD_Dif <- function(denl,vv,parms){
  require(kSamples)
  parms <- unique(parms)
  combo <- combn(nrow(parms),2)
  nc <- ncol(parms)-2
  #pb <- txtProgressBar(min = 0, max = ncol(combo), style = 3)
  #i <- 0
  require(plyr)
  mks <-adply(combo,2, function(x) {
    
    p1 <- parms[x[1],]
    p2 <- parms[x[2],]
    out<-NULL
    #i <<- i+1 
    #setTxtProgressBar(pb, i)
    # test the error!!!!!!!!!
    if(sum(p1[,1:nc]==p2[,1:nc])==nc) {
      d1 <- merge(denl,p1)
      d2 <- merge(denl,p2)
      if(nrow(d1)>2 & nrow(d2)>2) {
        ks <- ad.test(list(d1[,vv],d2[,vv]),method="simulated",Nsim=1000)
        out <-data.frame(p1,p2,stat=ks$ad[2,2],p.value=ks$ad[2,4],stringsAsFactors=F)
        ln <-length(names(p1))*2
        names(out)[1:(ln)]<-c(paste0(abbreviate(names(p1)),1),paste0(abbreviate(names(p1)),2))
      }
    }
    return(out)      
  })
  #close(pb)
  mks$p.adjust <- p.adjust(mks$p.value, method="hommel")
  return(mks)
}





# Read simulation output and change from wide to long format NO TIME
#
meltDensityOut_NT <- function(fname,num_sp){
  if(!grepl("Density.txt",fname)) fname <- paste0(fname,"Density.txt")

  den <- read.delim(fname,stringsAsFactors = F)
  names(den)[1:5]<-c("GrowthRate","MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate")
  
  # from 7 to 473 there are species densities
  # Put simpler names to variables to identify species
  names(den)[7:(6+num_sp)]<-as.character(1:num_sp)

  require(plyr)
  
  den <- ddply(den, 1:5, function(x){ i <- c(1:nrow(x)); data.frame(x,rep=i)})
  
  require(reshape2)
  
  den1 <- melt(den,id.vars=c("MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate","rep"),measure.vars=c(7:(6+num_sp)),variable.name="Species")
  den1 <- den1[den1$value!=0,] 
}

# Read simulation output and set variable names in wide format 
# Adds a variable Rep if there are several repeated simulations
#
readWideDensityOut <- function(fname){
  if(!grepl("Density.txt",fname)) fname <- paste0(fname,"Density.txt")
  den <- read.delim(fname,stringsAsFactors = F)
  names(den)[1:5]<-c("GrowthRate","MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate")
  names(den)[7] <- unlist(strsplit(names(den)[7],".",fixed=T))[1]
  if(nrow(den)<3)
    eTime <-1
  else
    eTime <- (max(den$Time)/(den$Time[3]-den$Time[2]))+1
  
  if(  eTime < nrow(den) ){
    den$Rep <- rep( 1:(nrow(den)/eTime),each=eTime)
  }
  return(den)
}

# Read simulation output and set variable names in wide format 
# Asumme that are all repetitions of the same set of parameters and 
#  Adds a variable Rep 
#
readClusterOut <- function(fname){
  if(!grepl("CluSizes.txt",fname)) fname <- paste0(fname,"CluSizes.txt")
  den <- read.delim(fname,stringsAsFactors = F)
  names(den)[1:6]<-c("GrowthRate","MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate","Time")
  #  
  #den$Rep <- c(1:nrow(den))
  #
  return(den)
}


# Calculates Dq from a data.frame read from the output of neutral model 
# auxiliar function for the next one
#
calcDq_frame <- function(pp)
{
  pp$Dq  <- with(pp,ifelse(q==1,alfa,Tau/(q-1)))
  pp$SD.Dq  <- with(pp,ifelse(q==1,SD.alfa,abs(SD.Tau/(q-1))))
  pp$R.Dq <- with(pp,ifelse(q==1,R.alfa,R.Tau))
  nc <- ncol(pp)
  return(pp[,c(1:6,(nc-2):nc)])
}              
# Reads the output of multifractal spectra of neutral model
# an calculates Dq 
#
readNeutral_calcDq <-function(fname)
{
  md1 <- read.table(fname,header=F,skip=1)
  md1 <- md1[,c(2:16)]
  names(md1)<-c("MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate","Time","q","Tau","alfa","f(alfa)","R.Tau","R.alfa","R.f","SD.Tau","SD.alfa","SD.f")
  
  md1 <-calcDq_frame(md1)
}


# Read Dq output from neutral model and calculate pairwise differences 
#
# Dqf: data frame from readNeutral_calcDq 
# qNumber: number of q used
#
compDq_frame <- function(Dqf,qNumber)
{
  if( !require(statmod) & !require(reshape2))
    stop("required statmod and reshape2")
  
  # Subset for testing 
  #
  #Dqf <- with(Dqf,Dqf[MortalityRate==.2 & DispersalDistance==0.04 & ColonizationRate==0.001, ])
  
  # Set the number of repetitions using nrow and number of q  
  # 
  Dqf$rep <- rep( 1:(nrow(Dqf)/qNumber),each=qNumber)
  
  # Build variable for comparisons
  Dqf$factor <- do.call(paste, c(Dqf[,1:4],sep="_"))
    
  # Prepare data.frame in wide format 
  #
  Dq2 <- melt(Dqf, id.vars=c("q","rep","factor"),measure.var="Dq")
  Dq2 <- dcast(Dq2, factor+rep~ q)
  
  # Compare SRS curves
  #
  c2 <- compareGrowthCurves(Dq2$factor,Dq2[,3:37],nsim=1000)
}

# Plot Dq with fixed parameters except ReplacementRate
# Calculate SD from repeated simulations
#
plotDq_ReplaceR <- function(Dqf,MortR,DispD,ColonR,tit="")
{
  require(plyr)
  c3 <- with(Dqf,Dqf[MortalityRate==MortR & DispersalDistance==DispD & ColonizationRate==ColonR, ])
  c3$factor <- do.call(paste, c(c3[,1:4],sep="_"))
  c3 <- ddply(c3, .(factor,q), summarize, SD.Dq=sd(Dq),Dq=mean(Dq))

  require(ggplot2)
  gp <- ggplot(c3, aes(x=q, y=Dq, colour=factor)) +
    geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1) +
    geom_point() + theme_bw() + ggtitle(tit)
  print(gp)

}

# Plot Dq by factor "grp" shows SD from data.frame
#
plotDq <- function(Dq1,grp) 
  {
  require(ggplot2)
  print(gp <- ggplot(Dq1, aes_string(x="q", y="Dq", colour=grp)) +
          geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1) +
          geom_point() + theme_bw()) 
  }





# Read a sed file in a matrix
#
# fname: file name of the sed file
#
read_sed <- function(fname)
{
  d <-read.table(fname, nrows=1,header=F)
  per <-data.matrix(read.table(fname, skip=2,header=F))
  if(d$V2!=nrow(per)) stop(paste("Incorrect formated sed file:",fname))
  return(per)
}

read_sed2xy <- function(fname)
{
  spa <- read_sed(fname)
  z <- 1:(nrow(spa)*ncol(spa))
  zpa <- data.frame(v=spa[z],x=trunc(z/ncol(spa)),y=1:nrow(spa))
}



# Read information of the fit of Dq (Zq) from t.file and q file
# Return a data.frame in long format
#
readZq <- function(fname,qname)
{
  zq <- read.table(fname, sep="\t",header=T)
  cna <- read_sed(qname)
  q <-t(cna)
  zq0 <- reshape(zq, timevar="q",times=q,v.names=c("logTr"),
                 varying=list(3:length(names(zq))),
                 direction="long")
}

# Function to plot Dq fit from t* files generated by mfSBA using ggplot2
# 
# zq0: data.frame with Zq, logTr 
# fac: factor to separate the lines 
#
plotDqFitG <- function(zq0,fac=3)
{
  require(ggplot2)
  require(dplyr)
  
  zq1 <- subset(zq0, q==1 | q==2 | q==3 | q==4 | q==5 | q==0 | q==-1 | q==-2 | q==-3 | q==-4 | q==-5 )
  zq1$logTr <- zq1$logTr+zq1$q/fac
  zq1 <- mutate(zq1,DqType=ifelse(grepl("SRS",Type),"DqSRS","DqSAD"), Type=ifelse(grepl("rnz",Type),"b) Randomized","a) Regular"))
  
#  g <- ggplot(zq1,aes(LogBox,logTr,colour=factor(q))) + geom_point(aes(shape=factor(q))) + 
#    scale_color_discrete(name="q") + 
#    
  
#  g <- ggplot(zq1,aes(LogBox,logTr,shape=factor(q))) + geom_point(aes(shape=factor(q)),size=1) + 
  g <- ggplot(zq1,aes(LogBox,logTr,colour=factor(q))) + geom_point(aes(shape=factor(q)),size=1) + 
      geom_smooth(method="lm",se=F)  
  g <- g + scale_shape_manual(values=c(0,1,2,3,4,5,6,8,15,16,17,21:24),name="q") 
  g <- g + scale_colour_brewer(palette="Set1",name="q")
  g <- g + ylab(expression(italic(paste("log ",  Z[q](epsilon) )))) + theme_bw() +
      xlab(expression(italic(paste("log ",epsilon))))+
    facet_wrap(Type ~ DqType, scales="free")
}


plotDqFitGQ <- function(zq0,qq,fac=3)
{
  require(ggplot2)
  require(dplyr)
  
  zq1 <- filter(zq0, q %in% qq  )
  zq1$logTr <- zq1$logTr+zq1$q/fac
  zq1 <- mutate(zq1,DqType=ifelse(grepl("SRS",Type),"SRS","SAD"), Type=ifelse(grepl("rnz",Type),"b) Randomized","a) Regular"))
  
  #  g <- ggplot(zq1,aes(LogBox,logTr,colour=factor(q))) + geom_point(aes(shape=factor(q))) + 
  #    scale_color_discrete(name="q") + 
  
  g <- ggplot(zq1,aes(LogBox,logTr,shape=factor(q))) + geom_point(aes(shape=factor(q)),size=1) + 
    geom_smooth(method="lm",se=F,colour="grey")  +
    scale_shape_manual(values=c(0,1,2,3,4,5,6,8,15,16,17,21:24),name="q") +
    ylab(expression(italic(paste("log ",  Z[q](epsilon) )))) + theme_bw() +
    xlab(expression(italic(paste("log ",epsilon))))+
    facet_wrap(Type ~ DqType, scales="free")
}

# Calculates theoretic Dq from pmodel
#
calcDqTeor <- function(q,p1) {
  if( q==1)
    q <- q+1e-10
  p2 <- p3 <- p4 <- (1-p1)/3
  f1 <- p1/(p1+p2+p3+p4)
  f2 <- p2/(p1+p2+p3+p4)
  f3 <- p3/(p1+p2+p3+p4)
  f4 <- p4/(p1+p2+p3+p4)
  dq <- log2(f1^q+f2^q+f3^q+f4^q)/(1-q)
  return(dq)
  }


# Save a matrix as a sed file with type BI (floating point)
#
save_matrix_as_sed <- function(mat,fname)
{
  header <- paste(nrow(mat),ncol(mat),"\nBI")
  write.table(header,file=fname,row.names=F,col.names=F,quote=F)
  write.table(mat,file=fname,row.names=F,col.names=F,quote=F,append=T)
}

calcDq_multiSBA <- function(fname,parms,pathBin="",recalc=FALSE)
{
  sname <- paste0("s.", fname)
  if((!file.exists(sname)) | recalc)
  {
    if(nchar(pathBin)==0)
    {
      syst.txt <- paste("./multiSpeciesSBA ",fname, parms)
    } else {
      syst.txt <- paste0(pathBin,"/multiSpeciesSBA ",fname," ",parms)
    }

    system(syst.txt)
  }
  pp <- read.delim(sname, header=T)
  
  for(nc in 1:ncol(pp)){
    if( class(pp[,nc ])=="factor") pp[,nc]<-as.numeric(as.character(pp[,nc])) 
  }

  pp$Dq  <- with(pp,ifelse(q==1,alfa,Tau/(q-1)))
  pp$SD.Dq  <- with(pp,ifelse(q==1,SD.alfa,abs(SD.Tau/(q-1))))
  pp$R.Dq <- with(pp,ifelse(q==1,R.alfa,R.Tau))
    
  return(pp[,c("q","Dq","SD.Dq","R.Dq")])
}


# Generate a banded image with nSp species and side = side
#
genUniformSAD_image <- function(nSp,side,rnd=F)
{
  if( side %% nSp != 0 )
    stop("Number of species [nSp] must divide [side]")
  if(rnd) {
    repeat {
      v <- rpois(nSp,side*side/nSp)
      if(side*side>sum(v[1:nSp-1])) break
    }
    v[nSp] <- side*side-sum(v[1:nSp-1])
    vv <- rep(1:nSp,times=v)
    matrix(vv,nrow=side)
    } else matrix(rep(1:nSp,each=side*side/nSp),nrow=side)
}

# Generate a regular image with Fisherian SAD with nsp species and side = side
# Requires untb package
# Returns a matrix representing the spatial distribution and a vector with proportions of each sp
#
genFisherSAD_image <- function(nsp,side)
{
  require(untb)
  N <- side*side
  repeat {
    ff<-fisher.ecosystem(N=N,S=nsp,nmax=N)
    if(nsp==nrow(ff)) {
      if(sum(ff)!=N) {
        c<- N/sum(ff)
        ff <- ceiling(ff*c)
      }
      m <- matrix(rep(1:nsp,ff),nrow=side)
      if(m[(sum(ff)+1)]==1) m[(sum(ff)+1):length(m)]<-0
      if(ncol(m)>=side) break
    }
  }
 
  prob <- ff/sum(ff)
  return(list("m"=m,"prob"=prob))
}


# Generate Logseries SAD with nsp*f species and side = side
#  The factor f is to compensate losses for use as a metacommunity
# Requires untb package
#
#
genFisherSAD <- function(nsp,side,f=1.33)
{
  require(untb)
  N <- side*side
#  alpha <-  fishers.alpha(N, nsp)
#  x <- N/(N + alpha)
#  j <- 1:(nsp)
#  prob <-  alpha * x^j/j  
  # Normalize
#  prob <- prob/sum(prob)  
  nsp <- ceiling(nsp*f)
  repeat {
    ff<-fisher.ecosystem(N=N,S=nsp,nmax=N)
    if(nsp==nrow(ff)) {
      if(sum(ff)!=N) {
        c<- N/sum(ff)
        ff <- ceiling(ff*c)
      }
      break
    }
  }
  prob <- ff/sum(ff)
  return(prob)
}

  

# Generate a Neutral model with Fisher logseries metacommunity SAD 
# randomize it and calculates SRS and DqSAD multifractal estimations
# Neutral model Hierarchical saturated
# Graph : plot graphs
# meta: metacommunity "L" logseries, any other uniform
#
compMethods_NeutralSAD <- function(nsp,side,simul=T,graph=T,meta="L") {
  if(!exists("mfBin")) stop("Variable mfBin not set (mfSBA binary)")
  if(!exists("neuBin")) stop("Variable neuBin not set (neutral binary)")
  if(!require(untb))  stop("Untb package not installed")

  if(simul){

    #N <- side*side   # The metacommunity have 100 times more individuals
    #alpha <-  fishers.alpha(N, nsp)
    #x <- N/(N + alpha)
    #j <- 1:(nsp)
    #prob <-  alpha * x^j/j  
    # Normalize
    #prob <- prob/sum(prob)
    if(toupper(meta)=="L") {
      prob <- genFisherSAD(nsp,side)
      neuParm <- "fishE"
      bname <- paste0("neuFish",nsp)
      sadName <- "Neutral"
    } else {
      prob <- rep(1/nsp,nsp)  
      neuParm <- "unifE"
      bname <- paste0("neuUnif",nsp)
      sadName <- "NeuUnif"
    }

    genNeutralParms(neuParm,side,prob,1,0.2,0.4,0.0001)


    # Delete old simulations
    system(paste0("rm ",bname,"*.txt"))

    par <- read.table("sim.par",quote="",stringsAsFactors=F)
    # Change base name

    par[par$V1=="nEvals",]$V2 <- 500
    par[par$V1=="inter",]$V2 <- 500 # interval to measure Density and Diversity
    par[par$V1=="init",]$V2 <- 500  # Firs time of measurement = interval
    par[par$V1=="modType",]$V2 <- 4 # Hierarchical saturated
    par[par$V1=="sa",]$V2 <- "S" # Save a snapshot of the model
    par[par$V1=="baseName",]$V2 <- bname# Time = 100 
    par[par$V1=="minBox",]$V2 <- 2
    par[par$V1=="pomac",]$V2 <- 0 # 0:one set of parms 
                                  # 1:several simulations with pomac.lin parameters 

    write.table(par, "sim.par",sep="\t",row.names=F,col.names=F,quote=F)

    system(paste(neuBin,"sim.par",paste0(neuParm,".inp")))
  }
  #fname <- paste0("neuFish",nsp,"Density.txt")
  #sad1 <- meltDensityOut_NT(fname,nsp)

  fname <- paste0(bname,"-0500.sed")
  spa <- read_sed(fname)

  sad1 <- data.frame(table(spa),Type="SAD",Side=side,NumSp=nsp,SAD=sadName)
  #sad1$alpha <- fishers.alpha(side*side, nrow(sad1))

  if(graph) plot_sed_image(spa,paste("Neutral T500",nsp),0,nsp,0)

  Dq1<- calcDq_multiSBA(fname,"q.sed 2 1024 20 S",mfBin,T)
  Dq1$Type <- "SRS"

  # Randomize the spatial distribution
  #
  spa <- matrix(sample(spa),nrow=side)
  fname1 <- paste0(bname,"-500rnz.sed")

  save_matrix_as_sed(spa,fname1)
  if(graph) plot_sed_image(spa,paste("Neutral T500 Rnz",nsp),0,nsp,0)

  Dq2<- calcDq_multiSBA(fname1,"q.sed 2 1024 20 S",mfBin,T)
  Dq2$Type <- "rnzSRS"
  Dq1<- rbind(Dq1,Dq2)
  if(graph) {
    plotDq(Dq1,"Type")

    bin <- range(Dq1$R.Dq)
    bin <- (bin[2]-bin[1])/10
    print(ggplot(Dq1, aes(x=R.Dq,fill=Type)) + geom_histogram(alpha=0.2,binwidth = bin))
  }
  
  # Now calculate DqSAD

  Dq3<- calcDq_multiSBA(fname,"q.sed 2 1024 20 E",mfBin,T)
  Dq3$Type <- "DqSAD"
  #Dq3<- rbind(Dq3,Dq2)

  Dq2<- calcDq_multiSBA(fname1,"q.sed 2 1024 20 E",mfBin,T)
  Dq2$Type <- "rnzDqSAD"
  Dq3<- rbind(Dq3,Dq2)

  if(graph) {
    plotDq(Dq3,"Type")

    plotDqFit(paste0("t.", fname1),"q.sed")

    bin <- range(Dq3$R.Dq)
    bin <- (bin[2]-bin[1])/10
    print(ggplot(Dq3, aes(x=R.Dq,fill=Type)) + geom_histogram(alpha=0.2,binwidth = bin))
  }  
  #require(pander)
  #pandoc.table(Dq3[Dq3$R.Dq<0.6,],caption="R2<0.6")
  
  Dqt <- rbind(Dq1,Dq3)
  Dqt$Side <-side
  Dqt$NumSp <-nsp
  Dqt$SAD <- sadName #Neutral with uniform metacommunity

  return(list("Dq"=Dqt,"SAD"=sad1))
}

# Calculates fractal infomation dimension and species area relationship
#
#
calc_CommunityDq <- function(nsp,side,disp,migr,repl,sims=10,meta="U"){
  require(plyr)
  require(dplyr)

  Dq <- ldply(repl,function(r){
    if(toupper(meta)=="L") {
      bname <- paste0("neuFish",nsp,"_",side,"R", r)
    } else {
      bname <- paste0("neuUnif",nsp,"_",side,"R", r)
    }
    fname <- paste0(bname,"mfOrd.txt")
    Dq1 <- readNeutral_calcDq(fname)
    Dq1$DqType <- "DqSRS"
    fname <- paste0(bname,"mfSAD.txt")
    Dq2 <- readNeutral_calcDq(fname)
    Dq2$DqType <- "DqSAD"
  
    return(rbind(Dq1,Dq2))
  
  })
  time <- unique(Dq$Time)
  iT <- max(time) - (time[3]-time[2])*10
  mDq1 <- Dq %>% group_by(ReplacementRate) %>% filter(Time>=iT, DqType=="DqSRS", q==1 ) %>% mutate(Rep=ntile(ReplacementRate,sims)) %>%
    group_by(MortalityRate,DispersalDistance,ColonizationRate,ReplacementRate,Rep ) %>% summarise(D1=mean(Dq))
  
  mDq0 <- Dq %>% group_by(ReplacementRate) %>% filter(Time>=iT, DqType=="DqSAD", q==0 ) %>% mutate(Rep=ntile(ReplacementRate,sims))  %>%
    group_by(MortalityRate,DispersalDistance,ColonizationRate,ReplacementRate,Rep ) %>% summarise(D0=mean(Dq))
  
  mDq <- inner_join(mDq1,mDq0,by=c("MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate","Rep")) %>% 
    mutate( Nsp=nsp, Side=side)
  
  return(mDq)
}

# Calculates Bray-Curtis dissimilarity
#             Kullback - Leiber divergence
#             Kolmogorov-smirnov distance
#
calc_CommunityDist <- function(nsp,side,disp,migr,repl,sims=10,meta="U"){
  require(plyr)
  require(dplyr)
  if(repl[1]==repl[2]) repl <- repl[1]
  den <- ldply(repl,function(r){
    if(toupper(meta)=="L") {
      bname <- paste0("neuFish",nsp,"_",side,"R", r)
    } else {
      bname <- paste0("neuUnif",nsp,"_",side,"R", r)
    }
    readWideDensityOut(bname)
    
  })

# Toma los ultimas 10 simulaciones para calcular
#
iT <- max(den$Time) - (den$Time[3]-den$Time[2])*10
d1 <- den %>% group_by(ReplacementRate,Rep) %>% filter(Time>=iT)  %>% select(starts_with("X")) %>% 
              summarise_each(funs(mean)) %>%
              ungroup()

#d1 <- den %>% group_by(ReplacementRate,Rep) %>% filter(Time==max(Time))  %>% select(starts_with("X")) %>% ungroup()
nc <-d1  %>% group_by(ReplacementRate) %>% summarise(n=n())
re <- nc$n[1]

parms <- unique(d1[,1:2])
combo <- combn(nrow(parms),2)
if(length(repl)>1){
  combo <- combo[,combo[2,]>re]
  combo <- combo[,combo[1,]<=re]
}
nc <- ncol(d1)

mks <-adply(combo,2, function(x) {
          p1 <- parms[x[1],]
          p2 <- parms[x[2],]
#          s1 <- merge(d1,p1)
#          s2 <- merge(d1,p2)
          s1 <- d1[x[1],]
          s2 <- d1[x[2],]
          k <- ks.test(as.numeric(s1[,3:nc]),as.numeric(s2[,3:nc]))
          out <-data.frame(p1,p2,bray=brayCurtis(s1[,3:nc],s2[,3:nc]),
                           KL=KLDiv(s1[,3:nc],s2[,3:nc]),
                           ks=k$statistic,ks.p=k$p.value,stringsAsFactors=F)
        })
}

# Calculates Bray-Curtis dissimilarity
#             Kullback - Leiber divergence
#             Kolmogorov-smirnov distance
#
calc_CommunityDist_1T <- function(nsp,side,disp,migr,repl,time,sims=10,meta="U"){
  require(plyr)
  require(dplyr)
  if(repl[1]==repl[2]) repl <- repl[1]
  den <- ldply(repl,function(r){
    if(toupper(meta)=="L") {
      bname <- paste0("neuFish",nsp,"_",side,"R", r, "T", time)
    } else {
      bname <- paste0("neuUnif",nsp,"_",side,"R", r, "T", time)
    }
    readWideDensityOut(bname)
    
  })

# Hace promedios de las densidades de las repeticiones para calcular
#
d1 <- den %>% group_by(ReplacementRate) %>% select(starts_with("X")) %>% 
              summarise_each(funs(mean)) %>%
              ungroup()

#d1 <- den %>% group_by(ReplacementRate,Rep) %>% filter(Time==max(Time))  %>% select(starts_with("X")) %>% ungroup()
nc <- d1 %>% group_by(ReplacementRate) %>% summarise(n=n())
re <- nc$n[1]
#re <-1
# 
# if(re>1) {
#   parms <- unique(d1[,1:2])
#   combo <- combn(nrow(parms),2)
#   if(length(repl)>1){
#     combo <- combo[,combo[2,]>re]
#     combo <- combo[,combo[1,]<=re]
#   }
#   nc <- ncol(d1)
#   ni <- 3                     # 2 first columns are parameters replacementRate Rep
# } else {
  parms <- unique(d1[,1])
  combo <- combn(nrow(parms),2)
  if(length(repl)>1){
    combo <- combo[,combo[2,]>re]
    combo <- combo[,combo[1,]<=re]
  }
  nc <- ncol(d1)
  ni <- 2                     # 1 first column is parameter   
#}

mks <-adply(combo,2, function(x) {
          p1 <- parms[x[1],]
          p2 <- parms[x[2],]
#          s1 <- merge(d1,p1)
#          s2 <- merge(d1,p2)
          s1 <- d1[x[1],]
          s2 <- d1[x[2],]
          k <- ks.test(as.numeric(s1[,ni:nc]),as.numeric(s2[,ni:nc]))
          out <-data.frame(p1,p2,bray=brayCurtis(s1[,ni:nc],s2[,ni:nc]),
                           KL=KLDiv(s1[,ni:nc],s2[,ni:nc]),
                           ks=k$statistic,ks.p=k$p.value,stringsAsFactors=F)
        })
}


# Kullback - Leiber divergence
#
KLDiv <- function(freqs1, freqs2)
{
  freqs1 = freqs1/sum(freqs1)
  freqs2 = freqs2/sum(freqs2)
  freqs1 <- freqs1[freqs2>0]
  freqs2 <- freqs2[freqs2>0]
  
  if (any(!(freqs2 > 0))) 
    warning("Vanishing value(s) in argument freqs2!")
  LR = ifelse(freqs1>0,log(freqs1/freqs2),0)
  
  KL = sum(freqs1 * LR)
  return(KL)
}

# Bray-Curtis dissimilarity
#
brayCurtis <- function (x, y) 
{
  bray <- 1 - sum(abs(x - y))/sum(x + y)
  return(bray)
}


# Simulate a time series of the neutral/hierarchical model
#
# disp: dispersal parameter
# migr: migration rate
# repl: replacement rate
# mortality fixed to 0.2
#
# simul: T=make simulations F=make plots
# sims:  Number of repetitions
# mf: calculate multifractal spectrum
# meta: U= uniform metacommunity
#       L= logseries metacommunity
#
simul_NeutralPlotTime <- function(nsp,side,disp,migr,repl,simul=T,time=1000,sims=10,mf="N",meta="U") {
  if(!exists("neuBin")) stop("Variable neuBin not set (neutral binary)")

  if(toupper(meta)=="L") {
    prob <- genFisherSAD(nsp,side)
    neuParm <- paste0("fishP",nsp,"_",side,"R", repl)
    bname <- paste0("neuFish",nsp,"_",side,"R", repl)
  } else {
    prob <- rep(1/nsp,nsp)  
    neuParm <- paste0("unifP",nsp,"_",side,"R", repl)
    bname <- paste0("neuUnif",nsp,"_",side,"R", repl)
  }
  pname <- paste0("pomacR",repl,".lin")

  if(simul){

    genNeutralParms(neuParm,side,prob,1,0.2,disp,migr,repl)

    # Delete old simulations
    system(paste0("rm ",bname,"m*.txt")) # MultiFractal mf
    system(paste0("rm ",bname,"D*.txt")) # Density
    system(paste0("rm ",bname,"C*.txt")) # Clusters


    par <- read.table("sim.par",quote="",stringsAsFactors=F)

    par[par$V1=="nEvals",]$V2 <- time
    par[par$V1=="inter",]$V2 <- 10 # interval to measure Density and Diversity
    par[par$V1=="init",]$V2 <- 1  # Firs time of measurement = interval
    par[par$V1=="modType",]$V2 <- 4 # Hierarchical saturated
    par[par$V1=="sa",]$V2 <- "N" # Save a snapshot of the model
    par[par$V1=="baseName",]$V2 <- bname# Time = 100 
    par[par$V1=="mfDim",]$V2 <- mf
    par[par$V1=="minBox",]$V2 <- 2
    par[par$V1=="pomac",]$V2 <- 1 # 0:one set of parms 
                                  # 1:several simulations with pomac.lin parameters 
    par[par$V1=="pomacFile",]$V2 <- pname # 0:one set of parms 
    par[par$V1=="minProp",]$V2 <- 0
    
    parfname <- paste0("sim",nsp,"_",side,"R", repl,".par")
    write.table(par, parfname, sep="\t",row.names=F,col.names=F,quote=F)

    genPomacParms(pname,1,c(0.2),disp,migr,repl,sims)
  
    # copy pomExp.lin to pomac.lin
    #system("cp pomExp.lin pomac.lin")
    s <- system("uname -a",intern=T)
    if(grepl("i686",s)) {
      system(paste(neuBin,parfname,paste0(neuParm,".inp")))
    } else {
      system(paste(neuBin64,parfname,paste0(neuParm,".inp")))
    }
  }
  den <-readWideDensityOut(bname)
  
  require(plyr)
  require(dplyr)
  
  if(!simul) {
    require(ggplot2)
    if(sims>10)  
      den <- filter(den,Rep %in%  sample(1:sims,10)) 

    print(ggplot(den, aes(x=Time, y=H,color=factor(Rep))) +
        geom_line() + theme_bw() +  ggtitle(paste(side,repl)))

#     print(ggplot(den, aes(x=Time, y=Richness,color=factor(Rep))) +
#         geom_line() + theme_bw() + ggtitle(paste(side,repl))) 
  }
 

  th <-function(x,tt,...){
    tt[which.max(x)]
  }
  
  d1 <- den %>% group_by(Rep) %>% summarize(MaxH=max(H),TMaxH=th(H,Time),
                  MaxRich=max(Richness),TMaxRich=th(Richness,Time))
  
  lInt <- max(den$Time)/(den$Time[3]-den$Time[2])/10
  
  d2 <- den %>%  group_by(Rep) %>% mutate(nt=ntile(Rep,lInt)) %>%  group_by(Rep,nt) %>% summarize(meanH=mean(H),sdH=sd(H),meanRich=mean(Richness),sdRich=sd(Richness)) %>% 
    mutate(nt=nt*100)

  if(!simul) {
    if(sims>10)  
      d3 <- filter(d2,Rep %in%  sample(1:sims,10)) 
    else
      d3 <- d2
    
    print(ggplot(d3, aes(x=nt, y=meanH, color=factor(Rep))) +
              geom_errorbar(aes(ymin=meanH-sdH, ymax=meanH+sdH), width=.1,colour="gray") +
              geom_line() + theme_bw() + ggtitle(paste(side,repl)))

    print(ggplot(d3, aes(x=nt, y=meanRich, color=factor(Rep))) +
              geom_errorbar(aes(ymin=meanRich-sdRich, ymax=meanRich+sdRich), width=.1,colour="gray") +
              geom_line() + theme_bw() + ggtitle(paste(side,repl)))
  }

  d2 <- d2 %>%  group_by(Rep) %>% summarise_each(funs(last))

  return(data.frame(MetaNsp=nsp,Side=side,Disp=disp,Migr=migr,Repl=repl,MetaType=meta,
             TMaxH=d1$TMaxH,MaxH=d1$MaxH,TMaxRich=d1$TMaxRich,MaxRich=d1$MaxRich,
             meanH=d2$meanH,sdH=d2$sdH,meanRich=d2$meanRich,sdRich=d2$sdRich,
             meanEven=d2$meanH/log(d2$meanRich)))
  
}

# Plot average H and Richnes of neutral simulations (Time=2900-3000) 
# 
plot_simul_timeRH <- function(CT,mct)
{
  require(ggplot2)
  print(ggplot(CT, aes(x=Repl, y=meanRich)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1)) + facet_grid(MetaType ~ . )) 
  print(ggplot(CT, aes(x=Repl, y=meanH)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))
  print(ggplot(CT, aes(x=Repl, y=meanEven)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))
  
  print(ggplot(mct, aes(x=Repl, y=meanH)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))
  print(ggplot(mct, aes(x=Repl, y=meanRich)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . )) 
  print(ggplot(mct, aes(x=Repl, y=meanEven)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . )) 
  
  print(ggplot(CT, aes(x=Repl, y=MaxH)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))
  print(ggplot(CT, aes(x=Repl, y=MaxRich)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))

  print(ggplot(mct, aes(x=Repl, y=MaxH)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))
  print(ggplot(mct, aes(x=Repl, y=MaxRich)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))

  print(ggplot(CT, aes(x=Repl, y=TMaxH)) + geom_point() + theme_bw() +scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))
  print(ggplot(CT, aes(x=Repl, y=TMaxRich)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))

  print(ggplot(mct, aes(x=Repl, y=TMaxH)) + geom_point() + theme_bw() +scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))
  print(ggplot(mct, aes(x=Repl, y=TMaxRich)) + geom_point() + theme_bw() + scale_x_log10(breaks=c(0.003,0.02,0.1,1))  + facet_grid(MetaType ~ . ))

}


# Simulations of the model with output of one time 
# if delo=T make new simulations and delete old ones
#    delo=F read files generated by previous simulations
#
simulNeutral_1Time <- function(nsp,side,disp,migr,repl,clus="S",time=1000,sims=10,simulate=T,mf="N",meta="U",delo=T) 
{

  if(toupper(meta)=="L") {
    prob <- genFisherSAD(nsp,side)
    neuParm <- paste0("fishP",nsp,"_",side,"R", repl)
    bname <- paste0("neuFish",nsp,"_",side,"R", repl)
  } else {
    prob <- rep(1/nsp,nsp)  
    neuParm <- paste0("unifP",nsp,"_",side,"R", repl)
    bname <- paste0("neuUnif",nsp,"_",side,"R", repl)
  }
  pname <- paste0("pomacR",repl,".lin")
  bname <- paste0(bname ,"T", time )
  

  # make simulations 
  if(simulate)
  {
    # Delete old simulations
    if(delo)
      system(paste0("rm ",bname,"*.txt"))

    genNeutralParms(neuParm,side,prob,1,0.2,disp,migr,repl)

    par <- read.table("sim.par",quote="",stringsAsFactors=F)

    par[par$V1=="nEvals",]$V2 <- time
    par[par$V1=="inter",]$V2 <- time          # interval to measure Density and Diversity
    par[par$V1=="init",]$V2 <- time           # Firs time of measurement = interval
    par[par$V1=="modType",]$V2 <- 4           # Hierarchical saturated
    par[par$V1=="sa",]$V2 <- "N"              # Save a snapshot of the model
    par[par$V1=="baseName",]$V2 <- bname      # Base name for output
    par[par$V1=="mfDim",]$V2 <- mf
    par[par$V1=="minBox",]$V2 <- 2
    par[par$V1=="pomac",]$V2 <- 1             # 0:one set of parms 
                                              # 1:several simulations with pomac.lin parameters 
    par[par$V1=="pomacFile",]$V2 <- pname     # pomac file name 
    par[par$V1=="minProp",]$V2 <- 0           # Calculate diversity for all species
    par[par$V1=="clusters",]$V2 <- clus       # Calculate max cluster size
      
    parfname <- paste0("sim",nsp,"_",side,"R", repl,".par")
    write.table(par, parfname, sep="\t",row.names=F,col.names=F,quote=F)

    genPomacParms(pname,1,c(0.2),disp,migr,repl,sims)

    # copy pomExp.lin to pomac.lin
    #system("cp pomExp.lin pomac.lin")
    s <- system("uname -a",intern=T)
    if(grepl("i686",s)) {
      system(paste(neuBin,parfname,paste0(neuParm,".inp")))
    } else {
      system(paste(neuBin64,parfname,paste0(neuParm,".inp")))
    }
  }

  #den <- meltDensityOut_NT(bname,nsp)
  require(plyr)
  require(dplyr)
  #den <-mutate(den, Species=substring(as.character(Species),2))
  #den <- group_by(den, rep) %>% select(Species,value) %>% top_n(n=1) 
  
  den <-readWideDensityOut(bname)
  den <-den[, c("GrowthRate","MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate","Time", 
                "Richness","H")]
  clu <-readClusterOut(bname)
  if(nrow(den)<nrow(clu))
    clu <- clu[1:nrow(den),]
  clu$Richness <- den$Richness 
  clu$H <- den$H
  clu <-rename(clu,MaxSpeciesAbund=TotalSpecies) %>%
               mutate( MetaNsp=nsp, Side=side, MetaType=as.character(meta),MaxClusterProp = MaxClusterSize/(side*side),MaxClusterSpProp=MaxClusterSize/MaxSpeciesAbund,
                SpanningClust=ifelse(SpanningSpecies>0,MaxClusterProp,0))

  
  clu <- clu[,c(14:16,1:13,17:19)]
  return(clu)
}

# Simulates up to time and saves a snapshot of the model 
#
#
simul_NeutralSAD <- function(nsp,side,time,meta="L",ReplRate=c(0,0.001,0.01,0.1,1)) {
  if(!exists("neuBin")) stop("Variable neuBin not set (neutral binary)")
  if(!require(untb))  stop("Untb package not installed")

  if(toupper(meta)=="L") {
    prob <- genFisherSAD(nsp,side)
    sadName <- "Neutral"
  } else {
  
  if(toupper(meta)=="L") {
    prob <- genFisherSAD(nsp,side)
    neuParm <- paste0("fishE",nsp,"_",side)
    bname <- paste0("neuFish",nsp,"_",side)
    sadName <- "Neutral"
  } else {
    prob <- rep(1/nsp,nsp)  
    neuParm <- paste0("unifE",nsp,"_",side)
    bname <- paste0("neuUnif",nsp,"_",side)
    sadName <- "NeuUnif"
  }
    sadName <- "NeuUnif"
  }
  sad <- data.frame()
  for(i in 1:length(ReplRate)) {
    if(toupper(meta)=="L") {
      neuParm <- paste0("fishE",nsp,"_",side,"R", ReplRate[i])
      bname <- paste0("neuFish",nsp,"_",side,"R", ReplRate[i])
    } else {
      neuParm <- paste0("unifE",nsp,"_",side,"R", ReplRate[i])
      bname <- paste0("neuUnif",nsp,"_",side,"R", ReplRate[i])
    }

    genNeutralParms(neuParm,side,prob,1,0.2,0.4,0.001,ReplRate[i])

    # Delete old simulations
    system(paste0("rm ",bname,"*"))

    par <- read.table("sim.par",quote="",stringsAsFactors=F)
    # Change base name
    par[par$V1=="nEvals",]$V2 <- time
    par[par$V1=="inter",]$V2 <- time # interval to measure Density and Diversity
    par[par$V1=="init",]$V2 <- time  # Firs time of measurement = interval
    par[par$V1=="modType",]$V2 <- 4 # Hierarchical saturated
    par[par$V1=="sa",]$V2 <- "S" # Save a snapshot of the model
    par[par$V1=="baseName",]$V2 <- bname 
    par[par$V1=="minBox",]$V2 <- 2
    par[par$V1=="pomac",]$V2 <- 0 # 0:one set of parms 
                                  # 1:several simulations with pomac.lin parameters 

    parfname <- paste0("sim",nsp,"_",side,"_",ReplRate[i],".par")
    write.table(par, parfname, sep="\t",row.names=F,col.names=F,quote=F)
    s <- system("uname -a",intern=T)
    if(grepl("i686",s)) {
      system(paste(neuBin,parfname,paste0(neuParm,".inp")))
    } else {
      system(paste(neuBin64,parfname,paste0(neuParm,".inp")))
    }

    #fname <- paste0("neuFish",nsp,"Density.txt")
    #sad1 <- meltDensityOut_NT(fname,nsp)

    
    fname <- paste0(bname,"-",formatC(time,width=4,flag=0),".sed")
    spa <- read_sed(fname)

    sad1 <- data.frame(table(spa),Type="SAD",Side=side,NumSp=nsp,SAD=sadName,RplRt=ReplRate[i])
    sad <- rbind(sad,sad1)
  }
  return(sad)
}







# Plot of multiespecies spatial pattern generated with diffent SAD
# type: U=uniform SAD L=Logseries SAD
# nsp: no of species
# side: side of the image
#
plotSAD_SpatPat<-function(nsp,side,type="U")
{
  require(ggplot2)
  #nsp<-64
  #side<-256
  if(type=="U") {
    fname <- paste0("unif",nsp,"_",side,".sed")
    fname1 <- paste0("unif",nsp,"_",side,"rnz.sed") 
  } else {
    fname <- paste0("fisher",nsp,"_",side,".sed")
    fname1 <- paste0("fisher",nsp,"_",side,"rnz.sed")
  }
  
  spa <-read_sed2xy(fname)
  spa$Type <- "a) Regular"
  #spa <- rbind(spa,spa[nrow(spa),])

  sp1 <- read_sed2xy(fname1)
  sp1$Type <- "b) Randomized"

  spa <-  rbind(spa,sp1)
  spa$SAD = ifelse(type=="U","Uniform","Logseries")

  if(type=="B") {
    fname <- paste0("unif",nsp,"_",side,".sed")
    fname1 <- paste0("unif",nsp,"_",side,"rnz.sed") 

    sp1 <-read_sed2xy(fname)
    sp1$Type <- "a) Regular"
    #spa <- rbind(spa,spa[nrow(spa),])

    sp2 <- read_sed2xy(fname1)
    sp2$Type <- "b) Randomized"

    sp1 <- rbind(sp1,sp2)   
    sp1$SAD <- "Uniform"

    spa <- rbind(spa,sp1)
  }

  #g <- ggplot(spa, aes(x, y, fill = factor(v))) + geom_raster(hjust = 0, vjust = 0) + 
  #  theme_bw() + coord_equal() 
  # g <- g + scale_fill_grey(guide=F) +
  g <- ggplot(spa, aes(x, y, fill = v)) + geom_raster(hjust = 0, vjust = 0) + 
    theme_bw() + coord_equal() 
  
  g <- g + scale_fill_gradient(low="red", high="green", guide=F) +
    scale_x_continuous(expand=c(.01,.01)) + 
    scale_y_continuous(expand=c(.01,.01)) +  
    labs(x=NULL, y=NULL) 

  if(type=="B")
  {
    g <- g + facet_grid(SAD ~ Type)
  } else {
    g <- g + facet_grid(. ~ Type) 
  }
  #g <- ggplot(spa, aes(x, y, fill = v)) + geom_raster(hjust = 0, vjust = 0) + theme_bw() + coord_equal() + facet_grid(. ~ Type)
  #g <- g + scale_fill_gradient(low="red", high="green", guide=F) +
  #    labs(x=NULL, y=NULL) 
}



plotDq_Side_Sp <- function(Dqq,side,nsp,sad="Uniform"){
  require(dplyr)
  require(ggplot2)
  require(grid)
  #mylabs <- list(bquote(D[q]^SRS),bquote(Rnz -~D[q]^SRS),bquote(D[q]^SAD),bquote(Rnz -~D[q]^SAD))
  mylabs <- list("Regular","Randomized")
  
  if(nsp!=0) {
    if(sad=="B"){
        mylabs <- list("Regular\nUniform","Randomized\nUniform","Regular\nLogseries","Randomized\nLogseries" )
        Dq1<- filter(Dqq,Side==side,NumSp==nsp,SAD=="Uniform" | SAD=="Logseries")
    } else
        Dq1<- filter(Dqq,SAD==sad,Side==side,NumSp==nsp)

    Dq1 <- group_by(Dq1,SAD,Type,q) %>% summarize(SD.Dq=sd(Dq),Dq=mean(Dq),count=n()) %>% 
#      mutate(DqType=ifelse(grepl("SRS",Type),"SRS","SAD"), SAD1=ifelse(grepl("Logseries",SAD),"a) Logseries","b) Uniform"))
      mutate(DqType=ifelse(grepl("SRS",Type),"SRS","SAD"), TypeSAD=paste(Type,SAD) )

    g <- ggplot(Dq1, aes(x=q, y=Dq, shape=TypeSAD)) +
              geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1,colour="gray") +
              geom_point() + theme_bw() + ylab(expression(D[q]))
    g <- g  + scale_shape_manual(values=c(21,24,21,24,3,4,3,4),guide=guide_legend(title=NULL),
                                                               breaks=c("DqSRS Uniform","rnzDqSRS Uniform","DqSRS Logseries","rnzDqSRS Logseries"),
                                                               labels=mylabs) 
    if(sad=="B")
      g <- g + facet_wrap( ~ DqType, scales="free")
    else
      g <- g + facet_wrap(~ DqType, scales="free")

  } else {
    if(sad=="B") {
      mylabs <- list("Regular\nUniform","Randomized\nUniform","Regular\nLogseries","Randomized\nLogseries" )
      Dq1<- filter(Dqq,Side==side,SAD=="Uniform" | SAD=="Logseries")
    } else
      Dq1<- filter(Dqq,SAD==sad,Side==side)

    Dq1 <- group_by(Dq1,SAD,Type,NumSp,q) %>% summarize(SD.Dq=sd(Dq),Dq=mean(Dq),count=n()) %>% 
#    mutate(DqType=ifelse(grepl("SRS",Type),"SRS","SAD"))
    mutate(DqType=ifelse(grepl("SRS",Type),"DqSRS","DqSAD"), TypeSAD=paste(Type,SAD), NumSpecies=paste("No. Species",NumSp))
    Dq1$NumSpecies <- factor(Dq1$NumSpecies, levels = c("No. Species 8", "No. Species 64", "No. Species 256"))

    g <- ggplot(Dq1, aes(x=q, y=Dq, shape=TypeSAD,colour=TypeSAD)) +
              geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1,colour="gray") +
              geom_point(size=1.3) + theme_bw() + ylab(expression(D[q]))
    g <- g  + scale_shape_manual(values=c(21,24,21,24,3,4,3,4),guide=guide_legend(title=NULL),
              breaks=c("DqSRS Uniform","rnzDqSRS Uniform","DqSRS Logseries","rnzDqSRS Logseries"),labels=mylabs)

    library(RColorBrewer)
    mc <- brewer.pal(5, "Set1")
    g <- g + scale_colour_manual(values=c(mc[1],mc[2],mc[1],mc[2],mc[3],mc[4],mc[3],mc[4]),
        guide=guide_legend(title=NULL),
        breaks=c("DqSRS Uniform","rnzDqSRS Uniform","DqSRS Logseries","rnzDqSRS Logseries"),labels=mylabs) 
 
    g <- g + facet_wrap(NumSpecies~ DqType, scales="free",ncol=2) + theme(legend.key.size = unit(1, "cm"))
  }
}



readDq_fit <- function(side,nsp,sad="U") {
  if(sad=="Uniform") {
    fname <- paste0("unif",nsp,"_",side,".sed")
    fname1 <- paste0("unif",nsp,"_",side,"rnz.sed") 
  } else {
    fname <- paste0("fisher",nsp,"_",side,".sed")
    fname1 <- paste0("fisher",nsp,"_",side,"rnz.sed")
  }  

  Dq1<- calcDq_multiSBA(fname,"q.sed 2 1024 20 S",mfBin,T)
  zq <- readZq(paste0("t.", fname),"q.sed")
  zq$Type <- "DqSRS"

  Dq1<- calcDq_multiSBA(fname1,"q.sed 2 1024 20 S",mfBin,T)
  zq1 <- readZq(paste0("t.", fname1),"q.sed")
  zq1$Type <- "rnzSRS"
  zq<- rbind(zq,zq1)

  Dq1<- calcDq_multiSBA(fname,"q.sed 2 1024 20 E",mfBin,T)
  zq1 <- readZq(paste0("t.", fname),"q.sed")
  zq1$Type <- "DqSAD"
  zq<- rbind(zq,zq1)

  Dq1<- calcDq_multiSBA(fname1,"q.sed 2 1024 20 E",mfBin,T)
  zq1 <- readZq(paste0("t.", fname1),"q.sed")
  zq1$Type <- "rnzDqSAD"
  zq<- rbind(zq,zq1)
}


# Plot of multiespecies spatial pattern generated with neutral model and logseries SAD
#
#
plotNeutral_SpatPat<-function(nsp,side,time,meta="L",ReplRate=c(0,0.001,0.01,0.1,1))
{
  require(ggplot2)
 
  spa <- data.frame()
  for(i in 1:length(ReplRate)) {
    if(toupper(meta)=="L") {
      bname <- paste0("neuFish",nsp,"_",side,"R", ReplRate[i])
    } else {
      bname <- paste0("neuUnif",nsp,"_",side,"R", ReplRate[i])
    }
    
    fname <- paste0(bname,"-",formatC(time,width=4,flag=0),".sed")
    
    sp1 <-read_sed2xy(fname)
    sp1$Type <- paste("Replacement:", ReplRate[i])
    sp1$Species <- paste("Species:",length(unique(sp1$v)))

    spa <-  rbind(spa,sp1)
  }
  #lvl <- unique(spa$Type)
  #spa$Type <- factor(spa$Type, levels = lvl)
  
#  g <- ggplot(spa, aes(x, y, fill = factor(v))) + geom_raster(hjust = 0, vjust = 0) + 
  g <- ggplot(spa, aes(x, y, fill = v)) + geom_raster(hjust = 0, vjust = 0) + 
    theme_bw() + coord_equal() 
  
  g <- g + scale_fill_gradient(low="red", high="green", guide=F) +
#  g <- g + scale_fill_grey(guide=F) +
    scale_x_continuous(expand=c(.01,.01)) + 
    scale_y_continuous(expand=c(.01,.01)) +  
    labs(x=NULL, y=NULL) 
  
  g <- g + facet_wrap( ~ Type +Species,ncol=2) 
  print(g)
  
  #g <- ggplot(spa, aes(x, y, fill = v)) + geom_raster(hjust = 0, vjust = 0) + theme_bw() + coord_equal() + facet_grid(. ~ Type)
  #g <- g + scale_fill_gradient(low="red", high="green", guide=F) +
  #    labs(x=NULL, y=NULL) 
}



# Plot of Dq multiespecies spatial pattern generated with neutral model and logseries SAD
#
#
plotNeutral_Dq<-function(nsp,side,time,meta="L",ReplRate=c(0,0.001,0.01,0.1,1))
{
  require(ggplot2)
  if(!is.null(ReplRate))
  { 
    Dq3 <- data.frame()
    for(i in 1:length(ReplRate)) {
      if(toupper(meta)=="L") {
        bname <- paste0("neuFish",nsp,"_",side,"R", ReplRate[i])
      } else {
        bname <- paste0("neuUnif",nsp,"_",side,"R", ReplRate[i])
      }
      
      fname <- paste0(bname,"-",formatC(time,width=4,flag=0),".sed")

      Dq1<- calcDq_multiSBA(fname,"q.sed 2 1024 20 S",mfBin,T)
      Dq1$DqType <- "DqSRS"
      Dq1$ReplacementRate <- ReplRate[i]
      Dq2<- calcDq_multiSBA(fname,"q.sed 2 1024 20 E",mfBin,T)
      Dq2$DqType <- "DqSAD"
      Dq2$ReplacementRate <- ReplRate[i]

      Dq3 <- rbind(Dq3,Dq1,Dq2)
    }
  } else {
    if(nsp==0){
      Dq3 <- data.frame()
      for(num in c(8,64,256))
      {
        Dq3 <- rbind(Dq3, plotNeutral_Dq_aux(num,side))  
      }
      Dq3 <- mutate(Dq3,spMeta = paste0("Metacommunity sp.",spMeta))
      Dq3$spMeta <- factor(Dq3$spMeta,levels=unique(Dq3$spMeta))
    } else {
      Dq3 <- plotNeutral_Dq_aux(nsp,side)  
    }
  }

  g <- ggplot(Dq3, aes(x=q, y=Dq, shape=factor(ReplacementRate),colour=factor(ReplacementRate))) +
            geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1,colour="gray") +
            geom_point(size=1.3) + theme_bw() + ylab(expression(D[q]))

  library(RColorBrewer)
  mc <- brewer.pal(6, "Set1")
  g <- g + scale_colour_manual(values=mc,guide=guide_legend(title="Replacement")) 
  g <- g + scale_shape_manual(values=c(21,24,4,25,3,8),guide=guide_legend(title="Replacement")) 
  
  if(nsp==0){
    g <- g + facet_wrap(spMeta ~ DqType, scales="free",ncol=2)
  } else {
    g <- g + facet_wrap(~ DqType, scales="free",ncol=2)
  }
  print(g)
}

plotNeutral_Dq_aux<-function(nsp,side,time=500,meta="L")
{
  if(toupper(meta)=="L") {
    bname <- paste0("neuFish",nsp,"_",side)
  } else {
    bname <- paste0("neuUnif",nsp,"_",side)
  }
  fname <- paste0(bname,"T",time,"mfOrd.txt")
  Dq1 <- readNeutral_calcDq(fname)
  Dq1$DqType <- "DqSRS"

  fname <- paste0(bname,"T",time,"mfSAD.txt")
  Dq2 <- readNeutral_calcDq(fname)
  Dq2$DqType <- "DqSAD"
  
  Dq1 <- rbind(Dq1,Dq2)
  require(dplyr)
  
  Dq3 <- filter(Dq1,MortalityRate==.2,DispersalDistance==0.4,ColonizationRate==0.001) %>% 
    group_by(DqType,ReplacementRate,q) %>% summarize(SD.Dq=sd(Dq),Dq=mean(Dq),count=n()) 
  Dq3$spMeta <- ceiling(as.numeric(nsp)*1.33)
  return(Dq3)
} 

R2Neutral_Dq<-function(side,time=500,meta="L")
{
  require(dplyr)
  hh <-function(x,c,...){
    length(x[x>c])/length(x)
  }

  Dq3 <- data.frame()
  for(nsp in c(8,64,256))
  {
    if(toupper(meta)=="L") {
      bname <- paste0("neuFish",nsp,"_",side)
    } else {
      bname <- paste0("neuUnif",nsp,"_",side)
    }
    fname <- paste0(bname,"T",time,"mfOrd.txt")
    Dq1 <- readNeutral_calcDq(fname)
    Dq1$DqType <- "DqSRS"

    fname <- paste0(bname,"T",time,"mfSAD.txt")
    Dq2 <- readNeutral_calcDq(fname)
    Dq2$DqType <- "DqSAD"
    
    Dq1 <- rbind(Dq1,Dq2)
    require(dplyr)
    
    Dq1 <- filter(Dq1,MortalityRate==.2,DispersalDistance==0.4,ColonizationRate==0.001) %>% 
      group_by(DqType,ReplacementRate) %>% 
      summarize(Freq60=hh(R.Dq,.6),Freq90=hh(R.Dq,.9)) 

    Dq1$spMeta <- ceiling(as.numeric(nsp)*1.33)


    fname <- paste0(bname,"T",time,"Density.txt")
    den1 <- meltDensityOut_NT(fname,unique(Dq1$spMeta))
    den1 <- filter(den1,MortalityRate==.2,DispersalDistance==0.4,ColonizationRate==0.001) %>% 
      group_by(rep,ReplacementRate) %>% summarize(nsp=n()) %>%
      group_by(ReplacementRate) %>% summarize(meanSp=mean(nsp))
    Dq1 <- left_join(Dq1,den1)   

    Dq3 <- rbind(Dq3,Dq1)
    
  }
  Dq3$Side <- side
  return(Dq3)
} 


# Plot RAD from simulations with only 1T output
#  local=T : plots the average of local RAD
#        F : plots metacommunity RAD
#
plotNeutral_SAD_1T<-function(nsp,side,repl,time,meta="U",local=T,cpal="")
{
  require(ggplot2)
  require(plyr)
  require(dplyr)
  den<- data.frame()
  if(nsp!=0){
    den1 <- data.frame()
    for(r in repl){
      if(toupper(meta)=="L") {
        bname <- paste0("neuFish",nsp,"_",side,"R", r)
      } else {
        bname <- paste0("neuUnif",nsp,"_",side,"R", r)
      }
      fname <- paste0(bname ,"T", time, "Density.txt")

      if(local) {
        den <- meltDensityOut_NT(fname,nsp)
        den <- calcRankSAD_by(den,"value",1:5)
        den <- group_by(den,ReplacementRate,Rank) %>% summarize(Freq=mean(value),count=n()) 
        den <- calcRankSAD_by(den,"Freq",1:1)
      } else {
        den <- group_by(den,ReplacementRate,Species) %>% summarize(Freq=mean(value),count=n()) 
        den <- calcRankSAD_by(den,"Freq",1:1)
      }
      den1 <- rbind(den1,den)
    }
    
    g <- ggplot(den1,aes(x=Rank,y=log(Freq),colour=factor(ReplacementRate))) +  theme_bw() + geom_point(size=1,shape=c(21))
    #library(RColorBrewer)
    #mc <- brewer.pal(6, "Set1")
    if(cpal=="")
      g <- g + scale_colour_discrete(name=bquote("  "~rho))  + geom_line()
    else
      g <- g + scale_colour_manual(values=cpal,name=bquote("  "~rho))  + geom_line()



#     g <- ggplot(den1,aes(x=Rank,y=log(Freq),colour=ReplacementRate)) +  theme_bw() + geom_point(size=1)
#     mc <- c("#b35806","#e08214","#fdb863","#fee0b6","#f7f7f7","#d8daeb","#b2abd2","#8073ac","#542788")
#     
#     g <- g + scale_colour_gradientn(colours=mc, guide="legend",trans="log",name=bquote("  "~rho)) +  geom_line()
    
  } else {
    for(nsp in c(8,64,256))
    {
      den <-rbind(den,plotNeutral_SAD_aux(nsp,side))
    }
  
    den <- mutate(den, metaLbl =paste0("Metacommunity sp.",nsp))
    ml <- unique(den$metaLbl)
    den$metaLbl <- factor(den$metaLbl,levels=c(ml[1],ml[2], ml[3]))
    g <- ggplot(den,aes(x=Rank,y=log(Freq),shape=factor(ReplacementRate),colour=factor(ReplacementRate))) +  theme_bw() + geom_point(size=1)
    library(RColorBrewer)
    mc <- brewer.pal(6, "Set1")
    g <- g + scale_colour_manual(values=mc,name=bquote("  "~rho)) 

    g <- g + scale_shape_manual(values=c(21,24,4,25,3,8),name=bquote("  "~rho)) +
      geom_smooth(se=F,span = 0.70) 
    g <- g + facet_wrap(~ metaLbl, scales="free",ncol=2)
  }
  print(g)
  return(den1)
}

plotNeutral_SAD_aux<-function(nsp,side,time=500,meta="L")
{
  if(toupper(meta)=="L") {
    bname <- paste0("neuFish",nsp,"_",side)
  } else {
    bname <- paste0("neuUnif",nsp,"_",side)
  }
  fname <- paste0(bname,"T",time,"Density.txt")
  spMeta <- ceiling(as.numeric(nsp)*1.33)

  den1 <- meltDensityOut_NT(fname,spMeta)

  den3 <- filter(den1,MortalityRate==.2,DispersalDistance==0.4,ColonizationRate==0.001) 
  den3 <- calcRankSAD_by(den3,"value",1:5)
  
  
  den3 <- group_by(den3,ReplacementRate,Rank) %>% summarize(Freq=mean(value),count=n()) 
  den3$spMeta <- spMeta
  return(den3)  
}

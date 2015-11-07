
# uses the CTFS R formated Tables

# From x,y,species,dbh data Generates a lattice (sed file) with species for dbh>=dbhmin 
# and calculates the clusters size distribution 
#
# bcidata       = dataframe with x,y data
# biomassXYname = name base for output files
# dbhmin        = min dbh for processing
#
# Result:       A list with H (Shannon Diversity) and S (SPecies richnes) before Hin,Sin 
#               and after the discretization Hout,Sout 
#               ad_p: the pvalue of Anderson-Darling test to check if SADs are different
#               a data.frame with cluster sizes and species. the firs record signals if an spannig cluster exists.
#
genSpSed_hk<- function(bcidata, biomassXYname, sizex,sizey,lc=1,dbhmin=0)
{
  require(kSamples)
  bci <-with(bcidata, subset(bcidata, status=='A' & !is.na(gx) & dbh>=dbhmin,select=c('gx','gy','sp','dbh')))
  bcin <- as.data.frame(table(bci$sp))
  bcin <-cbind(bcin[rev(order(bcin$Freq)),], seq(nrow(bcin)))
  names(bcin)[3]<-"SpNum"
  names(bcin)[1]<-"sp"
  bci$SpNum <- bcin$SpNum[match(bci$sp,bcin$sp)]
  
  require(vegan)
  answer <- list(Hin=diversity(bcin$Freq),Sin=nrow(bcin))

  fname.txt <- paste(biomassXYname ,".txt",sep="")
  fname.sed <- paste(biomassXYname ,".sed",sep="")
  fname.bin <- paste(biomassXYname ,".pat",sep="")
  system(paste0("rm ", fname.sed))
  system(paste0("rm ", fname.bin))
  
  write.table(bci[,c(1,2,5,4)],fname.txt,col.names=F,row.names=F)

  s <- system("uname -a",intern=T)
  if(!grepl("i686",s)) {
    ClusterizeBin <-paste0(ClusterizeBin,"64")
    hkBin <-paste0(hkBin,"64")
  }
  # Ejecuta XY2sed para generar archivo sed
  syst.txt <- paste(ClusterizeBin, fname.txt, fname.sed, sizex,lc,sizey,lc,"2")
  #syst.txt <- paste("./XY2sed", biomassXYname.txt, biomassXYname.sed, "1000 1 500 1 2")
  system(syst.txt)
  bci <- read_sed(fname.sed)
  
  bci1 <-as.data.frame(table(bci[bci>0]))
  
  ks <-ad.test(list(bcin$Freq,bci1$Freq),method="simulated",nsim=1000)
  
  answer <- c(answer,Hout=diversity(bci1$Freq),Sout=nrow(bci1),ad_p=ks$ad[[7]])
  # Runs hk to calculate clusters
  #

  syst.txt <- paste(hkBin,fname.sed,fname.bin, "SP")
  system(syst.txt)

  clu <-readClusterHk(fname.bin)  
  
  return( c(answer,Clusters=clu) )
  
  # Leer archivo 
}

# Read cluster sizes from hk program output 
#  
#
readClusterHk <- function(fname){
  den <- read.delim(fname,stringsAsFactors = F)
  return(den)
}



# Genera archivo sed con especie + abundante y calcula espectro MF 
#
genSpAbundSed_mfSBA <- function(bcidata, biomassXYname, dbhmin=0)
{
  attach(bcidata)
  # dbh >100 mm
  bci <- subset(bci.full1, status=='A' & !is.na(gx) & dbh>=dbhmin,select=c('gx','gy','sp','agb'))
  detach(bcidata)
  attach(bci)
  bcin <- as.data.frame(table(sp))
  bcin <-cbind(bcin[rev(order(bcin$Freq)),], seq(nrow(bcin)))
  names(bcin)[3]="SpNum"
  bci$SpNum <- bcin$SpNum[match(bci$sp,bcin$sp)]
  
  library(vegan)
  answer <- list(diversity(bcin$Freq),nrow(bcin))
  names(answer)[[1]]<-"Shannon H"
  names(answer)[[2]]<-"Number of species S"
  
  biomassXYname.txt <- paste(biomassXYname ,".txt",sep="")
  biomassXYname.sed <- paste(biomassXYname ,".sed",sep="")
  
  write.table(bci[bci$SpNum==1,c(1,2,5,4)],biomassXYname.txt,col.names=F,row.names=F)
  detach(bci)
  
  # Ejecuta XY2sed para generar archivo sed
  syst.txt <- paste("./XY2sed", biomassXYname.txt, biomassXYname.sed, "1000 1 500 1 2")
  system(syst.txt)
  
  # Ejecuta mfSBA para espectro multifractal 
  syst.txt <- paste("./mfSBA", biomassXYname.sed, "q21.sed 2  512 20 S")
  system(syst.txt)
  return(answer)
}


# 
# Fit 3 models using ML for patch sizes, exponential, power, and power with exponential cutoff
# and calculates DeltaAICc
# 

#
# return a data.frame with a field type: 
#  1=Spanning
#  2=Most Abundant
#  3=Other not most Abundant
#  4=Other not Spanning
#
BCIFit_Clusters <- function(clu) 
{

  require(plyr)
  require(dplyr)
  
 
  # Spanning species
  clu1 <-clu %>% slice(1:1) %>% rename(Spanning=ClusterSize)
  
  clu <- clu %>% slice(2:n()) %>% arrange(desc(ClusterSize))
  
  clu <- left_join(clu,clu1) %>% mutate(Spanning=ifelse(is.na(Spanning),0,1))
  
  # Select most abundant or spanning species !!!!!!!
  #
  clu1 <-clu %>% slice(1:1) %>% select(Species,Spanning)
  
  
  #if(clu1$Spanning[1]==0)
  #  clu1 <- clu %>%  group_by(Species) %>% summarize(totClu=sum(ClusterSize),Spanning=mean(Spanning)) %>% filter( totClu==max(totClu)) %>% select(Species,Spanning)
  
  clu2 <- inner_join(clu,clu1)
  
  # Fits only 1 species patches, most abundant or spanning
  #
  mdl <- fitNeutral_Clusters(clu2)
  mdl$Species <- clu$Species[1]
  mdl$ClusterSize <- clu$ClusterSize[1]
      
  clu2 <- anti_join(clu,clu1) # eliminates most abundant or spanning sp

  # if there is no spanning sp change Spanning =3, else Spanning =4
  #
  clu2 <- clu2 %>% mutate(Spanning=ifelse(clu1$Spanning[1]==0,3,4))
  mdl1 <- fitNeutral_Clusters(clu2)
  mdl1$Species <- 0
  mdl1$ClusterSize <- 0
  

  mdl <- bind_rows(mdl,mdl1)
  
  # Calculates DeltaAIC
  #
  
  mdl <- group_by(mdl,type) %>% mutate( DeltaAIC= AICc -min(AICc)) %>% arrange(DeltaAIC)
}


# Plot of cluster sizes
# clu: CLusters
# mdl: fitted models
#
BCIplot_Clusters <- function(clu,mdl,year) 
{
  options("scipen"=0, "digits"=4)
  
  require(plyr)
  require(dplyr)
  
  # Spanning species
  clu1 <-clu %>% slice(1:1) %>% rename(Spanning=ClusterSize)
  
  clu <- clu %>% slice(2:n()) %>% arrange(desc(ClusterSize))
  
  clu <- left_join(clu,clu1) %>% mutate(Spanning=ifelse(is.na(Spanning),0,1))
  
  # Select most abundant or spanning species !!!!!!!
  #
  clu1 <-clu %>% slice(1:1) %>% select(Species,Spanning)

  tpe<-1
  if(clu1$Spanning[1]==0)
  {
    tdes <- "Max Patch"
    tpe <- 2
  }
  else 
    tdes <- "Spanning"

  clu2 <- inner_join(clu,clu1)
  
  
  # Plots only 1 species patches, most abundant or spanning
  #
  mPow<-displ$new(clu2$ClusterSize)

  # select power law
  m0 <- filter(mdl,model=="Pow",Year==year,type==tpe)
  xmin <- m0$xmin
  mPow$setXmin(xmin)
  mPow$setPars(m0$alfa)
  
  mExp<-disexp$new(clu2$ClusterSize)
  m1 <- filter(mdl,model=="Exp",Year==year,type==tpe)
  mExp$setXmin(xmin)
  mExp$setPars(m1$alfa)
  tit <- paste("BCI",year,tdes)

  m2 <- filter(mdl,model=="PowExp",Year==year,type==tpe)
  freq_plot_displ_exp(clu2$ClusterSize,m0$alfa,m1$alfa,m2$alfa,m2$rate,xmin,tit)
  #cdfplot_displ_exp(clu2$ClusterSize,m0$alfa,m1$alfa,m2$alfa,m2$rate,xmin,tit)

  # Selects not spanning or not most abundant

  clu2 <- anti_join(clu,clu1) # eliminates most abundant or spanning sp

  # if there is no spanning sp change Spanning =3, else Spanning =4
  #
  clu2 <- clu2 %>% mutate(Spanning=ifelse(clu1$Spanning[1]==0,3,4))
  tdes <- ifelse(clu1$Spanning[1]==0,"Other Max Patch","Other Spanning")
  tpe <- ifelse(clu1$Spanning[1]==0,3,4)
  # Plots multi species patches, not most abundant or not spanning
  #
  mPow<-displ$new(clu2$ClusterSize)

  # select power law
  m0 <- filter(mdl,model=="Pow",Year==year,type==tpe)
  xmin <- m0$xmin
  mPow$setXmin(xmin)
  mPow$setPars(m0$alfa)
  
  mExp<-disexp$new(clu2$ClusterSize)
  m1 <- filter(mdl,model=="Exp",Year==year,type==tpe)
  mExp$setXmin(xmin)
  mExp$setPars(m1$alfa)
  tit <- paste("BCI",year,tdes)

  m2 <- filter(mdl,model=="PowExp",Year==year,type==tpe)
  freq_plot_displ_exp(clu2$ClusterSize,m0$alfa,m1$alfa,m2$alfa,m2$rate,xmin,tit)
  #cdfplot_displ_exp(clu2$ClusterSize,m0$alfa,m1$alfa,m2$alfa,m2$rate,xmin,tit)

}

BCI_fit_plot <-function(rdat,dat,spdat,year,x=1000,y=500,lc=1)
{

  setwd(BCIwd)

  load(rdat) 
  dat <- get(dat)
  sa <- genSpSed_hk(dat,spdat,x,y,lc)
  hs <- data.frame(Year=year,Hin=sa$Hin,Hout=sa$Hout,Sin=sa$Sin,Sout=sa$Sout,SpanSp=sa$Clusters.Species[1],
                         SpanClus=sa$Clusters.ClusterSize[1], AD_p=sa$ad_p)

  BCIClus <- data.frame(Species=sa$Clusters.Species,ClusterSize=sa$Clusters.ClusterSize)

  dl <- BCIFit_Clusters(BCIClus) 
  dl$Year <- year

  setwd(oldcd)

  setwd("figure")
  BCIplot_Clusters(BCIClus,dl,year)

  setwd(oldcd)
  return(list(BCIhsize=hs,BCIMdl=dl))
}
# Multifractal patterns for natural communities

## p-model

#I use the following to have only one parameter:
# $p2=p3=p4=(1-p1)/3$
#
#The parameter *iter* (default=9) gives the image size $lado = 2^iter$
#The parameter *rnd* is the how we assign p at each step S->random with reposition

# Calcula espectro multifractal y extrae Dq en un frame
#

calc_pmodel1 <- function(fname,p1,iter=9,rnd="S")
  {
  if(p1==0 | iter<2) stop("p1 == 0 or iter<2") 
  if(file.exists(fname))
    file.remove(fname)
  p2 <- p3 <- p4 <- (1-p1)/3
  syst.txt <- paste("./pmodel", p1,p2,p3,p4,iter,fname,rnd)
  system(syst.txt)
  }

plot_pmodel1 <- function(fname,p1,iter=9)
  {
  require(lattice)
  per <- read.table(fname, skip=2,header=F)
  per <-data.matrix(per)*(2^(2*iter))
  col.l <- colorRampPalette(c('white', 'green', 'purple', 'yellow', 'brown'))(64) 
  levelplot(per, scales = list(draw = FALSE),xlab =NULL, ylab = NULL,col.regions=col.l,
            useRaster=T,
            main=list( paste("p1=",p1),cex=1))
  }

# Plot an sed file
#
# fname: file name
# gname: graph title 
# dX: range in columns to plot
# col: vector of colors to make the color palette
# shf: shift in the vector of colors 
#
plot_image <- function(fname,gname,dX=0,col=0,shf=0)
  {
  require(lattice)
  require(RColorBrewer)
  per <-data.matrix(read.table(fname, skip=2,header=F))
  if(length(dX)>1) per <- per[,dX]
  mp = max(per)
  if(mp<50) {
    mp = 50
    seqat = seq( min(per),max(per),(max(per)-min(per))/50)
  }
  else
  {
    seqat = seq( min(per),mp,5)
  }
  if(length(col)==1) col.l <- colorRampPalette(c('white', 'green', 'purple', 'yellow', 'brown'))(mp) 
  else col.l <- colorRampPalette(col)(mp) 
  if( shf>0) col.l = col.l[shf:mp]
  levelplot(per, scales = list(draw = FALSE),xlab =NULL, ylab = NULL,col.regions=col.l,
            useRaster=T,at=seqat,
            main=list( gname,cex=1))
  }

read_sed <- function(fname)
{
  per <-data.matrix(read.table(fname, skip=2,header=F))
}

# Calculates the multifractal spectrum and 
# put the result in a frame
# 
calcDq_mfSBA <- function(fname,recalc=FALSE)
  {
  sname <- paste0("s.", fname)
  if((!file.exists(sname)) | recalc)
    {
    syst.txt <- paste("./mfSBA ",fname, "q.sed 2 512 20 S")
    system(syst.txt)
    }
  pp <- read.table(sname, header=T)

  pp$Dq  <- with(pp,ifelse(q==1,alfa,Tau/(q-1)))
  pp$SD.Dq  <- with(pp,ifelse(q==1,SD.alfa,abs(SD.Tau/(q-1))))
  pp$R.Dq <- with(pp,ifelse(q==1,R.alfa,R.Tau))
  return(pp[,c("q","Dq","SD.Dq","R.Dq")])
  }

# Calc mfSBA
# draw histogram
# add to a dataframe dq with a new variable site = siten
#
addDqRHist <- function(fname,dq=0,siten="")
{
  # Extract site label
  if( nchar(siten)==0 ){  
    siten <- strsplit(fname,".sed")[[1]][1]
  }
  
  # Make multifractal analysis
  
  dqnew <- calcDq_mfSBA(fname)
  dqnew$Site <- siten
  
  require(ggplot2)
  
  #print(ggplot(dqnew, aes(x=R.Dq)) + geom_histogram() + ggtitle(siten))
  hist(dqnew$R.Dq)
  if(is.data.frame(dq))
    return(rbind(dq,dqnew))
  else
    return(dqnew)
}



repeat_pmodel <- function(p1,fname,reps) {
  
  for(i in 1:reps)
  {
    calc_pmodel1(fname,p1)
    #plot_pmodel1(fname,p1)
    s_Dq <- calcDq_mfSBA(fname,T) # recalc = TRUE
    s_Dq$p1 <- p1
    s_Dq <-na.omit(s_Dq)
    a_Dq <-if(!exists("a_Dq")) s_Dq else rbind(a_Dq,s_Dq)
    }
  #plot_pmodel1(fname,p1)
  return(a_Dq)
}


# Function for boot
#
mediaDq <- function(dq,i) mean(dq[i])

# Calculates CI using boot 
#             sp: frame to put Dq|q
#             the following are to identify the output 
#             p1: the value of p1 used in p-model  
#             type: the type of the run 
calcDq_bootCI <- function(sp,p1,type)
  {
  require(boot)
  require(plyr)
  sp <- na.omit(sp)
  testBoot <- with(sp,by(Dq,q,boot, mediaDq, R=1000))
  # boot pg 331 
  bt <-boot.ci(testBoot[[1]])
  testCI <-  sapply(testBoot, boot.ci)
  compCI1<-data.frame(matrix(unlist(testCI),ncol=21,byrow=T)[,c(2,15,16)])
  names(compCI1)=c("Dq","LowCI","HighCI")
  compCI1$q<-as.numeric(names(testBoot))
  compCI1$p1 <-p1
  compCI1$Type <-type
  compCI1$SD.Dq <-ldply(testBoot, function(e) sd(e$t))$V1
  return(compCI1)
  }


# Función para graficar espectros con CI
# El error dinnames es que no todas las variables son numéricas
# Usa Type como variable de grupo y LowCI,HighCI para el CI
#
plot_DqCI <- function(DqCI, tit="Type")
  {
  require(lattice)
  require(Hmisc)
  with(DqCI,
    xYplot(Cbind(Dq,LowCI,HighCI) ~ q, data=DqCI, nx=FALSE, groups=Type, type="l",
           scales=list(tck=-1),
           cap=.01,
           ylim=c(min(DqCI$LowCI)-.01,max(DqCI$HighCI)+.01),
    	     panel=function(...){
  		     panel.abline(h=2, col="grey", lty=2)
  		     panel.xYplot(...)},
  		     ylab=expression(italic(D[q])),
  		     xlab=expression(italic(q)),
           label.curves=F,
           auto.key=list(x=.65,y=.9,title=tit,cex.title=.9, lines=T,points=F,cex=.7),
           )
       )
  }

# Calcula Dq teorico para pmodel
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

# Funcion para calcular CI a partir del SD 
#
calcCI <- function(dat,n,p) qt(1-p/2,n-1)*sqrt(dat*dat/n)

# Function to plot Dq fit from t* files generated by mfSBA 
# Beware! Is fixed for q interval -5 to 5
#
plotDqFit <- function(fname,wtitle="")
{
  
  zq <- read.table(fname, sep="\t", skip=1,
                   col.names=c("BoxSize","logBox","Xminus5","Xminus4.5","Xminus4","Xminus3.5","Xminus3","Xminus2.5","Xminus2","Xminus1.5","Xminus1","Xminus0.5","X0","X0.5","X1","X1.5","X2","X2.5","X3","X3.5","X4","X4.5","X5")
  )
  #cna=seq(-5,5,by=0.5)
  cna =c(-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)
  zq0=reshape(zq, timevar="q",times=cna,v.names=c("logTr"),
              varying=list(c("Xminus5","Xminus4.5","Xminus4","Xminus3.5","Xminus3","Xminus2.5","Xminus2","Xminus1.5","Xminus1","Xminus0.5","X0","X0.5","X1","X1.5","X2","X2.5","X3","X3.5","X4","X4.5","X5")),
              direction="long")
  
  
  library(lattice)
  zq1 <- subset(zq0, q==1 | q==2 | q==3 | q==4 | q==5 | q==0 | q==-1 | q==-2 | q==-3 | q==-4 | q==-5 )

#  oname <- paste("lsaravia_figS_W",sem,"_",wnro,".tif",sep="")
  
#  tiff(oname, width=4.86,height=4.86,units="in",res=600,compression=c("lzw"))
  #devAskNewPage(TRUE)
  
  trellis.par.set(superpose.symbol=list(pch=c(0,1,2,3,4,5,6,8,15,16,17)))
  trellis.par.set(superpose.symbol=list(cex=c(rep(0.6,11))))
  trellis.par.set(superpose.line=list(lty=3))
  
  #show.settings()
  
  print(xyplot(logTr~logBox , data =zq1, groups=q, type=c("r","p"), scales=list(tck=-1), 
               main=list(wtitle,cex=0.9),
               auto.key=list(space = "right",title=expression(italic("q")),cex.title=.7, points=TRUE,cex=.7),
               ylab=expression(italic(paste("log ",  Z[q](epsilon) ))) , xlab=expression(italic(paste("log ",epsilon)))  
  ))
#  dev.off()    
}


# Calculates Dq from a data.frame read from the output of neutral model 
# auxiliar function for the next one
#
calcDq_frame <- function(pp)
  {
  pp$Dq  <- with(pp,ifelse(q==1,alfa,Tau/(q-1)))
  pp$SD.Dq  <- with(pp,ifelse(q==1,SD.alfa,abs(SD.Tau/(q-1))))
  pp$R.Dq <- with(pp,ifelse(q==1,R.alfa,R.Tau))
  return(pp[,c(1:6,15:17)])
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


# Save a matrix as a sed file with type BI (floating point)
#
save_matrix_as_sed <- function(mat,fname)
{
  header <- paste(nrow(mat),ncol(mat),"BI")
  write.table(header,file=fname,row.names=F,col.names=F,quote=F)
  write.table(mat,file=fname,row.names=F,col.names=F,quote=F,append=T)
}


# Read time series of density richness (S) and Shannon (H) from neutral model
# with replacementRate
#
readNeutral_timeSeries <- function(fname)
{
  pp <- read.table(fname, skip=1,header=F)[,2:10]
  names(pp)<-c("MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate","Time","Tot.Dens","Tot.Num","S","H")
  return(pp)
}

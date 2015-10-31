
# Test clusters sizes read and fit

hkBin <- "~/Dropbox/cpp/SpatialAnalysis/Clusters/hk"

meta <- "L"
Repl <- 0.0002

p <-expand.grid(MetaType=meta,Repl=ReplRate) 

if(toupper(p$MetaType[i])=="L") {
  bname <- paste0("neuFish",nsp,"_",side,"R", p$Repl[i])
} else {
  bname <- paste0("neuUnif",nsp,"_",side,"R", p$Repl[i])
}

fname <- paste0(bname,"-",formatC(time,width=4,flag=0),".sed")

## assume the last is the one with sed
#
clu <-readClusterOut(bname)
spanSp <- clu$SpanningSpecies[nrow(clu)]
if(spanSp==0) # No tiene spanning cluster



#
h21=0.3;  vp1 <- 25
corAB=0.50
corAC=0.90
#
h22=0.1;vp2=100
corBC=0.50
#
h23=0.3;  vp3=100

corE=0


(va1=h21*vp1)
(va2=h22*vp2)
(va3=h23*vp3)

ve1=(vp1-va1)
ve2=(vp2-va2)
ve3=(vp3-va3)

covAB=corAB*sqrt(va1*va2)
covAC=corAC*sqrt(va1*va3)
covBC=corBC*sqrt(va2*va3)

covE=corE*sqrt(ve1*ve2)
#Va <- matrix(c(va1,covA,covA,va2),nrow=2,ncol=2,byrow=TRUE)
#Ve <- matrix(c(ve1,covE,covE,ve2),nrow=2,ncol=2,byrow=TRUE)

Va <- matrix(c(va1,covAB,covAC,covAB,va2,covBC,covAC,covBC,va3),nrow=3,ncol=3,byrow=TRUE)
Ve <- matrix(c(ve1,covE,covE,covE,ve2,covE,covE,covE,ve3),nrow=3,ncol=3,byrow=TRUE)
mu <- c(0.0025,100,1000)


source('mmtrait_simcode.R')
basedata <- makebasepop(nsires=200,ndams=1000,mu,Va,Ve)
cor(basedata[,-1:-5])
var(basedata[,-1:-5])


offdata <- makeoff(Numgen=10,basedata,nsires=200,ndams=1000,ls=3,
                   Va,Ve,sd='tbv/l',md='nested',trsel=1)




dataall <- rbind.data.frame(basedata,offdata)
Rgvals <- as.matrix(dataall[,c('G','TBV1','TBV2','TBV3')])
cor(Rgvals[,-1])
var(Rgvals[,-1])

Reserr <- as.matrix(dataall[,c('G','Res1','Res2','Res3')])
cor(Reserr[,-1])
var(Reserr[,-1])

for(i in unique(Rgvals[,1])){
  cat('#################################################################### \n')
  cat('... Generation ....',i,'\n')
  varres <- var(Rgvals[which(Rgvals[,1]==i),-1])
  print(diag(varres))
  cat('#################################################################### \n')
} 

dataall <- rbind.data.frame(dat,datoff2)
fams <- unique(dataall$Dam)[-1]

for(i in fams){
  datofffams <- dataall[which(dataall$Dam==i),c('Sire','Dam','TBV1','TBV2','Polygene')]
  datoffcor <- t(data.frame(R=0.5*rowSums(t(datofffams[,-1:-2]))))
  parents <- dataall[dataall$Progeny %in% c(datofffams$Sire,datofffams$Dam),c('TBV1','TBV2','Polygene')]
  parents <- t(data.frame(R=0.5*rowSums(t(parents))))
  cors <- cbind.data.frame(datoffcor,parents)
  if(i==fams[1]){
    corelstor <- cors
  } else {corelstor <- rbind.data.frame(corelstor,cors)}
}

cor(corelstor)
var(corelstor)

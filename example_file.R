#### requires package pedigree for inbreeding calculation
require('pedigree')

### single trait simulation (quantitative trait)
### define (co)variance structure
G <- matrix(c(va1=7.50),nrow=1,ncol=1,byrow=TRUE)
R <- matrix(c(ve1=17.50),nrow=1,ncol=1,byrow=TRUE)
mu <- c(250)

source('mmtraitsim.R')
basedata <- makebasepop(nsires=150,ndams=300,mu=mu,Va=G,Ve=R)
cor(basedata[,-1:-6])
var(basedata[,-1:-6])
offdata <- makeoff(Numgen=15,basedata,nsires=150,ndams=300,ls=20,
                    Va=G,Ve=R,sd='tbv/h',md='nested',trsel=1,selindex=c(0.55,0.45))
dataall <- rbind.data.frame(basedata,offdata)

gen <- dataall[which(dataall$G==1),]
table(gen$Dam)

########### checking simulated data for variance reudction, inbreeding, genetic progress
varGen <- aggregate.data.frame(dataall[,c('TBV1')],list(dataall$G),FUN='var')
colnames(varGen) <- c('G','TBV1')
varRes <- aggregate.data.frame(dataall[,c('Res1')],list(dataall$G),FUN='var')
colnames(varRes) <- c('G','Res1')
inbF <- aggregate.data.frame(dataall[,c('Fped')],list(dataall$G),FUN='mean')
colnames(inbF) <- c('G','F')
deltaG <- aggregate.data.frame(dataall[,c('TBV1')],list(dataall$G),FUN='mean')
colnames(deltaG) <- c('G','deltaG')

plot(y=varGen$TBV1,x=varGen$G,type = 'b',pch=20,ylim=range(varGen$TBV1)+c(-2,2),
     xlab='Generation',ylab='var Trait 1 (selected trait)')
plot(y=varRes$Res1,x=varRes$G,type = 'b',pch=20,ylim=range(varRes$Res1)+c(-2,2),col='red',
     xlab='Generation',ylab='Residual variance')
plot(y=inbF$F,x=inbF$G,type = 'b',pch=20,ylim=range(inbF$F)+c(0,0.05),
     xlab='Generation',ylab='Inbreeding')
plot(y=deltaG$deltaG,x=deltaG$G,type = 'b',pch=20,ylim=range(deltaG$deltaG)+c(-20,20),
     xlab='Generation',ylab='Genetic changes')

#### Multi-trait simulation
G <- matrix(c(va1=7.50,4.35,4.35,va2=10.00),nrow=2,ncol=2,byrow=TRUE)
R <- matrix(c(ve1=17.50,0.00,0.00,ve2=90),nrow=2,ncol=2,byrow=TRUE)
mu <- c(6,10)
source('mmtraitsim.R')
basedata <- makebasepop(nsires=150,ndams=300,mu=mu,Va=G,Ve=R)
offdata <- makeoff(Numgen=30,basedata,nsires=150,ndams=300,ls=15,
                   Va=G,Ve=R,sd='phen/h',md='rnd_ug',trsel=1,selindex=c(0.85,0.15))
dataall <- rbind.data.frame(basedata,offdata)

########### checking simulated data for variance reudction, inbreeding, genetic progress
varGen <- aggregate.data.frame(dataall[,c('TBV1','TBV2')],list(dataall$G),FUN='var')
colnames(varGen) <- c('G','TBV1','TBV2')
varRes <- aggregate.data.frame(dataall[,c('Res1','Res2')],list(dataall$G),FUN='var')
colnames(varRes) <- c('G','Res1','Res2')
inbF <- aggregate.data.frame(dataall[,c('Fped')],list(dataall$G),FUN='mean')
colnames(inbF) <- c('G','F')
deltaG <- aggregate.data.frame(dataall[,c('TBV1','TBV2')],list(dataall$G),FUN='mean')
colnames(deltaG) <- c('G','deltaG1','deltaG2')

plot(y=varGen$TBV1,x=varGen$G,type = 'b',pch=20,ylim=range(varGen[,-1])+c(-2,2),
     xlab='Generation',ylab='Genetic variance')
points(y=varGen$TBV2,x=varGen$G,type = 'b',pch=8,col='red')
plot(y=inbF$F,x=inbF$G,type = 'b',pch=20,ylim=range(inbF$F)+c(0,0.05),
     xlab='Generation',ylab='Inbreeding')
plot(y=deltaG$deltaG1,x=deltaG$G,type = 'b',pch=20,ylim=range(deltaG[,-1])+c(-20,20),
     xlab='Generation',ylab='Genetic changes')
points(y=deltaG$deltaG2,x=deltaG$G,type = 'b',pch=8,col='red')

(diff(c(0,inbF$F))/(1-inbF$F))*100

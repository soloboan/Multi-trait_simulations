G <- matrix(c(va1=7.50,4.35,4.35,va2=10.00),nrow=2,ncol=2,byrow=TRUE)
R <- matrix(c(ve1=17.50,0.00,0.00,ve2=90),nrow=2,ncol=2,byrow=TRUE)
mu <- c(6,250)

source('mmtraitsim.R')
basedata <- makebasepop(nsires=3000,ndams=3000,mu=mu,Va=G,Ve=R)
cor(basedata[,-1:-5])
var(basedata[,-1:-5])
offdata <- makeoff(Numgen=20,basedata,nsires=150,ndams=300,ls=15,
                    Va=G,Ve=R,sd='index/h',md='nested',trsel=1,selindex=c(0.55,0.45))

dataall <- rbind.data.frame(basedata,offdata)
varGen <- aggregate.data.frame(dataall[,c('TBV1','TBV2')],
                               list(dataall$G),FUN='var')
varRes <- aggregate.data.frame(dataall[,c('Res1','Res2')],
                               list(dataall$G),FUN='var')

plot(y=varGen$TBV1,x=varGen$Group.1,type = 'b',pch=20,
     ylim=range(varGen$TBV1)+c(-2,2),
     xlab='Generation',ylab='var Trait 1 (selected trait)')
plot(y=varGen$TBV2,x=varGen$Group.1,type = 'b',pch=20,
     ylim=range(varGen$TBV2)+c(-2,2),col='red',
     xlab='Generation',ylab='var Trait 2')

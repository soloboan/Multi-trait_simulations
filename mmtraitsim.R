#######################################################################
###### Making a base population for subsequent generations
#######################################################################
require(MASS)
require(pedigree)

makebasepop <- function(nsires,ndams,mu,Va,Ve){
  ID <- 1:sum(nsires,ndams)
  nanims <- sum(nsires,ndams)
  TBV <- data.frame(round(mvrnorm(nanims,mu,Va),6))
  TBV <- scale(TBV)
  sdtbv <- sqrt(diag(Va))
  for(s in 1:ncol(Va)){TBV[,s] <- mu[s] + (TBV[,s]*sdtbv[s])}
  E <- data.frame(round(mvrnorm(nanims,rep(0,length(mu)),Ve),6))
  E <- scale(E)
  sde <- sqrt(diag(Ve))
  for(s in 1:ncol(Ve)){E[,s] <- 0 + (E[,s]*sde[s])}
  pheno <- TBV + E
  datafile <- data.frame(TBV,E,pheno)
  colnames(datafile) <- c(paste('TBV',1:nrow(Va),sep=''),paste('Res',1:nrow(Va),sep=''),paste('Phen',1:nrow(Va),sep=''))
  Sex <- c(rep('M',nsires),rep('F',ndams))
  
  #### calculate inbreeding using pedigree package
  Fped <- calcInbreeding(data.frame(ID,sire=0,dam=0))
  basedata <- data.frame(G=0,ID,Sire=0,Dam=0,Sex,Fped,datafile)
  return(basedata) 
}

#######################################################################
###### Making offspring population from the base population
#######################################################################

makeoff <- function(Numgen,basedata,nsires,ndams,ls,Va,Ve,sd,md,trsel,selindex){
  for (m in 1:Numgen){
    if(m>1){basedata <- offspring}
    sires <- basedata[which(basedata$Sex=='M'),]
    dams <- basedata[which(basedata$Sex=='F'),]
    noff <- ndams*ls
    
    ############### selection design for parents ##################    
    if(tolower(sd)=='rnd'){
      s <- sample(x=sires$ID,size=nsires,replace=F)
      d <- sample(x=dams$ID,size=ndams,replace=F)
    } else if(tolower(sd)=='phen/h'){
      s <- sires[order(sires[,paste('Phen',trsel,sep='')],decreasing=T),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,paste('Phen',trsel,sep='')],decreasing=T),'ID']
      d <- d[1:ndams]
      s <- sample(x=s,size=nsires,replace=F)
      d <- sample(x=s,size=ndams,replace=F)
    } else if(tolower(sd)=='phen/l'){
      s <- sires[order(sires[,paste('Phen',trsel,sep='')],decreasing=F),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,paste('Phen',trsel,sep='')],decreasing=F),'ID']
      d <- d[1:ndams]
      s <- sample(x=s,size=nsires,replace=F)
      d <- sample(x=s,size=ndams,replace=F)
    } else if(tolower(sd)=='tbv/h'){
      s <- sires[order(sires[,paste('TBV',trsel,sep='')],decreasing=T),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,paste('TBV',trsel,sep='')],decreasing=T),'ID']
      d <- d[1:ndams]
      s <- sample(x=s,size=nsires,replace=F)
      d <- sample(x=s,size=ndams,replace=F)
    } else if(tolower(sd)=='tbv/l'){
      s <- sires[order(sires[,paste('TBV',trsel,sep='')],decreasing=F),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,paste('TBV',trsel,sep='')],decreasing=F),'ID']
      d <- d[1:ndams]
      s <- sample(x=s,size=nsires,replace=F)
      d <- sample(x=s,size=ndams,replace=F)
    } else if(tolower(sd)=='index/h'){
      indexW <- matrix(selindex/sqrt(diag(Va)),nrow=nrow(Va),ncol=1)
      Sireindex <- as.matrix(sires[,paste('TBV',1:nrow(Va),sep='')])
      sires$index <- Sireindex  %*% indexW
      Damindex <- as.matrix(dams[,paste('TBV',1:nrow(Va),sep='')])
      dams$index <- Damindex  %*% indexW
      s <- sires[order(sires[,c('index')],decreasing=T),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,c('index')],decreasing=T),'ID']
      d <- d[1:ndams]
      s <- sample(x=s,size=nsires,replace=F)
      d <- sample(x=s,size=ndams,replace=F)
    } else if(tolower(sd)=='index/l'){
      indexW <- matrix(selindex/sqrt(diag(Va)),nrow=nrow(Va),ncol=1)
      Sireindex <- as.matrix(sires[,paste('TBV',1:nrow(Va),sep='')])
      sires$index <- Sireindex  %*% indexW
      Damindex <- as.matrix(dams[,paste('TBV',1:nrow(Va),sep='')])
      dams$index <- Damindex  %*% indexW
      s <- sires[order(sires[,c('index')],decreasing=T),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,c('index')],decreasing=T),'ID']
      d <- d[1:ndams]
      s <- sample(x=s,size=nsires,replace=F)
      d <- sample(x=s,size=ndams,replace=F)
    } else if(tolower(sd)=='phenindex/h'){
      indexW <- matrix(selindex/sqrt(diag(Va)+diag(Ve)),nrow=nrow(Va),ncol=1)
      Sireindex <- as.matrix(sires[,paste('Phen',1:nrow(Va),sep='')])
      sires$index <- Sireindex  %*% indexW
      Damindex <- as.matrix(dams[,paste('Phen',1:nrow(Va),sep='')])
      dams$index <- Damindex  %*% indexW
      s <- sires[order(sires[,c('index')],decreasing=T),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,c('index')],decreasing=T),'ID']
      d <- d[1:ndams]
      s <- sample(x=s,size=nsires,replace=F)
      d <- sample(x=s,size=ndams,replace=F)
    } else if(tolower(sd)=='phenindex/l'){
      indexW <- matrix(selindex/sqrt(diag(Va)+diag(Ve)),nrow=nrow(Va),ncol=1)
      Sireindex <- as.matrix(sires[,paste('Phen',1:nrow(Va),sep='')])
      sires$index <- Sireindex  %*% indexW
      Damindex <- as.matrix(dams[,paste('Phen',1:nrow(Va),sep='')])
      dams$index <- Damindex  %*% indexW
      s <- sires[order(sires[,c('index')],decreasing=T),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,c('index')],decreasing=T),'ID']
      d <- d[1:ndams]
      s <- sample(x=s,size=nsires,replace=F)
      d <- sample(x=s,size=ndams,replace=F)
    }
    
    ################## mating design  ############
    if(tolower(md)=='rnd_ug'){
      use.sires <- rep(s,length.out=noff)
      use.dams <- sample(x=rep(d,length.out=noff),size=noff,replace=F)
    } else if(tolower(md)=='nested'){
      use.sires <- rep(s,length.out=noff)
      use.dams <- sample(x=rep(d,length.out=noff),size=noff,replace=F)
    } 
    #  else if(substr(md,1,3)=='fac'){
    #   partialfacnumber <- as.numeric(gsub(x=unlist(strsplit(md,split='\\['))[2],pattern='\\]',replacement=''))
    #   use.sires <- rep(rep(sort(s),each=ls),length.out=noff*partialfacnumber)
    #   use.dams <- sort(rep(d,each=ls*partialfacnumber))
    # }
    
    ################# making pedigree  ##################
    parent <- cbind.data.frame(Sire=use.sires,Dam=use.dams)
    parent <- parent[order(parent$Sire,parent$Dam),]
    ID <- (1:noff)+tail(basedata,1)[,2]
    offspring <- cbind.data.frame(ID,parent)
    
    ######## compute average inbreeding ##########
    avF <- mean(basedata[basedata$ID %in% c(use.sires,use.dams),'Fped'])
    
    ######## sampling MS and Residuals ##########
    MS <- data.frame(round(mvrnorm(noff,rep(0,length(mu)),0.5*Va*(1-avF)),6))
    MS <- scale(MS)
    sdms <- sqrt(0.5*diag(Va)*(1-avF))
    for(j in 1:ncol(Va)){MS[,j] <- 0 + (MS[,j]*sdms[j])}
    colnames(MS) <- paste('MS',1:nrow(Va),sep='')
    E <- data.frame(round(mvrnorm(noff,rep(0,length(mu)),Ve),6))
    E <- scale(E)
    sde <- sqrt(diag(Ve))
    for(j in 1:ncol(Ve)){E[,j] <- 0 + (E[,j]*sde[j])}
    colnames(E) <- paste('Res',1:nrow(Ve),sep='')
    
    ########### computing BV (Parents average) + MS  ######################  
    ebvsire <- merge(offspring[,c('ID','Sire')],basedata[,-c(1,3,4,5)],by.x='Sire',by.y='ID')[,c('ID',paste('TBV',1:nrow(Va),sep=''))]
    ebvdam <- merge(offspring[,c('ID','Dam')],basedata[,-c(1,3,4,5)],by.x='Dam',by.y='ID')[,c('ID',paste('TBV',1:nrow(Va),sep=''))]
    ebvparents <- merge(ebvsire,ebvdam,by='ID')[,-1]
    sirecol=c(1:ncol(Ve))
    damcol=c(sirecol+ncol(Ve))
    ebvoff <- 0.5*ebvparents[,sirecol] + 0.5*ebvparents[,damcol] + MS
    colnames(ebvoff) <- paste('TBV',1:nrow(Va),sep='')
    
    ############## making phenotypes  ##################
    pheno <- ebvoff + E
    colnames(pheno) <- paste('Phen',1:nrow(Va),sep='')
    
    ################ assign sex #########################
    Sex <- sample(rep(c('M','F'),noff/2),size=noff,replace=F)
    
    ######## final datafile containg all columns ############
    offspring <- cbind.data.frame(offspring,Sex)
    offspring <- cbind.data.frame(G=m,offspring,Fped=0,ebvoff,E,pheno)
    if(m==1){
      offspringgen <- offspring
      #### calculate inbreeding using pedigree package
      offspringgen$Fped <- round(calcInbreeding(data.frame(offspringgen[,c('ID','Sire','Dam')])),4)
    } else {
      offspringgen <- rbind.data.frame(offspringgen,offspring)
      #### calculate inbreeding using pedigree package
      offspringgen$Fped <- round(calcInbreeding(data.frame(offspringgen[,c('ID','Sire','Dam')])),4)
    }
    cat('... generation ...',m,' ... completed ...\n')
  }
  return(offspringgen)
} 

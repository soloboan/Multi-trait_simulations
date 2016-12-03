#######################################################################
###### To get base population
#######################################################################
library(MASS)
makebasepop <- function(nsires=100,ndams=100,mu=c(0.0025,100,1000),Va,Ve){
  ID <- 1:sum(nsires,ndams)
  nanims <- sum(nsires,ndams)
  TBV <- data.frame(round(mvrnorm(nanims,mu,Va),6))
  TBV <- scale(TBV)
  sdtbv <- sqrt(diag(Va))
  for(s in 1:ncol(Va)){
    TBV[,s] <- mu[s] + (TBV[,s]*sdtbv[s])
  }
  E <- data.frame(round(mvrnorm(nanims,rep(0,length(mu)),Ve),6))
  E <- scale(E)
  sde <- sqrt(diag(Ve))
  for(s in 1:ncol(Ve)){
    E[,s] <- 0 + (E[,s]*sde[s])
  }
  pheno <- TBV + E
  datafile <- data.frame(TBV,E,pheno)
  colnames(datafile) <- c(paste('TBV',1:nrow(Va),sep=''),paste('Res',1:nrow(Va),sep=''),paste('Phen',1:nrow(Va),sep=''))
  Sex <- c(rep('M',nsires),rep('F',ndams))
  basedata <- data.frame(G=0,ID,Sire=0,Dam=0,Sex,datafile)
  return(basedata) 
}

makeoff <- function(Numgen,basedata,nsires=50,ndams=1000,ls=5,Va,Ve,sd='tbv/h',md='rnd_ug',trsel=2) {
  for (m in 1:Numgen){
    if(m>1){basedata <- offspring}
    sires <- basedata[which(basedata$Sex=='M'),]
    dams <- basedata[which(basedata$Sex=='F'),]
    noff <- ndams*ls
    
    if(sd=='rnd'){
      s <- sort(sample(x=sires$ID,size=nsires,replace=F))
      d <- sort(sample(x=dams$ID,size=ndams,replace=F))
      #d <- sample(x=rep(dams$ID,ls),size=ndams*ls,replace=F)
    } else if(sd=='tbv/h'){
      s <- sires[order(sires[,trsel],decreasing=T),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,trsel],decreasing=T),'ID']
      d <- d[1:ndams]
    } else if(sd=='tbv/l'){
      s <- sires[order(sires[,trsel],decreasing=F),'ID']
      s <- s[1:nsires]
      d <- dams[order(dams[,trsel],decreasing=F),'ID']
      d <- d[1:ndams]
    }
    ################## mating design  ############
    if(md=='rnd_ug'){
      use.sires <- sort(rep(s,length.out=noff))
      use.dams <- sample(x=rep(d,length.out=noff),size=noff,replace=F)
    } else if(md=='nested'){
      use.sires <- sort(rep(s,length.out=noff))
      use.dams <- sort(sample(x=rep(d,length.out=noff),size=noff,replace=F))
    }
    
    ################# making pedigree  ##################
    parent <- cbind.data.frame(Sire=use.sires,Dam=use.dams)
    parent <- parent[order(parent$Sire,parent$Dam),]
    ID <- (1:noff)+tail(basedata,1)[,2]
    offspring <- cbind.data.frame(ID,parent)
    
    ######## sampleing MS and Residuals ##########
    MS <- data.frame(round(mvrnorm(noff,rep(0,length(mu)),0.5*Va),6))
    MS <- scale(MS)
    sdms <- sqrt(0.5*diag(Va))
    for(j in 1:ncol(Va)){MS[,j] <- 0 + (MS[,j]*sdms[j])}
    colnames(MS) <- paste('MS',1:nrow(Va),sep='')
    E <- data.frame(round(mvrnorm(noff,rep(0,length(mu)),Ve),6))
    E <- scale(E)
    sde <- sqrt(diag(Ve))
    for(j in 1:ncol(Ve)){E[,j] <- 0 + (E[,j]*sde[j])}
    colnames(E) <- paste('Res',1:nrow(Ve),sep='')
    
    ########### computing BV (Parents average) + MS 
    ebvsire <- merge(basedata,data.frame(ID=use.sires),by='ID',sort=F)[,paste('TBV',1:nrow(Va),sep='')]
    ebvdam <- merge(basedata,data.frame(ID=use.dams),by='ID',sort=F)[,paste('TBV',1:nrow(Va),sep='')]
    ebvoff <- 0.5*ebvsire + 0.5*ebvdam + MS
    colnames(ebvoff) <- paste('TBV',1:nrow(Va),sep='')
    
    ############## making phenotypes  ##################
    pheno <- ebvoff + E
    colnames(pheno) <- paste('Phen',1:nrow(Va),sep='')
    
    ################ assign sex #########################
    if(ls==1){
      Sex=sample(c('M','F'),size=ndams,replace=T)
    } else {
      fams <- unique(offspring$Dam)
      fams <- fams[fams!=0]
      for(l in fams){
        offfams <- data.frame(ID=offspring[which(offspring$Dam==l),'ID'])
        offfams$Sex <- sort(rep(c('M','F'),length.out=nrow(offfams)),decreasing=T)
        if(l==fams[1]){offSex <- offfams} else {offSex <- rbind.data.frame(offSex,offfams)}
      }
    }
    ######## final datafile containg all columns ############
    offspring <- merge(offspring,offSex,by='ID')
    offspring <- cbind.data.frame(G=m,offspring,ebvoff,E,pheno)
    if(m==1){offspringgen <- offspring} else {offspringgen <- rbind.data.frame(offspringgen,offspring)}
    cat('... generation ...',m,' ... completed ...\n')
  }
  return(offspringgen)
} 


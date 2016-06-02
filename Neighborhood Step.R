## Put Data Back Into original Image
COUNT <- 0
IMAGE <- matrix(NA,512,512)
for (ii in 1:512) {
  for (jj in 1:512) {
    if (MASK[ii,jj]==1) {
      COUNT <- COUNT + 1
      IMAGE[ii,jj] <- which.max(TAU[COUNT,])
    }
  }
}
par( mar=c(2,2,2,2))
image.plot(1:512,1:512,IMAGE,col=tim.colors(17),xlab='',ylab='',main='Stage 1')

## Set d Neighborhood Structure
DIST <- 2

NEIGHBOURS <- list()
COUNT <- 0
for (ii in 1:512) {
  for (jj in 1:512) {
    if (MASK[ii,jj]==1) {
      COUNT <- COUNT + 1
      NEIGHBOURS[[COUNT]] <- c(IMAGE[max((ii-DIST),1):min((ii+DIST),512),max((jj-DIST),1):min((jj+DIST),512)])
    }
  }
}

## Get Eta values
ETA <- matrix(NA,NN,GG)
for (nn in 1:NN) {
  for (gg in 1:GG) {
    ETA[nn,gg] <- mean(NEIGHBOURS[[nn]]==gg,na.rm=T) 
  }
}
## Get data
CC <- matrix(NA,NN,GG)
for (nn in 1:NN) {
  for (gg in 1:GG) {
    CC[nn,gg] <- as.numeric(which.max(TAU[nn,])==gg)
  }
}

## Make Parameters
PARA <- matrix(0,2,GG)
for (gg in 1:(GG-1)) {
  PARA[,gg] <- runif(2)
}

### Load C Code
Rcpp::sourceCpp('Pseudolike.cpp')

PL <- PSEUDO_FUN(CC,ETA,PARA,NN,GG)

OptiFUN <- function(LP) {
  PARA1 <- matrix(c(LP,0,0),2,GG)
  -PSEUDO_FUN(CC,ETA,PARA1,NN,GG)
}

OPTIM <- optim(rep(0,2*GG-2),OptiFUN,method='BFGS')
PARAOUT <- matrix(c(OPTIM$par,0,0),2,GG)
-2*(-OPTIM$value) + log(NN)*(2*GG-2)

PROBS <- matrix(NA,NN,GG)
for (nn in 1:NN) {
  for (gg in 1:GG) {
    PROBS[nn,gg] <- exp(sum(c(1,ETA[nn,gg])*PARAOUT[,gg]))    
  }
}

CLUSTER <- c()
for (nn in 1:NN) {
  CLUSTER[nn] <- which.max(PROBS[nn,])
}

## Put CLUSTER Back Into original Image
COUNT <- 0
SMOOTH <- matrix(NA,512,512)
for (ii in 1:512) {
  for (jj in 1:512) {
    if (MASK[ii,jj]==1) {
      COUNT <- COUNT + 1
      SMOOTH[ii,jj] <- CLUSTER[COUNT]
    }
  }
}
par( mar=c(2,2,2,2))
image.plot(1:512,1:512,SMOOTH,col=tim.colors(17),xlab='',ylab='',main='Stage 2')

plot(YY[[200]],type='l')

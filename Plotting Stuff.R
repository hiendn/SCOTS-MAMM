COUNT <- 0
for (ii in 1:512) {
  for (jj in 1:512) {
    if (MASK[ii,jj]==1) {
      COUNT <- COUNT + 1
      while (min(YY[[COUNT]])< -5.5) {
        YY[[COUNT]][which(YY[[COUNT]]< -5.5)] <- YY[[COUNT]][sample(which(YY[[COUNT]] >= -5.5),length(which(YY[[COUNT]]< -5.5)))] 
      }
 
    }
  }
}

CLUSTERING <- apply(TAU,1,which.max)
par(mfrow=c(3,6))

for (gg in 1:17)
{
  WHICH <- which(CLUSTERING==gg)
  HOLDER <- matrix(NA,length(WHICH),481)
  for (ii in 1:length(WHICH)) {
    HOLDER[ii,] <- YY[[WHICH[ii]]][10:490]
  }
  MED <- apply(HOLDER,2,quantile,0.5)
  UPPER <- apply(HOLDER,2,quantile,0.975)
  LOWER <- apply(HOLDER,2,quantile,0.025)
  UNLIST <- unlist(HOLDER)
  plot(seq(min(UNLIST),max(UNLIST),length.out=500),main=paste('Cluster',gg),ylab='Y',xlab='',type='n')
  lines(10:490,MED,col='black')
  lines(10:490,UPPER,col='blue')
  lines(10:490,LOWER,col='blue')
}

library(fields)
TABLE_PLOT <- matrix(NA,25,20)
for (ii in 1:25) {
  for (jj in 1:20) {
    if (ii + jj <= 26) {
      TABLE_PLOT[ii,jj] <- TABLE[ii,jj]
    }
  }
}
image.plot(1:25,1:20,TABLE_PLOT,xlab='g',ylab='p',col=tim.colors(10000))
points(17,3,pch=20,col='white',cex=2)

## Put CLUSTER Back Into original Image
COUNT <- 0
SMOOTH <- matrix(NA,512,512)
for (ii in 1:512) {
  for (jj in 1:512) {
    if (MASK[ii,jj]==1) {
      COUNT <- COUNT + 1
      SMOOTH[ii,jj] <- YMEAN[COUNT]
    }
  }
}
image.plot(SMOOTH)

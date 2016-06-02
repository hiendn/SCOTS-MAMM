### Simulation
SIM1 <- matrix(NA,100,100)
SIM2 <- matrix(NA,100,100)
SIM3 <- matrix(NA,100,100)
SIM4 <- matrix(NA,100,100)
for (ii in 1:100) {
  for (jj in 1:100) {
    if (floor((ii-1/2)/25) %in% c(0,2)) {
      if (floor((jj-1/2)/25) %in% c(0,2)) {
        SIM1[ii,jj] <- 1
      } else {
        SIM1[ii,jj] <- 2
      }
    } else {
      if (floor((jj-1/2)/25) %in% c(1,3)) {
        SIM1[ii,jj] <- 1
      } else {
        SIM1[ii,jj] <- 2
      }
    }
  }
}

for (ii in 1:100) {
  for (jj in 1:100) {
    if (floor((ii-1/2)/20) %in% c(0,2,4)) {
      if (floor((jj-1/2)/20) %in% c(0,2,4)) {
        SIM2[ii,jj] <- 1
      } else {
        SIM2[ii,jj] <- 2
      }
    } else {
      if (floor((jj-1/2)/20) %in% c(1,3)) {
        SIM2[ii,jj] <- 1
      } else {
        SIM2[ii,jj] <- 2
      }
    }
  }
}

for (ii in 1:100) {
  for (jj in 1:100) {
    if (floor((ii-1/2)/25) %in% c(0,2)) {
      if (floor((jj-1/2)/25) %in% c(0,2)) {
        SIM3[ii,jj] <- 1
      } else {
        SIM3[ii,jj] <- 2
      }
    } else {
      if (floor((jj-1/2)/25) %in% c(1,3)) {
        SIM3[ii,jj] <- 3
      } else {
        SIM3[ii,jj] <- 4
      }
    }
  }
}

for (ii in 1:100) {
  for (jj in 1:100) {
    if (floor((ii-1/2)/20) %in% c(0,2,4)) {
      if (floor((jj-1/2)/20) %in% c(0,2,4)) {
        SIM4[ii,jj] <- 1
      } else {
        SIM4[ii,jj] <- 2
      }
    } else {
      if (floor((jj-1/2)/20) %in% c(1,3)) {
        SIM4[ii,jj] <- 3
      } else {
        SIM4[ii,jj] <- 4
      }
    }
  }
}

par(mfrow=c(2,2))
image(1:100,1:100,SIM1,col=tim.colors(4),xlab='',ylab='',main='S1')
image(1:100,1:100,SIM2,col=tim.colors(4),xlab='',ylab='',main='S2')
image(1:100,1:100,SIM3,col=tim.colors(4),xlab='',ylab='',main='S3')
image(1:100,1:100,SIM4,col=tim.colors(4),xlab='',ylab='',main='S4')

CASES <- list()
CASES[[1]] <- c(0,0.25)
CASES[[2]] <- c(0,-0.25)
CASES[[3]] <- c(0.25,0)
CASES[[4]] <- c(-0.25,0)

SS1 <- array(NA,c(100,100,100))
SS2 <- array(NA,c(100,100,100))
for (ii in 1:100) {
  for (jj in 1:100) {
    SS1[ii,jj,] <- arima.sim(n=100,list(ar = CASES[[SIM1[ii,jj]]]))
    SS2[ii,jj,] <- arima.sim(n=100,list(ar = CASES[[SIM2[ii,jj]]]))
  }
}

par(mfrow=c(4,1))
plot(seq(-3.5,3.5,length=100),type='n',main='C1',ylab='Y')
lines(arima.sim(n=100,list(ar = CASES[[1]])),lty=1,col=c('blue','red','black')[1])
lines(arima.sim(n=100,list(ar = CASES[[1]])),lty=2,col=c('blue','red','black')[2])
lines(arima.sim(n=100,list(ar = CASES[[1]])),lty=3,col=c('blue','red','black')[3])
plot(seq(-3.5,3.5,length=100),type='n',main='C2',ylab='Y')
lines(arima.sim(n=100,list(ar = CASES[[2]])),lty=1,col=c('blue','red','black')[1])
lines(arima.sim(n=100,list(ar = CASES[[2]])),lty=2,col=c('blue','red','black')[2])
lines(arima.sim(n=100,list(ar = CASES[[2]])),lty=3,col=c('blue','red','black')[3])
plot(seq(-3.5,3.5,length=100),type='n',main='C3',ylab='Y')
lines(arima.sim(n=100,list(ar = CASES[[3]])),lty=1,col=c('blue','red','black')[1])
lines(arima.sim(n=100,list(ar = CASES[[3]])),lty=2,col=c('blue','red','black')[2])
lines(arima.sim(n=100,list(ar = CASES[[3]])),lty=3,col=c('blue','red','black')[3])
plot(seq(-3.5,3.5,length=100),type='n',main='C4',ylab='Y')
lines(arima.sim(n=100,list(ar = CASES[[4]])),lty=1,col=c('blue','red','black')[1])
lines(arima.sim(n=100,list(ar = CASES[[4]])),lty=2,col=c('blue','red','black')[2])
lines(arima.sim(n=100,list(ar = CASES[[4]])),lty=3,col=c('blue','red','black')[3])

par(mfrow=c(2,2))
image(1:100,1:100,IMAGE1,col=tim.colors(4),xlab='',ylab='',main='S1')
image(1:100,1:100,IMAGE2,col=tim.colors(4),xlab='',ylab='',main='S2')
image(1:100,1:100,IMAGE3,col=tim.colors(4),xlab='',ylab='',main='S3')
image(1:100,1:100,IMAGE4,col=tim.colors(4),xlab='',ylab='',main='S4')

par(mfrow=c(2,2))
image(1:100,1:100,SMOOTH1,col=tim.colors(4),xlab='',ylab='',main='S1')
image(1:100,1:100,SMOOTH2,col=tim.colors(4),xlab='',ylab='',main='S2')
image(1:100,1:100,SMOOTH3,col=tim.colors(4),xlab='',ylab='',main='S3')
image(1:100,1:100,SMOOTH4,col=tim.colors(4),xlab='',ylab='',main='S4')

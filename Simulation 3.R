library(mclust)

### Simulation
SIM1 <- matrix(NA,100,100)
SIM2 <- matrix(NA,100,100)
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
        SIM1[ii,jj] <- 3
      } else {
        SIM1[ii,jj] <- 4
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
        SIM2[ii,jj] <- 3
      } else {
        SIM2[ii,jj] <- 4
      }
    }
  }
}

CASES <- list()
CASES[[1]] <- c(0,0.25)
CASES[[2]] <- c(0,-0.25)
CASES[[3]] <- c(0.25,0)
CASES[[4]] <- c(-0.25,0)


TABLE <- matrix(NA,100,3)







for (scc in 1:100) {
  
  SSS <- SIM2
  
  SS1 <- array(NA,c(100,100,100))
  #SS2 <- array(NA,c(100,100,100))
  for (ii in 1:100) {
    for (jj in 1:100) {
      SS1[ii,jj,] <- arima.sim(n=100,list(ar = CASES[[SSS[ii,jj]]]))
      #SS2[ii,jj,] <- arima.sim(n=100,list(ar = CASES[[SIM2[ii,jj]]]))
    }
  }
  
  for (gpos in 4:4) {
    for (ppos in 2:2) {
      
      
      
      DFRAME <- SS1
      MASK <- matrix(1,100,100)
      
      GG <- gpos
      PP <- ppos
      NN <- sum(MASK)
      MM <- 100
      
      YY <- list()
      COUNT <- 0
      for (ii in 1:100) {
        for (jj in 1:100) {
          if (MASK[ii,jj]==1) {
            COUNT <- COUNT + 1
            YY[[COUNT]] <- DFRAME[ii,jj,]      
          }
        }
      }
      
      # COUNT <- 0
      # for (ii in 1:512) {
      #   for (jj in 1:512) {
      #     if (MASK[ii,jj]==1) {
      #       COUNT <- COUNT + 1
      #       while (min(YY[[COUNT]])< -5.1) {
      #         YY[[COUNT]][which(YY[[COUNT]]< -5.1)] <- YY[[COUNT]][sample(which(YY[[COUNT]] >= -5.1),length(which(YY[[COUNT]]< -5.1)))] 
      #       }
      #  
      #     }
      #   }
      # }
      
      rm(list=c('DFRAME'))
      
      # COUNT <- 0
      # for (ii in 1:512) {
      #   for (jj in 1:512) {
      #     COUNT <- COUNT + 1
      #     YY[[COUNT]][which(YY[[COUNT]]==-Inf)] <- sample(YY[[COUNT]][which(YY[[COUNT]]!=-Inf)],sum(YY[[COUNT]]==-Inf),replace=T)
      #   }
      # }
      XX <- list()
      COUNT <- 0
      for (ii in 1:100) {
        for (jj in 1:100) {
          if (MASK[ii,jj]==1) {
            COUNT <- COUNT + 1
            XX[[COUNT]] <- rep(1,MM-PP)
            for (pp in PP:1) {
              XX[[COUNT]] <- cbind(XX[[COUNT]],YY[[COUNT]][pp:(MM-(PP-pp+1))])
            }   
          }
          
        }
      }
      for (ss in 1:NN) {
        COUNT <- ss
        YY[[COUNT]] <- YY[[COUNT]][(PP+1):MM]
      }
      
      YMEAN <- c()
      for (nn in 1:NN) {
        YMEAN[nn] <- mean(YY[[nn]])
      }
      QMEAN <- quantile(YMEAN,(1:GG-1/2)/GG)
      
      PI <- rep(1/GG,GG)
      BETA <- matrix(mean(YY[[1]]),PP+1,GG)
      for (gg in 1:GG) {
        BETA[,gg] <- c(QMEAN[gg],rep(0,PP))
      }
      VAR <- rep(var(YY[[1]]),GG)
      for (gg in 1:GG) {
        VAR[gg] <- var(YY[[gg]])
      }
      
      Rcpp::sourceCpp('Loglike.cpp')
      Rcpp::sourceCpp('TauCalc.cpp')
      Rcpp::sourceCpp('Var_Calc.cpp')
      Rcpp::sourceCpp('BetaCalc.cpp')
      
      LL <- LOGLIKE_FUN(XX,YY,BETA,PI,VAR,MM,NN,GG,PP)
      LL
      
      TAU <- matrix(1/GG,NN,GG)
      TAU_OLD <- TAU
      BETA_OLD <- BETA
      VAR_OLD <- VAR
      
      for (ii in 1:10) {
        
        TAU <- TAU_FUN(TAU,XX,YY,BETA,PI,VAR,MM,NN,GG,PP)
        for (ss in 1:NN) {
          if ( is.na(sum(TAU[ss,])) ) {
            TAU[ss,] <- PI
          }
        }
        
        for (gg in 1:GG) {
          PI[gg] <- sum(TAU[,gg])/NN
        }
        
        BETA <- BETA_FUN(TAU,XX,YY,BETA,PI,VAR,MM,NN,GG,PP)
        
        VAR <- VAR_FUN(TAU,XX,YY,BETA,PI,VAR,MM,NN,GG,PP)
        
        #   print(c(ii,PI))
        #   print(BETA)
        #   print(VAR)
        
        LL <- LOGLIKE_FUN(XX,YY,BETA,PI,VAR,MM,NN,GG,PP)
        #   while (LL %in% c(-Inf)) {
        #     VAR <- VAR*2
        #     LL <- LOGLIKE_FUN(XX,YY,BETA,PI,VAR,MM,NN,GG,PP)
        #     print(LL)
        #   }
        print(c(ii,LL))
      }
      
      COUNT <- 0
      IMAGE <- matrix(NA,100,100)
      for (ii in 1:100) {
        for (jj in 1:100) {
          if (MASK[ii,jj]==1) {
            COUNT <- COUNT + 1
            IMAGE[ii,jj] <- which.max(TAU[COUNT,])
          }
        }
      }
      image(IMAGE,col=rainbow(5))      
      

      
      
      
      
      
      
      TABLE[scc,1] <- adjustedRandIndex(c(IMAGE),c(SSS))
      
      PLIC <- 100000000
      PLOLD <- Inf
      DIST <- 0
      while (PLIC < PLOLD) {
        ## Set d Neighborhood Structure
        DIST <- DIST + 1
        PLOLD <- PLIC
        
        NEIGHBOURS <- list()
        COUNT <- 0
        for (ii in 1:100) {
          for (jj in 1:100) {
            if (MASK[ii,jj]==1) {
              COUNT <- COUNT + 1
              NEIGHBOURS[[COUNT]] <- c(IMAGE[max((ii-DIST),1):min((ii+DIST),100),max((jj-DIST),1):min((jj+DIST),100)])
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
        SMOOTH <- matrix(NA,100,100)
        for (ii in 1:100) {
          for (jj in 1:100) {
            if (MASK[ii,jj]==1) {
              COUNT <- COUNT + 1
              SMOOTH[ii,jj] <- CLUSTER[COUNT]
            }
          }
        }
        image(SMOOTH,col=rainbow(6))
        PLIC <- 2*OPTIM$value
        print(c(DIST,PLIC))
      }
      
      TABLE[scc,2] <- adjustedRandIndex(c(SMOOTH),c(SSS))
      TABLE[scc,3] <- DIST
      
      save(TABLE,file='TT4.rdata')
      print(c(scc))
      print(TABLE[scc,])
    }
  }
}

      

## Copyright Hien Duy Nguyen - University of Queensland 2016/06/02

library(AnalyzeFMRI)
rm(list=ls())

DATA <- f.read.nifti.volume('T50-sml-res.nii')
MASK <- f.read.nifti.volume('T50-max-z-mask.nii')
MASK <- MASK[,,1,1]
DFRAME <- log(DATA[,,seq(1,2000,by=4),1])

rm(list=c('DATA'))

GG <- 17
PP <- 3
NN <- sum(MASK)
MM <- 500

YY <- list()
COUNT <- 0
for (ii in 1:512) {
  for (jj in 1:512) {
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
#       while (min(YY[[COUNT]])< 0.006) {
#         YY[[COUNT]][which(YY[[COUNT]]< 0.006)] <- YY[[COUNT]][sample(which(YY[[COUNT]] >= 0.006),length(which(YY[[COUNT]]< 0.006)))] 
#       }
#  
#     }
#   }
# }

# COUNT <- 0
# for (ii in 1:512) {
#   for (jj in 1:512) {
#     COUNT <- COUNT + 1
#     YY[[COUNT]][which(YY[[COUNT]]==-Inf)] <- sample(YY[[COUNT]][which(YY[[COUNT]]!=-Inf)],sum(YY[[COUNT]]==-Inf),replace=T)
#   }
# }
XX <- list()
COUNT <- 0
for (ii in 1:512) {
  for (jj in 1:512) {
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
  VAR[gg] <- var(YY[[1]])*10
}

Rcpp::sourceCpp('Loglike.cpp')
Rcpp::sourceCpp('TauCalc.cpp')
Rcpp::sourceCpp('Var_Calc.cpp')
Rcpp::sourceCpp('BetaCalc.cpp')

LL <- LOGLIKE_FUN(XX,YY,BETA,PI,VAR,MM,NN,GG,PP,20)
LL

TAU <- matrix(1/GG,NN,GG)
TAU_OLD <- TAU
BETA_OLD <- BETA
VAR_OLD <- VAR

for (ii in 1:10) {
  
  TAU <- TAU_FUN(TAU,XX,YY,BETA,PI,VAR,MM,NN,GG,PP,20)
  for (ss in 1:NN) {
    if ( is.na(sum(TAU[ss,])) ) {
      TAU[ss,] <- PI
    }
  }
  
  for (gg in 1:GG) {
    PI[gg] <- sum(TAU[,gg])/NN
  }
  
  BETA <- BETA_FUN(TAU,XX,YY,BETA,PI,VAR,MM,NN,GG,PP,20)
  
  VAR <- VAR_FUN(TAU,XX,YY,BETA,PI,VAR,MM,NN,GG,PP,20)
  
#   print(c(ii,PI))
#   print(BETA)
#   print(VAR)
  
  LL <- LOGLIKE_FUN(XX,YY,BETA,PI,VAR,MM,NN,GG,PP,20)
#   while (LL %in% c(-Inf)) {
#     VAR <- VAR*2
#     LL <- LOGLIKE_FUN(XX,YY,BETA,PI,VAR,MM,NN,GG,PP)
#     print(LL)
#   }
  print(c(ii,LL))
}

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
image(IMAGE,col=rainbow(17))

-2*LL + (GG*(PP+3)-1)*log(NN)

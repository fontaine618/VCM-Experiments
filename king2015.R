library(mlbench)
library(grplasso)
data(BostonHousing2)

# dp function
constant<-3.7
dp<-function(lambda,theta)
{
  dp<-lambda*((theta<=lambda)+(constant*lambda-theta)*(constant*lambda>theta)/(constant-1)/lambda*(theta>lambda))
}

originalx<-sqrt(BostonHousing2[,19])

originalync<-BostonHousing2[,6]

n<-length(originalync)

tempZ1<-rep(1,n)
tempZ2<-BostonHousing2[,7]
tempZ3<-BostonHousing2[,12]
tempZ4<-BostonHousing2[,16]
tempZ5<-BostonHousing2[,11]

Z1nc<-tempZ1
Z2nc<-tempZ2
Z3nc<-tempZ3
Z4nc<-tempZ4
Z5nc<-tempZ5



originaly<-originalync-rep(mean(originalync),n)
Z1<-Z1nc
Z2<-(Z2nc-rep(mean(Z2nc),n))/sd(Z2nc)
Z3<-(Z3nc-rep(mean(Z3nc),n))/sd(Z3nc)
Z4<-(Z4nc-rep(mean(Z4nc),n))/sd(Z4nc)
Z5<-(Z5nc-rep(mean(Z5nc),n))/sd(Z5nc)


p<-1

# test grid 
testno<-50

tempgrid<-seq(min(originalx),max(originalx),by=(max(originalx)-min(originalx))/(200-1))

sequence<-seq(1,197,4)

grid<-tempgrid[sequence]

# cal number
calno<-100

# cal points
calpoints<-seq(min(originalx),max(originalx),by=(max(originalx)-min(originalx))/(calno-1))

Z<-cbind(Z1,Z2,Z3,Z4,Z5)

cvsplit <- function(n,folds,myseed)
{set.seed(myseed)
  subn = floor(n/folds)
  splitsize = subn*rep(1,folds);
  splitsize[folds]=n-(folds-1)*subn;
  permid = sample(1:n, replace=F)
  ll = n-splitsize[folds]+1
  findex = matrix(0,splitsize,folds-1)
  for (kk in 1:(folds-1))
  {findex[,kk]=permid[(kk-1)*subn+(1:subn)]}
  lastindex = permid[(subn*(folds-1)+1):n]
  return(list(findex=findex, lastindex=lastindex))
}

myseed<-2011
fold<-5
splits<-cvsplit(n,fold, myseed)

splitx<-splity<-splitz1<-splitz2<-splitz3<-splitz4<-splitz5<-list()
splittrainx<-splittrainy<-splittrainz1<-splittrainz2<-splittrainz3<-splittrainz4<-splittrainz5<-list()
for (kk in 1:(fold-1)){
  splitx[[kk]]<-originalx[splits$findex[,kk]]
  splity[[kk]]<-originaly[splits$findex[,kk]]
  splitz1[[kk]]<-Z1[splits$findex[,kk]]
  splitz2[[kk]]<-Z2[splits$findex[,kk]]
  splitz3[[kk]]<-Z3[splits$findex[,kk]]
  splitz4[[kk]]<-Z4[splits$findex[,kk]]
  splitz5[[kk]]<-Z5[splits$findex[,kk]]
  splittrainx[[kk]]<-originalx[-splits$findex[,kk]]
  splittrainy[[kk]]<-originaly[-splits$findex[,kk]]
  splittrainz1[[kk]]<-Z1[-splits$findex[,kk]]
  splittrainz2[[kk]]<-Z2[-splits$findex[,kk]]
  splittrainz3[[kk]]<-Z3[-splits$findex[,kk]]
  splittrainz4[[kk]]<-Z4[-splits$findex[,kk]]
  splittrainz5[[kk]]<-Z5[-splits$findex[,kk]]
}
splitx[[fold]]<-originalx[splits$lastindex]
splity[[fold]]<-originaly[splits$lastindex]
splitz1[[fold]]<-Z1[splits$lastindex]
splitz2[[fold]]<-Z2[splits$lastindex]
splitz3[[fold]]<-Z3[splits$lastindex]
splitz4[[fold]]<-Z4[splits$lastindex]
splitz5[[fold]]<-Z5[splits$lastindex]
splittrainx[[fold]]<-originalx[-splits$lastindex]
splittrainy[[fold]]<-originaly[-splits$lastindex]
splittrainz1[[fold]]<-Z1[-splits$lastindex]
splittrainz2[[fold]]<-Z2[-splits$lastindex]
splittrainz3[[fold]]<-Z3[-splits$lastindex]
splittrainz4[[fold]]<-Z4[-splits$lastindex]
splittrainz5[[fold]]<-Z5[-splits$lastindex]


T<-matrix(0,5,10)
T[1,1]<-T[2,3]<-T[3,5]<-T[4,7]<-T[5,9]<-1

hmin<-(max(originalx)-min(originalx))/10
hmax<-(max(originalx)-min(originalx))/2

hvalue<-exp(seq(log(hmin),log(hmax),length.out=12))

lengthh<-length(hvalue)

cverror<-rep(0,lengthh)

for (hs in 1:lengthh) {
  h<-hvalue[hs]
  
  K<-function(x) 
  {
    K<-(1/h)*(1-(x/h)^2)*((1-(x/h)^2)>0)
  }
  
  temperror<-rep(NA,fold)
  for (kk in 1:fold) {
    trainx<-splittrainx[[kk]]
    testx<-splitx[[kk]]
    trainy<-splittrainy[[kk]]
    testy<-splity[[kk]]
    trainz1<-splittrainz1[[kk]]
    testz1<-splitz1[[kk]]
    trainz2<-splittrainz2[[kk]]
    testz2<-splitz2[[kk]]
    trainz3<-splittrainz3[[kk]]
    testz3<-splitz3[[kk]]
    trainz4<-splittrainz4[[kk]]
    testz4<-splitz4[[kk]]
    trainz5<-splittrainz5[[kk]]
    testz5<-splitz5[[kk]]
    
    
    pilotdesignX1<-function(x) {
      pilotdesignX1<-matrix(0,length(trainx),(p+3))
      for (i in 1:length(trainx)){
        for (j in 1:(3+p)) {
          pilotdesignX1[i,j]<-trainz1[i]*(trainx[i]-x)^(j-1)
        }
      }
      return(pilotdesignX1)
    }
    
    pilotdesignX2<-function(x) {
      pilotdesignX2<-matrix(0,length(trainx),(p+3))
      for (i in 1:length(trainx)){
        for (j in 1:(3+p)) {
          pilotdesignX2[i,j]<-trainz2[i]*(trainx[i]-x)^(j-1)
        }
      }
      return(pilotdesignX2)
    }
    
    pilotdesignX3<-function(x) {
      pilotdesignX3<-matrix(0,length(trainx),(p+3))
      for (i in 1:length(trainx)){
        for (j in 1:(3+p)) {
          pilotdesignX3[i,j]<-trainz3[i]*(trainx[i]-x)^(j-1)
        }
      }
      return(pilotdesignX3)
    }
    
    pilotdesignX4<-function(x) {
      pilotdesignX4<-matrix(0,length(trainx),(p+3))
      for (i in 1:length(trainx)){
        for (j in 1:(3+p)) {
          pilotdesignX4[i,j]<-trainz4[i]*(trainx[i]-x)^(j-1)
        }
      }
      return(pilotdesignX4)
    }
    
    pilotdesignX5<-function(x) {
      pilotdesignX5<-matrix(0,length(trainx),(p+3))
      for (i in 1:length(trainx)){
        for (j in 1:(3+p)) {
          pilotdesignX5[i,j]<-trainz5[i]*(trainx[i]-x)^(j-1)
        }
      }
      return(pilotdesignX5)
    }
    
    
    pilotdesignX<-function(x) {
      pilotdesignX<-cbind(pilotdesignX1(x),pilotdesignX2(x),pilotdesignX3(x),pilotdesignX4(x),pilotdesignX5(x))
    }
    
    W<-function(x) {
      tempvector<-trainx-rep(x,length(trainx))
      W<-diag(apply(as.matrix(tempvector),1,K))
      return(W)
    }
    
    pilotsolution<-function(x) {
      temp<-t(pilotdesignX(x))%*%W(x)%*%pilotdesignX(x)
      if (det(temp)<=10^(-10))
        temp<-temp+10^(-8)*diag(1,5*(p+3))
      pilotsolution<-solve(temp)%*%t(pilotdesignX(x))%*%W(x)%*%trainy
    }
    
    pilottestestimate<-apply(as.matrix(testx),1,pilotsolution)
    
    temperror[kk]<-t(testy-testz1*pilottestestimate[1,]-testz2*pilottestestimate[5,]-testz3*pilottestestimate[9,]-testz4*pilottestestimate[13,]-testz5*pilottestestimate[17,])%*%(testy-testz1*pilottestestimate[1,]-testz2*pilottestestimate[5,]-testz3*pilottestestimate[9,]-testz4*pilottestestimate[13,]-testz5*pilottestestimate[17,])
  }
  cverror[hs]<-sum(temperror)
}
index<-which.min(cverror)
h<-hvalue[index]

pilotdesignX1<-function(x) {
  pilotdesignX1<-matrix(0,n,(p+3))
  for (i in 1:n){
    for (j in 1:(3+p)) {
      pilotdesignX1[i,j]<-Z1[i]*(originalx[i]-x)^(j-1)
    }
  }
  return(pilotdesignX1)
}

pilotdesignX2<-function(x) {
  pilotdesignX2<-matrix(0,n,(p+3))
  for (i in 1:n){
    for (j in 1:(3+p)) {
      pilotdesignX2[i,j]<-Z2[i]*(originalx[i]-x)^(j-1)
    }
  }
  return(pilotdesignX2)
}

pilotdesignX3<-function(x) {
  pilotdesignX3<-matrix(0,n,(p+3))
  for (i in 1:n){
    for (j in 1:(3+p)) {
      pilotdesignX3[i,j]<-Z3[i]*(originalx[i]-x)^(j-1)
    }
  }
  return(pilotdesignX3)
}

pilotdesignX4<-function(x) {
  pilotdesignX4<-matrix(0,n,(p+3))
  for (i in 1:n){
    for (j in 1:(3+p)) {
      pilotdesignX4[i,j]<-Z4[i]*(originalx[i]-x)^(j-1)
    }
  }
  return(pilotdesignX4)
}

pilotdesignX5<-function(x) {
  pilotdesignX5<-matrix(0,n,(p+3))
  for (i in 1:n){
    for (j in 1:(3+p)) {
      pilotdesignX5[i,j]<-Z5[i]*(originalx[i]-x)^(j-1)
    }
  }
  return(pilotdesignX5)
}


pilotdesignX<-function(x) {
  pilotdesignX<-cbind(pilotdesignX1(x),pilotdesignX2(x),pilotdesignX3(x),pilotdesignX4(x),pilotdesignX5(x))
}

W<-function(x) {
  tempvector<-originalx-rep(x,length(originalx))
  W<-diag(apply(as.matrix(tempvector),1,K))
  return(W)
}

pilotsolution<-function(x) {
  temp<-t(pilotdesignX(x))%*%W(x)%*%pilotdesignX(x)
  if (det(temp)<=10^(-10))
    temp<-temp+10^(-8)*diag(1,5*(p+3))
  pilotsolution<-solve(temp)%*%t(pilotdesignX(x))%*%W(x)%*%originaly
}

designX1<-function(x) {
  designX1<-matrix(0,n,(p+1))
  for (i in 1:n){
    for (j in 1:(1+p)) {
      designX1[i,j]<-Z1[i]*(originalx[i]-x)^(j-1)}
  }
  return(designX1)
}

designX2<-function(x) {
  designX2<-matrix(0,n,(p+1))
  for (i in 1:n){
    for (j in 1:(1+p)) {
      designX2[i,j]<-Z2[i]*(originalx[i]-x)^(j-1)}
  }
  return(designX2)
}

designX3<-function(x) {
  designX3<-matrix(0,n,(p+1))
  for (i in 1:n){
    for (j in 1:(1+p)) {
      designX3[i,j]<-Z3[i]*(originalx[i]-x)^(j-1)}
  }
  return(designX3)
}

designX4<-function(x) {
  designX4<-matrix(0,n,(p+1))
  for (i in 1:n){
    for (j in 1:(1+p)) {
      designX4[i,j]<-Z4[i]*(originalx[i]-x)^(j-1)}
  }
  return(designX4)
}

designX5<-function(x) {
  designX5<-matrix(0,n,(p+1))
  for (i in 1:n){
    for (j in 1:(1+p)) {
      designX5[i,j]<-Z5[i]*(originalx[i]-x)^(j-1)}
  }
  return(designX5)
}


designX<-function(x) {
  designX<-cbind(designX1(x),designX2(x),designX3(x),designX4(x),designX5(x))
}


pilotestimate<-apply(as.matrix(originalx),1,pilotsolution)

calestimate<-apply(as.matrix(calpoints),1,pilotsolution)

sigmaestimate<-rep(0,calno)

tau<-bias<-cov<-omega<-list()
for (cals in 1:calno) {
  calx<-calpoints[cals]
  designweightmatrix<-W(calx)
  
  tempX<-pilotdesignX(calx)
  
  tempdetmatrixX<-t(tempX)%*%designweightmatrix%*%tempX
  #if (det(tempdetmatrixX)<=10^(-15))
  #tempdetmatrixX<-tempdetmatrixX+10^(-6)*diag(1,ncol(tempX))
  tempmatrixinverse<-solve(tempdetmatrixX)
  tempmatrixsquare<-t(tempX)%*%designweightmatrix%*%designweightmatrix%*%tempX
  tempmatrix2<-tempmatrixinverse%*%tempmatrixsquare                                   
  trace<-sum(diag(designweightmatrix))-sum(diag(tempmatrix2))
  
  totalsum<-t(originaly)%*%(designweightmatrix-designweightmatrix%*%tempX%*%tempmatrixinverse%*%t(tempX)%*%designweightmatrix)%*%originaly
  
  sigmaestimate[cals]<-totalsum/trace
  
  caldesignmatrix<-pilotdesignX(calx)
  tempindex<-c(3,4,7,8,11,12,15,16,19,20)
  calfit<-caldesignmatrix[,tempindex]%*%calestimate[tempindex,cals]
  tau[[cals]]<-calfit
  omega[[cals]]<-(t(Z)%*%designweightmatrix%*%Z)/sum(diag(designweightmatrix))
}

hmin<-(max(originalx)-min(originalx))/10
hmax<-(max(originalx)-min(originalx))/2

finalhvalue<-exp(seq(log(hmin),log(hmax),length.out=20))
lengthfinalh<-length(finalhvalue)

IMSE<-rep(0,lengthfinalh)
for (fhs in 1:lengthfinalh){
  h<-finalhvalue[fhs]
  
  
  
  
  K<-function(x) 
  {
    K<-(1/h)*(1-(x/h)^2)*((1-(x/h)^2)>0)
  }
  
  W<-function(x) {
    tempvector<-originalx-rep(x,n)
    W<-diag(apply(as.matrix(tempvector),1,K))
    return(W)
  }
  
  tempMSE<-rep(0,calno)
  for (cals in 1:calno) {
    calx<-calpoints[cals]
    
    designweightmatrix<-W(calx)
    designmatrix<-designX(calx)
    temporarymatrix<-t(designmatrix)%*%designweightmatrix%*%designmatrix
    if (det(temporarymatrix)<=1)
      temporarymatrix<-temporarymatrix+10^(-8)*diag(1,10)
    
    matrixinverse<-solve(temporarymatrix)
    matrixsquare<-t(designmatrix)%*%designweightmatrix%*%designweightmatrix%*%designmatrix
    
    # covariance
    cov[[cals]]<-T%*%matrixinverse%*%matrixsquare%*%matrixinverse%*%t(T)*sigmaestimate[cals]
    
    # bias
    bias[[cals]]<-T%*%matrixinverse%*%t(designmatrix)%*%designweightmatrix%*%tau[[cals]]
    tempMSE[cals]<-t(bias[[cals]])%*%omega[[cals]]%*%bias[[cals]]+sum(diag(omega[[cals]]%*%cov[[cals]]))
  }
  IMSE[fhs]<-mean(tempMSE)
}
selectindex<-which.min(IMSE)

optimalindex<-which.min(IMSE)
optimalh<-finalhvalue[optimalindex]

h<-optimalh


K<-function(x) 
{
  K<-0.75*(1/h)*(1-(x/h)^2)*((1-(x/h)^2)>0)
}

W<-function(x) {
  tempvector<-originalx-rep(x,n)
  W<-diag(apply(as.matrix(tempvector),1,K))
  return(W)
}

solution<-function(x) {
  tempsolvematrix<-t(designX(x))%*%W(x)%*%designX(x)
  if (det(tempsolvematrix)<1)
    tempsolvematrix<-tempsolvematrix+10^(-8)*diag(1,10)
  solution<-solve(tempsolvematrix)%*%t(designX(x))%*%W(x)%*%originaly
}

# estimate of the design points
originalgridestimate1<-originalgridestimate2<-originalgridestimate3<-originalgridestimate4<-originalgridestimate5<-rep(0,testno)
for (ii in 1:testno) {
  tempestimate<-solution(grid[ii])
  originalgridestimate1[ii]<-tempestimate[1]
  originalgridestimate2[ii]<-tempestimate[p+2]
  originalgridestimate3[ii]<-tempestimate[p+4]
  originalgridestimate4[ii]<-tempestimate[p+6]
  originalgridestimate5[ii]<-tempestimate[p+8]
}

# testcode
#determinant<-rep(0,testno)
#for (ii in 1:testno) {
#tempsolvematrix<-t(designX(grid[ii]))%*%W(grid[ii])%*%designX(grid[ii])
#determinant[ii]<-det(tempsolvematrix)
#}

penalizedfestimate1<-penalizedfestimate2<-penalizedfestimate3<-penalizedfestimate4<-penalizedfestimate5<-rep(0,testno)
temprecordlambda<-rep(0,testno)

for (ii in 1:testno){
  fixx<-grid[ii]
  
  kernelestimate<-rep(0,n)
  for (i in 1:n) {
    kernelestimate[i]<-K(fixx-originalx[i])
  }
  
  samplen<-sum(kernelestimate)/0.75*h
  K0<-0.75/h
  
  OLSestimateoriginal<-solution(fixx)
  
  OLSestimate1<-OLSestimateoriginal[1:2]
  OLSestimate2<-OLSestimateoriginal[(p+2):(p+3)]
  OLSestimate3<-OLSestimateoriginal[(p+4):(p+5)]
  OLSestimate4<-OLSestimateoriginal[(p+6):(p+7)]
  OLSestimate5<-OLSestimateoriginal[(p+8):(p+9)]
  OLSestimate<-c(OLSestimate1,OLSestimate2,OLSestimate3,OLSestimate4,OLSestimate5)
  
  Wmatrix<-W(fixx)/K0
  scaleW<-sqrt(Wmatrix)
  weight<-rep(0,n)
  for (i in 1:n) {
    weight[i]<-scaleW[i,i]
  }
  
  Wmatrix<-W(fixx)/K0
  
  grplassoweight<-rep(0,n)
  for (i in 1:n) {
    grplassoweight[i]<-Wmatrix[i,i]
  }
  
  tempX<-designX(fixx)
  
  tempmatrix<-(1/sqrt(2))*diag(c(1,1,1,1,1,1,1,1,1,1))
  
  X<-tempX
  
  
  lambdavalue<-seq(-5,4,0.05)
  lambdalength<-length(lambdavalue)
  
  pseudoy<-scaleW%*%originaly
  pseudoX<-scaleW%*%X
  
  newOLSestimate<-as.vector(OLSestimate%*%solve(tempmatrix))
  newOLSestimate1<-newOLSestimate[1:2]
  newOLSestimate2<-newOLSestimate[(p+2):(p+3)]
  newOLSestimate3<-newOLSestimate[(p+4):(p+5)]
  newOLSestimate4<-newOLSestimate[(p+6):(p+7)]
  newOLSestimate5<-newOLSestimate[(p+8):(p+9)]
  standardsd<-apply(pseudoX,2,sd)
  newX<-pseudoX%*%solve(diag(standardsd))%*%tempmatrix
  
  BIC<-RSS<-rep(0,lambdalength)
  
  
  for (vv in 1:lambdalength){
    lambdatune<-exp(lambdavalue[vv])
    
    w1<-1
    
    w2<-dp(lambdatune,sqrt(sum(newOLSestimate2^2)))
    
    w3<-dp(lambdatune,sqrt(sum(newOLSestimate3^2)))
    
    w4<-dp(lambdatune,sqrt(sum(newOLSestimate4^2)))
    
    w5<-dp(lambdatune,sqrt(sum(newOLSestimate5^2)))
    
    indexvector<-c(rep(NA,p+1),rep(2,p+1),rep(3,p+1),rep(4,p+1),rep(5,p+1))
    
    w<-c(w1,w1,w2,w2,w3,w3,w4,w4,w5,w5)
    
    # test code
    if (w2==0) {
      indexvector[3:4]<-c(NA,NA)
      w[3:4]<-c(1,1)
    }
    
    # test code
    if (w3==0) {
      indexvector[5:6]<-c(NA,NA)
      w[5:6]<-c(1,1)
    }
    
    # test code
    if (w4==0) {
      indexvector[7:8]<-c(NA,NA)
      w[7:8]<-c(1,1)
    }
    
    # test code
    if (w5==0) {
      indexvector[9:10]<-c(NA,NA)
      w[9:10]<-c(1,1)
    }
    
    if (sum(w)==10){
      betasolve<-newOLSestimate
    }
    
    if (sum(w)!=10){
      
      tempnewX<-matrix(0,nrow(newX),ncol(newX))
      for (j in 1:ncol(newX)){
        tempnewX[,j]<-newX[,j]/w[j]
      }
      grouplasso <- grplasso(tempnewX, pseudoy,index=indexvector, coef.init=newOLSestimate*w,lambda=samplen, model = LinReg(),
                             center = FALSE, standardize = FALSE, control = grpl.control(trace=0, max.iter = 5000))
      betasolve<-rep(0,ncol(newX))
      for (j in 1:ncol(newX)){
        betasolve[j]<-grouplasso$coefficients[j]/w[j]
      }
    }
    
    betasolve1<-betasolve[1:2]
    betasolve2<-betasolve[3:4]
    betasolve3<-betasolve[5:6]
    betasolve4<-betasolve[7:8]
    betasolve5<-betasolve[9:10]
    
    
    
    nonzerono<-0
    if (min(abs(betasolve2))>0 )
      nonzerono<-nonzerono+1
    if (min(abs(betasolve3))>0 )
      nonzerono<-nonzerono+1
    if (min(abs(betasolve4))>0 )
      nonzerono<-nonzerono+1
    if (min(abs(betasolve5))>0 )
      nonzerono<-nonzerono+1
    
    
    
    
    df<-nonzerono+p*(sum(betasolve2^2))/(sum(newOLSestimate2^2))+p*(sum(betasolve3^2))/(sum(newOLSestimate3^2))+p*(sum(betasolve4^2))/(sum(newOLSestimate4^2))+p*(sum(betasolve5^2))/(sum(newOLSestimate5^2))
    
    RSS[vv]<-t(pseudoy-newX%*%betasolve)%*%(pseudoy-newX%*%betasolve)
    BIC[vv]<-samplen*log(RSS[vv]/samplen)+log(samplen)*df
    
    # testing code
    #t(pseudoy)%*%pseudoy
    #t(pseudoy-newX[,1:2]%*%betasolve[1:2])%*%(pseudoy-newX[,1:2]%*%betasolve[1:2])
    #t(pseudoy-newX[,3:4]%*%betasolve[3:4])%*%(pseudoy-newX[,3:4]%*%betasolve[3:4])
    #t(pseudoy-newX[,5:6]%*%betasolve[5:6])%*%(pseudoy-newX[,5:6]%*%betasolve[5:6])
    #t(pseudoy-newX[,1:2]%*%betasolve[1:2]-newX[,3:4]%*%betasolve[3:4]-newX[,5:6]%*%betasolve[5:6])%*%(pseudoy-newX[,1:2]%*%betasolve[1:2]-newX[,3:4]%*%betasolve[3:4]-newX[,5:6]%*%betasolve[5:6])
  }
  
  selectno<-which.min(BIC)
  lambda<-exp(lambdavalue[selectno])
  
  temprecordlambda[ii]<-lambda
  
  w1<-1
  
  w2<-dp(temprecordlambda[ii],sqrt(sum(newOLSestimate2^2)))
  
  w3<-dp(temprecordlambda[ii],sqrt(sum(newOLSestimate3^2)))
  
  w4<-dp(temprecordlambda[ii],sqrt(sum(newOLSestimate4^2)))
  
  w5<-dp(temprecordlambda[ii],sqrt(sum(newOLSestimate5^2)))
  
  
  w<-c(w1,w1,w2,w2,w3,w3,w4,w4,w5,w5)
  
  # test code
  if (w2==0) {
    indexvector[3:4]<-c(NA,NA)
    w[3:4]<-c(1,1)
  }
  
  # test code
  if (w3==0) {
    indexvector[5:6]<-c(NA,NA)
    w[5:6]<-c(1,1)
  }
  
  # test code
  if (w4==0) {
    indexvector[7:8]<-c(NA,NA)
    w[7:8]<-c(1,1)
  }
  
  # test code
  if (w5==0) {
    indexvector[9:10]<-c(NA,NA)
    w[9:10]<-c(1,1)
  }
  
  if (sum(w)==10){
    betasolve<-newOLSestimate
  }
  
  if (sum(w)!=10){
    
    tempnewX<-matrix(0,nrow(newX),ncol(newX))
    for (j in 1:ncol(newX)){
      tempnewX[,j]<-newX[,j]/w[j]
    }
    grouplasso <- grplasso(tempnewX, pseudoy,index=indexvector, coef.init=newOLSestimate*w,lambda=samplen, model = LinReg(),
                           center = FALSE, standardize = FALSE, control = grpl.control(trace=0, max.iter = 5000))
    betasolve<-rep(0,ncol(newX))
    for (j in 1:ncol(newX)){
      betasolve[j]<-grouplasso$coefficients[j]/w[j]
    }
  }
  
  finalbeta<-t(solve(diag(standardsd))%*%betasolve)%*%tempmatrix
  penalizedfestimate1[ii]<-finalbeta[1]
  penalizedfestimate2[ii]<-finalbeta[p+2]
  penalizedfestimate3[ii]<-finalbeta[p+4]
  penalizedfestimate4[ii]<-finalbeta[p+6]
  penalizedfestimate5[ii]<-finalbeta[p+8]
}

finalresult<-rbind(penalizedfestimate1, penalizedfestimate2, penalizedfestimate3, penalizedfestimate4, penalizedfestimate5, originalgridestimate1,originalgridestimate2,originalgridestimate3,originalgridestimate4,originalgridestimate5)
save(finalresult,file="Bostonhousingallintscale.dat")

#load("NMMAPSvaryingtemperaturehumnew.dat")
#originalgridestimate1<-finalresult[3,]
#penalizedfestimate1<-finalresult[1,]
#originalgridestimate2<-finalresult[4,]
#penalizedfestimate2<-finalresult[2,]

postscript(file = "BostonhousingCRIMintscalefunction150.eps")
plot(grid, penalizedfestimate1,type="l",lty=1,col="black")
title("Function 1")
dev.off()

postscript(file = "BostonhousingCRIMintscalefunction250.eps")
plot(grid, penalizedfestimate2,type="l", lty=1,col="black")
title("Function 2")
dev.off()

postscript(file = "BostonhousingCRIMintscalefunction350.eps")
plot(grid, penalizedfestimate3,type="l", lty=1,col="black")
title("Function 3")
dev.off()

postscript(file = "BostonhousingCRIMintscalefunction450.eps")
plot(grid, penalizedfestimate4,type="l", lty=1,col="black")
title("Function 4")
dev.off()

postscript(file = "BostonhousingCRIMintscalefunction550.eps")
plot(grid, penalizedfestimate5,type="l", lty=1,col="black")
title("Function 5")
dev.off()
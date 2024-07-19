rm(list=ls(all=TRUE))
xx <- read.table(file = "D:\\Data\\DeHartSimplified.csv",
header=TRUE , sep= ",", na.strings = " ")
y=xx$numall
y[is.na(y)]<-0
yy=y
  x1=xx$negevent
  x2=xx$nrel

fit <-glm(formula =y ~x1+x2,family=poisson)

##############################################
################################### ESTIMATION

################################## POISSON
lfpois=function(par) {
  lap=exp(par[1]+par[2]*x1+par[3]*x2)
  dp <- dpois(y,lap)           
  -sum(log(dp))
}

 init=c(0.1,0.1,0.1)
 outp=optim(init, lfpois)
################################## k=0
lpmf0 <- function(par) {
  ww0=(1-par[1])
  la0=exp(par[2]+par[3]*x1+par[4]*x2)
  d0 <- (par[1] * (y == 0) +  ww0 * dpois(y,la0))            
  -sum(log(d0))
}


 c_i=c(-1,0)
 u_i=rbind(c(-1,0,0,0),c(1,0,0,0))
 init=c(.2,.1,-.1,-.1)
 out0=constrOptim(init, lpmf0, NULL, ui=u_i, ci=c_i)

################################## k=1
lpmf1 <- function(par) {
  ww1=1-(par[1]+par[2])
  la1=exp(par[3]+par[4]*x1+par[5]*x2)
  d1 <- (par[1] * (y == 0) + par[2] * (y == 1)+  ww1 * dpois(y,la1))            
  -sum(log(d1))
}


 c_i=c(-1,0,0)
 u_i=rbind(c(-1,-1,0,0,0),c(1,0,0,0,0),c(0,1,0,0,0))
 init=c(.2,.2,.1,-.1,-.1)
 out1=constrOptim(init, lpmf1, NULL, ui=u_i, ci=c_i)

################################## k=2


lpmf2 <- function(par) {
  ww2=1-(par[1]+par[2]+par[3])
  la2=exp(par[4]+par[5]*x1+par[6]*x2)
  d2 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
  ww2 * dpois(y,la2))            
  -sum(log(d2))
}


 c_i=c(-1,0,0,0)
 u_i=rbind(c(-1,-1,-1,0,0,0),c(1,0,0,0,0,0),c(0,1,0,0,0,0),c(0,0,1,0,0,0))
 init=c(.2,.2,.2,.1,-.1,-.1)
 out2=constrOptim(init, lpmf2, NULL, ui=u_i, ci=c_i)

################################## k=3

lpmf3 <- function(par) {
  ww3=1-(par[1]+par[2]+par[3]+par[4])
  la3=exp(par[5]+par[6]*x1+par[7]*x2)
  d3 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
 par[4] * (y == 3)+ ww3 * dpois(y,la3))            
  -sum(log(d3))
}


 c_i=c(-1,0,0,0,0)
 u_i=rbind(c(-1,-1,-1,-1,0,0,0),c(1,0,0,0,0,0,0),c(0,1,0,0,0,0,0),
c(0,0,1,0,0,0,0),c(0,0,0,1,0,0,0))
 init=c(.2,.2,.2,.2,.1,-.1,-.1)
 out3=constrOptim(init, lpmf3, NULL, ui=u_i, ci=c_i)
################################## k=4

lpmf4 <- function(par) {
  ww4=1-(par[1]+par[2]+par[3]+par[4]+par[5])
  la4=exp(par[6]+par[7]*x1+par[8]*x2)
  d4 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
 par[4] * (y == 3) + par[5] * (y == 4)+ ww4 * dpois(y,la4))            
  -sum(log(d4))
}


 c_i=c(-1,0,0,0,0,0)
 u_i=rbind(c(-1,-1,-1,-1,-1,0,0,0),c(1,0,0,0,0,0,0,0),c(0,1,0,0,0,0,0,0),
c(0,0,1,0,0,0,0,0),c(0,0,0,1,0,0,0,0),c(0,0,0,0,1,0,0,0))
 init=c(.1,.1,.1,.1,.1,.1,-.1,-.1)
 out4=constrOptim(init, lpmf4, NULL, ui=u_i, ci=c_i)
################################## k=5

lpmf5 <- function(par) {
  ww5=1-(par[1]+par[2]+par[3]+par[4]+par[5]+par[6])
  la5=exp(par[7]+par[8]*x1+par[9]*x2)
  d5 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
 par[4] * (y == 3) + par[5] * (y == 4)+ par[6] * (y == 5)+ ww5 * dpois(y,la5))            
  -sum(log(d5))
}


 c_i=c(-1,0,0,0,0,0,0)
 u_i=rbind(c(-1,-1,-1,-1,-1,-1,0,0,0),c(1,0,0,0,0,0,0,0,0),
c(0,1,0,0,0,0,0,0,0),c(0,0,1,0,0,0,0,0,0),c(0,0,0,1,0,0,0,0,0),
c(0,0,0,0,1,0,0,0,0),c(0,0,0,0,0,1,0,0,0))
 init=c(.1,.1,.1,.1,.1,.1,.1,-.1,-.1)
 out5=constrOptim(init, lpmf5, NULL, ui=u_i, ci=c_i)
################################## k=6

lpmf6 <- function(par) {
  ww6=1-(par[1]+par[2]+par[3]+par[4]+par[5]+par[6]+par[7])
  la6=exp(par[8]+par[9]*x1+par[10]*x2)
  d6 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
 par[4] * (y == 3) + par[5] * (y == 4)+ par[6] * (y == 5) + 
 par[7] * (y == 6)+ ww6 * dpois(y,la6))            
  -sum(log(d6))
}


 c_i=c(-1,0,0,0,0,0,0,0)
 u_i=rbind(c(-1,-1,-1,-1,-1,-1,-1,0,0,0),c(1,0,0,0,0,0,0,0,0,0),
c(0,1,0,0,0,0,0,0,0,0),c(0,0,1,0,0,0,0,0,0,0),c(0,0,0,1,0,0,0,0,0,0),
c(0,0,0,0,1,0,0,0,0,0),c(0,0,0,0,0,1,0,0,0,0),c(0,0,0,0,0,0,1,0,0,0))
 init=c(.1,.1,.1,.1,.1,.1,.1,.1,-.1,-.1)
 out6=constrOptim(init, lpmf6, NULL, ui=u_i, ci=c_i)
################################## k=7

lpmf7 <- function(par) {
  ww7=1-(par[1]+par[2]+par[3]+par[4]+par[5]+par[6]+par[7]++par[8])
  la7=exp(par[9]+par[10]*x1+par[11]*x2)
  d7 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
 par[4] * (y == 3) + par[5] * (y == 4)+ par[6] * (y == 5) + 
 par[7] * (y == 6) +  par[8] * (y == 7) + ww7 * dpois(y,la7))            
  -sum(log(d7))
}


 c_i=c(-1,0,0,0,0,0,0,0,0)
 u_i=rbind(c(-1,-1,-1,-1,-1,-1,-1,-1,0,0,0),c(1,0,0,0,0,0,0,0,0,0,0),
c(0,1,0,0,0,0,0,0,0,0,0),c(0,0,1,0,0,0,0,0,0,0,0),c(0,0,0,1,0,0,0,0,0,0,0),
c(0,0,0,0,1,0,0,0,0,0,0),c(0,0,0,0,0,1,0,0,0,0,0),c(0,0,0,0,0,0,1,0,0,0,0),
c(0,0,0,0,0,0,0,1,0,0,0))
 init=c(.1,.1,.1,.1,.1,.1,.1,.1,.1,-.1,-.1)
 out7=constrOptim(init, lpmf7, NULL, ui=u_i, ci=c_i)
################################## k=8

lpmf8 <- function(par) {
  ww8=1-(par[1]+par[2]+par[3]+par[4]+par[5]+par[6]+par[7]+par[8]+par[9])
  la8=exp(par[10]+par[11]*x1+par[12]*x2)
  d8 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
 par[4] * (y == 3) + par[5] * (y == 4)+ par[6] * (y == 5) + 
 par[7] * (y == 6) +  par[8] * (y == 7) +  par[9] * (y == 8) +
 ww8 * dpois(y,la8))            
  -sum(log(d8))
}


 c_i=c(-1,0,0,0,0,0,0,0,0,0)
 u_i=rbind(
c(-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0),c(1,0,0,0,0,0,0,0,0,0,0,0),
c(0,1,0,0,0,0,0,0,0,0,0,0),c(0,0,1,0,0,0,0,0,0,0,0,0),
c(0,0,0,1,0,0,0,0,0,0,0,0),c(0,0,0,0,1,0,0,0,0,0,0,0),
c(0,0,0,0,0,1,0,0,0,0,0,0),c(0,0,0,0,0,0,1,0,0,0,0,0),
c(0,0,0,0,0,0,0,1,0,0,0,0),c(0,0,0,0,0,0,0,0,1,0,0,0))
 init=c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,-.1,-.1)
 out8=constrOptim(init, lpmf8, NULL, ui=u_i, ci=c_i)
################################## k=9

lpmf9 <- function(par) {
  ww9=1-(par[1]+par[2]+par[3]+par[4]+par[5]+par[6]+par[7]+par[8]+par[9]+
      par[10])
  la9=exp(par[11]+par[12]*x1+par[13]*x2)
  d9 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
 par[4] * (y == 3) + par[5] * (y == 4)+ par[6] * (y == 5) + 
 par[7] * (y == 6) +  par[8] * (y == 7) +  par[9] * (y == 8) +
 par[10] * (y == 9) + ww9 * dpois(y,la9))            
  -sum(log(d9))
}


 c_i=c(-1,0,0,0,0,0,0,0,0,0,0)
 u_i=rbind(
c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0),c(1,0,0,0,0,0,0,0,0,0,0,0,0),
c(0,1,0,0,0,0,0,0,0,0,0,0,0),c(0,0,1,0,0,0,0,0,0,0,0,0,0),
c(0,0,0,1,0,0,0,0,0,0,0,0,0),c(0,0,0,0,1,0,0,0,0,0,0,0,0),
c(0,0,0,0,0,1,0,0,0,0,0,0,0),c(0,0,0,0,0,0,1,0,0,0,0,0,0),
c(0,0,0,0,0,0,0,1,0,0,0,0,0),c(0,0,0,0,0,0,0,0,1,0,0,0,0),
c(0,0,0,0,0,0,0,0,0,1,0,0,0))
 init=c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.1,-.1,-.1)
 out9=constrOptim(init, lpmf9, NULL, ui=u_i, ci=c_i)
################################## k=10

lpmf10 <- function(par) {
  ww10=1-(par[1]+par[2]+par[3]+par[4]+par[5]+par[6]+par[7]+par[8]+par[9]+
      par[10]+par[11])
  la10=exp(par[12]+par[13]*x1+par[14]*x2)
  d10 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
 par[4] * (y == 3) + par[5] * (y == 4)+ par[6] * (y == 5) + 
 par[7] * (y == 6) +  par[8] * (y == 7) +  par[9] * (y == 8) +
 par[10] * (y == 9)+ par[11] * (y == 10) + ww10 * dpois(y,la10))            
  -sum(log(d10))
}


 c_i=c(-1,0,0,0,0,0,0,0,0,0,0,0)
 u_i=rbind(
c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0),c(1,0,0,0,0,0,0,0,0,0,0,0,0,0),
c(0,1,0,0,0,0,0,0,0,0,0,0,0,0),c(0,0,1,0,0,0,0,0,0,0,0,0,0,0),
c(0,0,0,1,0,0,0,0,0,0,0,0,0,0),c(0,0,0,0,1,0,0,0,0,0,0,0,0,0),
c(0,0,0,0,0,1,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,1,0,0,0,0,0,0,0),
c(0,0,0,0,0,0,0,1,0,0,0,0,0,0),c(0,0,0,0,0,0,0,0,1,0,0,0,0,0),
c(0,0,0,0,0,0,0,0,0,1,0,0,0,0),c(0,0,0,0,0,0,0,0,0,0,1,0,0,0))
 init=c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.1,-.1,-.1)
 out10=constrOptim(init, lpmf10, NULL, ui=u_i, ci=c_i)

################################## k=11

lpmf11 <- function(par) {
  ww11=1-(par[1]+par[2]+par[3]+par[4]+par[5]+par[6]+par[7]+par[8]+par[9]+
      par[10]+par[11]+par[12])
  la11=exp(par[13]+par[14]*x1+par[15]*x2)
  d11 <- (par[1] * (y == 0) + par[2] * (y == 1)+ par[3] * (y == 2) +
 par[4] * (y == 3) + par[5] * (y == 4)+ par[6] * (y == 5) + 
 par[7] * (y == 6) +  par[8] * (y == 7) +  par[9] * (y == 8) +
 par[10] * (y == 9)+ par[11] * (y == 10) +
+ par[12] * (y == 11) + ww11 * dpois(y,la11))            
  -sum(log(d11))
}


 c_i=c(-1,0,0,0,0,0,0,0,0,0,0,0,0)
 u_i=rbind(
c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0),c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0))
 init=c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.1,-.1,-.1)
 out11=constrOptim(init, lpmf11, NULL, ui=u_i, ci=c_i)


############################################
#################################### COPARE
# AIC
aic=c()
aic[1]=2*3+2*outp$value
aic[2]=2*4+2*out0$value
aic[3]=2*5+2*out1$value
aic[4]=2*6+2*out2$value
aic[5]=2*7+2*out3$value
aic[6]=2*8+2*out4$value
aic[7]=2*9+2*out5$value
aic[8]=2*10+2*out6$value
aic[9]=2*11+2*out7$value
aic[10]=2*12+2*out8$value
aic[11]=2*13+2*out9$value
aic[12]=2*14+2*out10$value
aic[13]=2*15+2*out11$value
k=(1:13)
plot(k,aic)

#########################################################################
############################################Randomitze Quantile Residuals

################################## POISSON
n=length(y)
lap=exp(outp$par[1]+outp$par[2]*x1+outp$par[3]*x2)

rqp=ppois(y-1,lap) + dpois(y,lap) * runif(n)

################################## k=0
par=out0$par
pmf0 <- function(y,la) {
  ww0=(1-par[1])
  d0 <- (par[1] * (y == 0) +  ww0 * dpois(y,la))            
  d0
}
CDF0=function(y,la){
 n=length(y)
 cdf0=rep(0,n)
 for(i in 1:n)
  {
   cdf0[i]=0
   if(y[i]<0) cdf0[i]=0
   else if(y[i]>-1){
     for(j in 0:y[i]) cdf0[i]=cdf0[i]+pmf0(j,la[i])
    }
   }
 cdf0
}
la0=exp(par[2]+par[3]*x1+par[4]*x2)
rq0=CDF0(y-1,la0) + pmf0(y,la0) * runif(n)


################################## k=1
par1=out1$par
pmf1 <- function(y,la) {
  ww1=1-(par1[1]+par1[2])
  d1 <- (par1[1] * (y == 0) + par1[2] * (y == 1)+  ww1 * dpois(y,la))            
  d1
}

CDF1=function(y,la){
 n=length(y)
 cdf1=rep(0,n)
 for(i in 1:n)
  {
   cdf1[i]=0
   if(y[i]<0) cdf1[i]=0
   else if(y[i]>-1){
     for(j in 0:y[i]) cdf1[i]=cdf1[i]+pmf1(j,la[i])
    }
   }
 cdf1
}

la1=exp(par1[3]+par1[4]*x1+par1[5]*x2)
rq1=CDF1(y-1,la1) + pmf1(y,la1) * runif(n)


################################## k=2

par2=out2$par
pmf2 <- function(y,la) {
  ww2=1-(par2[1]+par2[2]+par2[3])
  d2 <- (par2[1] * (y == 0) + par2[2] * (y == 1)+ par2[3] * (y == 2) +
  ww2 * dpois(y,la))            
  d2
}

CDF2=function(y,la){
 n=length(y)
 cdf2=rep(0,n)
 for(i in 1:n)
  {
   cdf2[i]=0
   if(y[i]<0) cdf2[i]=0
   else if(y[i]>-1){
     for(j in 0:y[i]) cdf2[i]=cdf2[i]+pmf2(j,la[i])
    }
   }
 cdf2
}

la2=exp(par2[4]+par2[5]*x1+par2[6]*x2)
rq2=CDF2(y-1,la2) + pmf2(y,la2) * runif(n)

################################## k=3

par3=out3$par
pmf3 <- function(y,la) {
  ww3=1-(par3[1]+par3[2]+par3[3]+par3[4])
  d3 <- (par3[1] * (y == 0) + par3[2] * (y == 1)+ par3[3] * (y == 2) +
 par3[4] * (y == 3)+ ww3 * dpois(y,la))            
  d3
}

CDF3=function(y,la){
 n=length(y)
 cdf3=rep(0,n)
 for(i in 1:n)
  {
   cdf3[i]=0
   if(y[i]<0) cdf3[i]=0
   else if(y[i]>-1){
     for(j in 0:y[i]) cdf3[i]=cdf3[i]+pmf3(j,la[i])
    }
   }
 cdf3
}

la3=exp(par3[5]+par3[6]*x1+par3[7]*x2)
rq3=CDF3(y-1,la3) + pmf3(y,la3) * runif(n)

################################## k=4

par4=out4$par
pmf4 <- function(y,la) {
  ww4=1-(par4[1]+par4[2]+par4[3]+par4[4]+par4[5])
  d4 <- (par4[1] * (y == 0) + par4[2] * (y == 1)+ par4[3] * (y == 2) +
 par4[4] * (y == 3) + par4[5] * (y == 4)+ ww4 * dpois(y,la))            
  d4
}

CDF4=function(y,la){
 n=length(y)
 cdf4=rep(0,n)
 for(i in 1:n)
  {
   cdf4[i]=0
   if(y[i]<0) cdf4[i]=0
   else if(y[i]>-1){
     for(j in 0:y[i]) cdf4[i]=cdf4[i]+pmf4(j,la[i])
    }
   }
 cdf4
}

la4=exp(par4[6]+par4[7]*x1+par4[8]*x2)
rq4=CDF4(y-1,la4) + pmf4(y,la4) * runif(n)

################################## k=5

par5=out5$par
pmf5 <- function(y,la) {
  ww5=1-(sum(par5[1:6]))
  d5 <- (par5[1] * (y == 0) + par5[2] * (y == 1)+ par5[3] * (y == 2) +
 par5[4] * (y == 3)+ par5[5] * (y == 4) + par5[6] * (y == 5)+ ww5 * dpois(y,la))            
  d5
}

CDF5=function(y,la){
 n=length(y)
 cdf5=rep(0,n)
 for(i in 1:n)
  {
   cdf5[i]=0
   if(y[i]<0) cdf5[i]=0
   else if(y[i]>-1){
     for(j in 0:y[i]) cdf5[i]=cdf5[i]+pmf5(j,la[i])
    }
   }
 cdf5
}

la5=exp(par5[7]+par5[8]*x1+par5[9]*x2)
rq5=CDF5(y-1,la5) + pmf4(y,la5) * runif(n)

################################## k=6

par6=out6$par
pmf6 <- function(y,la) {
  ww6=1-(sum(par6[1:7]))
  d6 <- (par6[1] * (y == 0) + par6[2] * (y == 1)+ par6[3] * (y == 2) +
 par6[4] * (y == 3) + par6[5] * (y == 4)+
 par6[6] * (y == 5) + par6[7] * (y == 6) + ww6 * dpois(y,la))            
  d6
}

CDF6=function(y,la){
 n=length(y)
 cdf6=rep(0,n)
 for(i in 1:n)
  {
   cdf6[i]=0
   if(y[i]<0) cdf6[i]=0
   else if(y[i]>-1){
     for(j in 0:y[i]) cdf6[i]=cdf6[i]+pmf6(j,la[i])
    }
   }
 cdf6
}

la6=exp(par6[8]+par6[9]*x1+par6[10]*x2)
rq6=CDF6(y-1,la6) + pmf6(y,la6) * runif(n)

################################## k=7

par7=out7$par
pmf7 <- function(y,la) {
  ww7=1-(sum(par7[1:8]))
  d7 <- (par7[1] * (y == 0) + par7[2] * (y == 1)+ par7[3] * (y == 2) +
 par7[4] * (y == 3) + par7[5] * (y == 4)+
 par7[6] * (y == 5) + par7[7] * (y == 6) + par7[8] * (y == 7)+
 ww7 * dpois(y,la))            
  d7
}

CDF7=function(y,la){
 n=length(y)
 cdf7=rep(0,n)
 for(i in 1:n)
  {
   cdf7[i]=0
   if(y[i]<0) cdf7[i]=0
   else if(y[i]>-1){
     for(j in 0:y[i]) cdf7[i]=cdf7[i]+pmf7(j,la[i])
    }
   }
 cdf7
}

la7=exp(par7[9]+par7[10]*x1+par7[11]*x2)
rq7=CDF7(y-1,la7) + pmf7(y,la7) * runif(n)

################################## k=8

par8=out8$par
pmf8 <- function(y,la) {
  ww8=1-(sum(par8[1:9]))
  d8 <- (par8[1] * (y == 0) + par8[2] * (y == 1)+ par8[3] * (y == 2) +
 par8[4] * (y == 3) + par8[5] * (y == 4)+
 par8[6] * (y == 5) + par8[7] * (y == 6)+
 par8[8] * (y == 7)+  par8[9] * (y == 8) + ww8 * dpois(y,la))            
  d8
}

CDF8=function(y,la){
 n=length(y)
 cdf8=rep(0,n)
 for(i in 1:n)
  {
   cdf8[i]=0
   if(y[i]<0) cdf8[i]=0
   else if(y[i]>-1){
     for(j in 0:y[i]) cdf8[i]=cdf8[i]+pmf8(j,la[i])
    }
   }
 cdf8
}

la8=exp(par8[10]+par8[11]*x1+par8[12]*x2)
rq8=CDF8(y-1,la8) + pmf8(y,la8) * runif(n)

####################Randomized Quantile
plot(qnorm(rqp),xlab="Observation Number",ylab="Randomized Quantile",main="Poisson")
plot(qnorm(rq0),xlab="Observation Number",ylab="Randomized Quantile",main="k=0")
plot(qnorm(rq1),xlab="Observation Number",ylab="Randomized Quantile",main="k=1")
plot(qnorm(rq2),xlab="Observation Number",ylab="Randomized Quantile",main="k=2")
plot(qnorm(rq3),xlab="Observation Number",ylab="Randomized Quantile",main="k=3")
plot(qnorm(rq4),xlab="Observation Number",ylab="Randomized Quantile",main="k=4")
plot(qnorm(rq5),xlab="Observation Number",ylab="Randomized Quantile",main="k=5")
plot(qnorm(rq6),xlab="Observation Number",ylab="Randomized Quantile",main="k=6")
plot(qnorm(rq7),xlab="Observation Number",ylab="Randomized Quantile",main="k=7")
plot(qnorm(rq8),xlab="Observation Number",ylab="Randomized Quantile",main="k=8")

hist(rqp,xlab="p-value",main="Randomized Quantile")
hist(rq0,xlab="p-value",main="Randomized Quantile")
hist(rq1,xlab="p-value",main="Randomized Quantile")
hist(rq2,xlab="p-value",main="Randomized Quantile")
hist(rq3,xlab="p-value",main="Randomized Quantile")
hist(rq4,xlab="p-value",main="Randomized Quantile")
hist(rq5,xlab="p-value",main="Randomized Quantile")
hist(rq6,xlab="p-value",main="Randomized Quantile")
hist(rq7,xlab="p-value",main="Randomized Quantile")
hist(rq8,xlab="p-value",main="Randomized Quantile")

#######################################
################################QQ-plot

qqnorm(qnorm(rqp),main="Randomized Quantile,Poisson",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rqp),main="Randomized Quantile",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

qqnorm(qnorm(rq0),main="Randomized Quantile, k=0",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rq0),main="Randomized Quantile, k=0",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

qqnorm(qnorm(rq1),main="Randomized Quantile, k=1",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rq1),main="Randomized Quantile, k=1",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

qqnorm(qnorm(rq2),main="Randomized Quantile, k=2",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rq2),main="Randomized Quantile, k=2",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

qqnorm(qnorm(rq3),main="Randomized Quantile, k=3",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rq3),main="Randomized Quantile, k=3",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

qqnorm(qnorm(rq4),main="Randomized Quantile, k=4",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rq4),main="Randomized Quantile",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

qqnorm(qnorm(rq5),main="Randomized Quantile, k=5",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rq5),main="Randomized Quantile, k=5",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

qqnorm(qnorm(rq6),main="Randomized Quantile, k=6",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rq6),main="Randomized Quantile, k=6",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

qqnorm(qnorm(rq7),main="Randomized Quantile, k=7",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rq7),main="Randomized Quantile, k=7",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

qqnorm(qnorm(rq8[rq8!=1]),main="Randomized Quantile, k=8",ylab="Sample Quantiles",xlab="Theoretical Quantiles")
qqline(qnorm(rq8),main="Randomized Quantile, k=8",ylab="Sample Quantiles",xlab="Theoretical Quantiles")

p=out7$par


rzkip <- function(n,par) {
la=(exp(par[9]+par[10]*x1+par[11]*x2))
  p=c(par[1:8],(1-sum(par[1:8])))
  out=rep(0,n)
  for(i in 1:n)
   {
    ran=sample(c(0,1,2,3,4,5,6,7,8), size=1, prob=p)
    if(ran==0) out[i]=0
    else if(ran==1) out[i]=1
    else if(ran==2) out[i]=2
    else if(ran==3) out[i]=3
    else if(ran==4) out[i]=4
    else if(ran==5) out[i]=5
    else if(ran==6) out[i]=6
    else if(ran==7) out[i]=7
    else if(ran==8) out[i]=rpois(1,la[i])
   }
return(out)
}
table(rzkip(length(yy),p))

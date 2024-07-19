rm(list=ls(all=TRUE))
library(EnvStats)
library(extraDistr)
library(VGAM)
library(rando)
library(data.table)

data = read.csv('d://Data/pap.csv')
y = data$PAP6YR1 # Number of Pap tests, last 6 years
x1 = data$R_MARITL # Marital Status
x2 = data$HPVPAP # Had HPV test with Pap test  
x3 = data$REDMETNO # Freq eating red meat during the past month: # of units  
x4 = data$WTIA_SA # Weight - Interim Annual 
dat = data.frame(y=y,x1=x1,x2=x2,x3=x3,x4=x4)
da = na.omit(dat)

yy=y=da$y
###################################################
#####################################zero-2 poisson
lf_poiss<- function(par) {
  ww2=1-(par[1])
  d2 <- (par[1] * (y == 6) + ww2 * dpois(y,par[2]))            
  -sum(log(d2))
}


 c_i=c(-1,0,0)
 u_i=rbind(c(-1,0),c(1,0),c(0,1))
 init=c(.1,.5)
 out2_poiss=constrOptim(init, lf_poiss, NULL, ui=u_i, ci=c_i)
###################################################
#####################################poisson
lf_poiss_w<- function(par) {
  d2 <- (dpois(y,par))            
  -sum(log(d2))
}

out2_poiss_w=optimize(lf_poiss_w, c(0, 1000))
###################################################
#####################################zero-2 discrete weibull
lf_dw<- function(par) {
  ww2=1-(par[1])
  d2 <- (par[1] * (y == 6) +  ww2 * ddweibull(y,par[2],par[3]))            
  -sum(log(d2))
}


 c_i=c(-1,0,0,0,-1)
 u_i=rbind(c(-1,0,0),c(1,0,0),c(0,1,0),c(0,0,1),c(0,-1,0))
 init=c(.1,.5,1)
 out2_dw=constrOptim(init, lf_dw, NULL, ui=u_i, ci=c_i)


#####################################discrete weibull
lf_dw_w<- function(par) {
  d2 <- ( ddweibull(y,par[1],par[2]))            
  -sum(log(d2))
}


 c_i=c(0,0,-1)
 u_i=rbind(c(1,0),c(0,1),c(-1,0))
 init=c(.5,1)
 out2_dw_w=constrOptim(init, lf_dw_w, NULL, ui=u_i, ci=c_i)


###################################################
#####################################zero-2 Geometric
lf_dgeo<- function(par) {
  ww2=1-(par[1])
  d2 <- (par[1] * (y == 6) + ww2 * dgeom(y,par[2]))            
  -sum(log(d2))
}


 c_i=c(-1,0,0,-1)
 u_i=rbind(c(-1,0),c(1,0),c(0,1),c(0,-1))
 init=c(.1,.1)
 out2_geo=constrOptim(init, lf_dgeo, NULL, ui=u_i, ci=c_i)

#####################################Geometric
lf_dgeo_w<- function(par) {
  d2 <- ( dgeom(y,par))            
  -sum(log(d2))
}

out2_geo_w=optimize(lf_dgeo_w, c(0, 1))
###################################################
#####################z_2 Discrete modified Weibull

dmw=function(x,the,c,q){
return(q^(x^the * c^x)-q^((x+1)^the * c^(x+1)) )
}

dF_dmw = function(x,the,c,q){
 return(1-q^((x^the) * (c^x)))
 }

r_dmw= function(n,par){
y=r_cdf(dF_dmw, n=n,  min = 0,max = 20, the=par[1],c=par[2],q=par[3])
return(floor(y))
}
#table((r_dmw(1000,c(2,2,.9))))

lf_dmw=function(par){
  ww2=1-(par[1])
  d2 <- (par[1] * (y == 6) + ww2 * dmw(y,par[2],par[3],par[4]))            
  -sum(log(d2))
}

 c_i=c(-1,0,0,0,1,-1)
 u_i=rbind(c(-1,0,0,0),c(1,0,0,0),c(0,1,0,0),c(0,0,0,1),c(0,0,1,0),c(0,0,0,-1))
 init=c(.5,1,1.01,.75)
 out2_dmw=constrOptim(init, lf_dmw, NULL, ui=u_i, ci=c_i)

#####################Discrete modified Weibull
lf_dmw_w=function(par){
  d2 <- (dmw(y,par[1],par[2],par[3]))            
  -sum(log(d2))
}

 c_i=c(0,0,1,-1)
 u_i=rbind(c(1,0,0),c(0,0,1),c(0,1,0),c(0,0,-1))
 init=c(1,1.01,.75)
 out2_dmw_w=constrOptim(init, lf_dmw_w, NULL, ui=u_i, ci=c_i)
################################################
################################################
################################################
################################################

out2_poiss_w$objective
out2_dw_w$value
out2_geo_w$objective
out2_dmw_w$value

out2_poiss$value
out2_dw$value
out2_geo$value
out2_dmw$value

######################################
n=length(y)
aic=bic=aicc=caic=c()
s=5 # Number of parameter

aic[1]=2*(s-4)+2*(out2_poiss_w$objective)
bic[1]=log(n)*(s-4)+2*(out2_poiss_w$objective)


aic[2]=2*(s-3)+2*(out2_dw_w$value)
bic[2]=log(n)*(s-3)+2*(out2_dw_w$value)


aic[3]=2*(s-4)+2*(out2_geo_w$objective)
bic[3]=log(n)*(s-4)+2*(out2_geo_w$objective)


aic[4]=2*(s-2)+2*(out2_dmw_w$value)
bic[4]=log(n)*(s-2)+2*(out2_dmw_w$value)

s=3
aic[5]=2*(s-1)+2*(out2_poiss$value)
bic[5]=log(n)*(s-1)+(s-1)*(out2_poiss$value)


aic[6]=2*(s)+2*(out2_dw$value)
bic[6]=log(n)*(s)+2*(out2_dw$value)


aic[7]=2*(s-1)+2*(out2_geo$value)
bic[7]=log(n)*(s-1)+2*(out2_geo$value)


aic[8]=2*(s+1)+2*(out2_dmw$value)
bic[8]=log(n)*(s+1)+2*(out2_dmw$value)


aic
bic


table(na.omit(y)>6)
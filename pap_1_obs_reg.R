rm(list=ls(all=TRUE))
library(EnvStats)
library(extraDistr)
library(VGAM)
library(rando)
library(data.table)

data = read.csv('d://Data/pap.csv')
y = data$PAP6YR1 # Number of Pap tests, last 6 years
x1 = data$R_MARITL #
x2 = data$HPVPAP # Marital Status 
x3 = data$REDMETNO # Freq eating red meat during the past month: # of units  
x4 = data$WTIA_SA # Weight - Interim Annual 
dat = data.frame(y=y,x1=x1,x2=x2,x3=x3,x4=x4)
da = na.omit(dat)

yy=y=da$y
x1=da$x1
x2=da$x2
x3=da$x3
x4=da$x4

j=jj=1
ind1=c()
for(i in 1:(length(y))){
 if(y[i]<=19) {
  ind1[jj]=j
  jj=jj+1
 }
 j=j+1
}
y=y[ind1]
x1=x1[ind1]
x2=x2[ind1]
x3=x3[ind1]
x4=x4[ind1]

###################################################
#####################################zero-2 poisson
lf_poiss<- function(par) {
 w0=exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)/
    (1+exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)) 
 ww2=1-(w0)
  d2 <- (w0 * (y == 6) + ww2 * dpois(y,par[6]))            
  -sum(log(d2))
}


 c_i=c(0)
 u_i=rbind(c(0,0,0,0,0,1))
 init=c(.01,.001,.001,.001,.0001,.1)
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
 w0=exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)/
    (1+exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)) 
 ww2=1-(w0)
  d2 <- (w0 * (y == 6) +  ww2 * ddweibull(y,par[6],par[7]))            
  -sum(log(d2))
}


 c_i=c(0,0,-1)
 u_i=rbind(c(0,0,0,0,0,1,0),c(0,0,0,0,0,0,1),c(0,0,0,0,0,-1,0))
 init=c(.01,.001,.001,.001,.0001,.5,1)
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
 w0=exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)/
    (1+exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)) 
 ww2=1-(w0)
  d2 <- (w0 * (y == 6) + ww2 * dnbinom(y,1,par[6]))            
  -sum(log(d2))
}


 c_i=c(0,-1)
 u_i=rbind(c(0,0,0,0,0,1),c(0,0,0,0,0,-1))
 init=c(.01,.001,.001,.001,.0001,.1)
 out2_geo=constrOptim(init, lf_dgeo, NULL, ui=u_i, ci=c_i)
 out2_geo$value
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
 w0=exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)/
    (1+exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)) 
 ww2=1-(w0)
  d2 <- (w0 * (y == 6) + ww2 * dmw(y,par[6],par[7],par[8]))            
  -sum(log(d2))
}

 c_i=c(0,0,1,-1)
 u_i=rbind(c(0,0,0,0,0,1,0,0),c(0,0,0,0,0,0,0,1),
c(0,0,0,0,0,0,1,0),c(0,0,0,0,0,0,0,-1))
 init=c(.01,.001,.001,.001,.0001,1,1.01,.75)
 out2_dmw=constrOptim(init, lf_dmw, NULL, ui=u_i, ci=c_i)

#####################Discrete modified Weibull
lf_dmw_w=function(par){
  d2 <- (dmw(y,par[1],par[2],par[3]))            
  -sum(log(d2))
}

 c_i=c(0,0,1,-1)
 u_i=rbind(c(1,0,0),c(0,0,1),c(0,1,0),c(0,0,-1))
 init=c(1,1.01,.75)
 out2_dmw_w=constrOptim(init, lf_dmw_w, NULL, ui=u_i, ci=c_i)################################################
################################################
################################################
###########################################Poiss
pmf_poiss <- function(par,x) {
 w0=exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)/
    (1+exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)) 
 ww2=1-(w0)
  d2 <- (w0 * (x == 6) + ww2 * dpois(x,par[6]))            
  d2
}

################################################
################################## Dw

pmf_dw <- function(par,y) {
 w0=exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)/
    (1+exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)) 
 ww2=1-(w0)
  d2 <- (w0 * (y == 6) +  ww2 * ddweibull(y,par[6],par[7]))            
  d2
}

################################################
################################## GEO

pmf_geo <- function(par,y) {
 w0=exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)/
    (1+exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)) 
 ww2=1-(w0)
  d2 <- (w0 * (y == 6) +  ww2 * dgeom(y,par[6]))            
  d2
}

################################################
################################## Dmw

pmf_dmw <- function(par,y) {
 w0=exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)/
    (1+exp(par[1]+par[2]*x1+par[3]*x2+par[4]*x3+par[5]*x4)) 
 ww2=1-(w0)
  d2 <- (w0 * (y == 6) +  ww2 * dmw(y,par[6],par[7],par[8]))            
  d2
}

################################################
################################################
z=y
nn=length(y)

m=length(table(z))
k=7
o=c(table(y)[1:(k)],sum(table(y)[(k+1):(m)]))
################################################
##########################################poisson
ee_poiss=c()
for(i in 0:(m-1)){
ee_poiss[i+1]=mean(pmf_poiss(out2_poiss$par,i)*nn)
}
e6=c(ee_poiss[1:k],sum(ee_poiss[(k+1):m]))
round(e6,2)

ABE6=sum(abs(o-e6))
round(ABE6,5)

KI6=sum(((o-e6)^2)/e6)
round(KI6,5)
##########################################poisson_w
ee_poiss_w=c()
for(i in 0:(m-1)){
ee_poiss_w[i+1]=mean(dpois(i,out2_poiss_w$minimum)*nn)
}
e6=c(ee_poiss_w[1:k],sum(ee_poiss_w[(k+1):m]))
round(e6,2)

ABE6=sum(abs(o-e6))
round(ABE6,5)

KI6=sum(((o-e6)^2)/e6)
round(KI6,5)
################################################
##########################################DW
ee_dw=c()
for(i in 0:(m-1)){
ee_dw[i+1]=mean(pmf_dw(out2_dw$par,i)*nn)
}
e6=c(ee_dw[1:k],sum(ee_dw[(k+1):m]))
round(e6,2)

ABE6=sum(abs(o-e6))
round(ABE6,5)

KI6=sum(((o-e6)^2)/e6)
round(KI6,5)
##########################################DW_w
ee_dw_w=c()
for(i in 0:(m-1)){
ee_dw_w[i+1]=mean(ddweibull(i,out2_dw_w$par[1],out2_dw_w$par[2])*nn)
}
e6=c(ee_dw_w[1:k],sum(ee_dw_w[(k+1):m]))
round(e6,2)

ABE6=sum(abs(o-e6))
round(ABE6,5)

KI6=sum(((o-e6)^2)/e6)
round(KI6,5)

################################################
##########################################GEO
ee_geo=c()
for(i in 0:(m-1)){
ee_geo[i+1]=mean(pmf_geo(out2_geo$par,i)*nn)
}
e6=c(ee_geo[1:k],sum(ee_geo[(k+1):m]))
round(e6,2)

ABE6=sum(abs(o-e6))
round(ABE6,5)

KI6=sum(((o-e6)^2)/e6)
round(KI6,5)
##########################################GEO_w
ee_geo_w=c()
for(i in 0:(m-1)){
ee_geo_w[i+1]=mean(dgeom(i,out2_geo_w$minimum)*nn)
}
e6=c(ee_geo_w[1:k],sum(ee_geo_w[(k+1):m]))
round(e6,2)

ABE6=sum(abs(o-e6))
round(ABE6,5)

KI6=sum(((o-e6)^2)/e6)
round(KI6,5)

################################################
##########################################DMW
ee_dmw=c()
for(i in 0:(m-1)){
ee_dmw[i+1]=mean(pmf_dmw(out2_dmw$par,i)*nn)
}
e6=c(ee_dmw[1:k],sum(ee_dmw[(k+1):m]))
round(e6,2)

ABE6=sum(abs(o-e6))
round(ABE6,5)

KI6=sum(((o-e6)^2)/e6)
round(KI6,5)

##########################################DMW_w
ee_dmw_w=c()
for(i in 0:(m-1)){
ee_dmw_w[i+1]=mean(dmw(i,out2_dmw_w$par[1],
out2_dmw_w$par[2],out2_dmw_w$par[3])*nn)
}
e6=c(ee_dmw_w[1:k],sum(ee_dmw_w[(k+1):m]))
round(e6,2)

ABE6=sum(abs(o-e6))
round(ABE6,5)

KI6=sum(((o-e6)^2)/e6)
round(KI6,5)
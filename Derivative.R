rm(list=ls(all=TRUE))
library(EnvStats)
library(extraDistr)
library(VGAM)
library(rando)
library("Deriv")
DD <- function(fname, x, order = 1) {
   if(order < 1) stop("'order' must be >= 1")
   if(order == 1) Deriv(fname,x)
   else DD(Deriv(fname, x), x, order - 1)
}

lf_dw<- function(par1,par2,par3,par4,par5,par6,par7){
 w0=exp(par1+par2*x1+par3*x2+par4*x3+par5*x4)/
    (1+exp(par1+par2*x1+par3*x2+par4*x3+par5*x4))
 ww2=1-(w0)
  d2 <- (w0 * (y == 6) +  ww2 * (par6^y^par7-par6^(y+1)^par7))            
  -sum(log(d2))
}

##############################################beta0
x=c("par1")#vector of unknown parameters
DD(lf_dw,x,2)

dd_beta0 = function (par1, par2, par3, par4, par5, par6, par7) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2 + par4 * x3 + par5 * 
        x4)
    .e2 <- 1 + .e1
    .e3 <- 1 - .e1/.e2
    .e4 <- y == 6
    .e5 <- ddweibull(y, par[6], par[7])
    .e6 <- .e1 * .e4
    .e7 <- (.e3 * .e5 + .e6/.e2) * .e2
    -sum((.e3/.e7 - .e6/.e7^2) * (.e4 - .e5) * .e3 * .e1)
}

##############################################beta1
x=c("par2")#vector of unknown parameters
DD(lf_dw,x,2)

dd_beta1 = function (par1, par2, par3, par4, par5, par6, par7) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2 + par4 * x3 + par5 * 
        x4)
    .e2 <- 1 + .e1
    .e3 <- 1 - .e1/.e2
    .e4 <- y == 6
    .e5 <- ddweibull(y, par[6], par[7])
    .e6 <- .e1 * .e4
    .e7 <- (.e3 * .e5 + .e6/.e2) * .e2
    -sum(x1^2 * (.e3/.e7 - .e6/.e7^2) * (.e4 - .e5) * .e3 * .e1)
}
##############################################beta2
x=c("par3")#vector of unknown parameters
DD(lf_dw,x,2)

dd_beta2 = function (par1, par2, par3, par4, par5, par6, par7) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2 + par4 * x3 + par5 * 
        x4)
    .e2 <- 1 + .e1
    .e3 <- 1 - .e1/.e2
    .e4 <- y == 6
    .e5 <- ddweibull(y, par[6], par[7])
    .e6 <- .e1 * .e4
    .e7 <- (.e3 * .e5 + .e6/.e2) * .e2
    -sum(x2^2 * (.e3/.e7 - .e6/.e7^2) * (.e4 - .e5) * .e3 * .e1)
}
##############################################beta3
x=c("par4")#vector of unknown parameters
DD(lf_dw,x,2)

dd_beta3 = function (par1, par2, par3, par4, par5, par6, par7) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2 + par4 * x3 + par5 * 
        x4)
    .e2 <- 1 + .e1
    .e3 <- 1 - .e1/.e2
    .e4 <- y == 6
    .e5 <- ddweibull(y, par[6], par[7])
    .e6 <- .e1 * .e4
    .e7 <- (.e3 * .e5 + .e6/.e2) * .e2
    -sum(x3^2 * (.e3/.e7 - .e6/.e7^2) * (.e4 - .e5) * .e3 * .e1)
}
##############################################beta4
x=c("par5")#vector of unknown parameters
DD(lf_dw,x,2)

dd_beta4 = function (par1, par2, par3, par4, par5, par6, par7) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2 + par4 * x3 + par5 * 
        x4)
    .e2 <- 1 + .e1
    .e3 <- 1 - .e1/.e2
    .e4 <- y == 6
    .e5 <- ddweibull(y, par[6], par[7])
    .e6 <- .e1 * .e4
    .e7 <- (.e3 * .e5 + .e6/.e2) * .e2
    -sum(x4^2 * (.e3/.e7 - .e6/.e7^2) * (.e4 - .e5) * .e3 * .e1)
}

##############################################alpha
x=c("par6")#vector of unknown parameters
DD(lf_dw,x,2)

dd_alpha = function (par1, par2, par3, par4, par5, par6, par7) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2 + par4 * x3 + par5 * 
        x4)
    .e2 <- (1 + y)^par7
    .e3 <- y^par7
    .e4 <- 1 + .e1
    .e5 <- 1 - .e1/.e4
    .e7 <- .e5 * (par6^.e3 - par6^.e2) + .e1 * (y == 6)/.e4
    .e8 <- .e2 - 1
    .e9 <- .e3 - 1
    -sum(.e5 * (par6^(.e3 - 2) * .e3 * .e9 - (.e5 * (par6^.e9 * 
        .e3 - par6^.e8 * .e2)^2/.e7 + par6^(.e2 - 2) * .e8 * 
        .e2))/.e7)
}

##############################################gamma
x=c("par7")#vector of unknown parameters
DD(lf_dw,x,2)

dd_gamma = function (par1, par2, par3, par4, par5, par6, par7) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2 + par4 * x3 + par5 * 
        x4)
    .e2 <- 1 + y
    .e3 <- .e2^par7
    .e4 <- y^par7
    .e5 <- 1 + .e1
    .e6 <- par6^.e3
    .e7 <- par6^.e4
    .e8 <- 1 - .e1/.e5
    .e9 <- log(par6)
    .e11 <- .e8 * (.e7 - .e6) + .e1 * (y == 6)/.e5
    .e12 <- .e8 * .e9
    .e13 <- 2 * par7
    .e14 <- log(y+1)
    .e15 <- log1p(y)
    .e16 <- .e6 * .e3
    .e17 <- .e7 * .e4
    -sum(.e12 * (.e14^2 * (.e7 * y^.e13 * .e9 + .e17) - (.e12 * 
        (.e17 * .e14 - .e16 * .e15)^2/.e11 + .e15^2 * (.e6 * 
        .e2^.e13 * .e9 + .e16)))/.e11)
}


###################################################
###################################################
###################################################
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
x4 = log(x4)
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
par = out2_dw$par

#SE(beta0)
round(1/dd_beta0(par[1],par[2],par[3],par[4],par[5],par[6],par[7]),4)

#SE(beta1)
round(1/dd_beta1(par[1],par[2],par[3],par[4],par[5],par[6],par[7]),4)

#SE(beta2)
round(1/dd_beta2(par[1],par[2],par[3],par[4],par[5],par[6],par[7]),4)

#SE(beta3)
round(1/dd_beta3(par[1],par[2],par[3],par[4],par[5],par[6],par[7]),4)

#SE(beta4)
round(1/dd_beta4(par[1],par[2],par[3],par[4],par[5],par[6],par[7]),4)

#SE(alpha)
round(1/dd_alpha(par[1],par[2],par[3],par[4],par[5],par[6],par[7]),4)

#SE(gamma)
round(1/dd_gamma(par[1],par[2],par[3],par[4],par[5],par[6],par[7]),4)


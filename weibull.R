rm(list=ls(all=TRUE))

library(survival)

data(flchain)
data = flchain
#hist(log(flchain$futime),breaks=100)
# creating training dataset
mse_model = mse_lm = pre_model = pre_lm = c()
l_model = l_lm = aic_model = bic_model = aic_lm = bic_lm = c()
for(iii in 1:1){
sample <- sample(c(TRUE,FALSE), nrow(data),  
                 replace=TRUE, prob=c(0.7,0.3)) 
train  <- data[sample, ] 
  
# creating testing dataset 
test  <- data[!sample, ] 

y=train$futime
y=log(y+2)
x1=train$lambda
x2=train$age
med=median(y)
y_te=test$futime
y_te=log(y_te+2)
x1_te=test$lambda
x2_te=test$age
######################################################
##################### Gamma

lf_norm=function(par){
  w0 = (exp(par[1]+par[2]*x1+par[3]*x2))/
(exp(par[1]+par[2]*x1+par[3]*x2)+1)
  etww=1-w0
  d <- (w0 * ((y>8)&(y<9)) + etww * dweibull(y,par[4],par[5]))            
  -sum(log(d))
}
 c_i=c(0,0)
 u_i=rbind(c(0,0,0,1,0),c(0,0,0,0,1))
 init=c(.01,.01,.01,10,10)
 out=constrOptim(init, lf_norm, NULL, ui=u_i, ci=c_i)
#out$value

lf_lm=function(par){
  d <- (y-par[1]-par[2]*x1-par[3]*x2)^2    
  sum(d)
}
 c_i=c(-100)
 u_i=rbind(c(0,0,1))
 init=c(.01,.01,.01)
 out_lm=constrOptim(init, lf_lm, NULL, ui=u_i, ci=c_i)

pred = function(x1,x2){
  par = out$par
  w0 = (exp(par[1]+par[2]*x1+par[3]*x2))/
   (exp(par[1]+par[2]*x1+par[3]*x2)+1)
  p=c(w0,(1-w0))
  out=c()
  for(i in 1:100){
  ran=sample(c(0,1), size=1, prob=p)
    if(ran==0) out[i]=med
    else if(ran==1) out[i]= rweibull(1,par[4],par[5])
 }
   mean(out)
}
predictt = function(x1,x2){
n=length(x1)
 predd=c()
 for(i in 1:n){
  predd[i] = pred(x1[i],x2[i])

 }
 predd
}
predd = predictt(x1_te,x2_te)
mse_model[iii] = mean((y_te-predd)^2)

pa = out_lm$par
fit = pa[1] + pa[2] * x1_te + pa[3] * x2_te
fit_1 = pa[1] + pa[2] * x1 + pa[3] * x2
mse_lm[iii] = sum((fit - y_te)^2)

m1 = table((round(predd))==(round(y_te)))
n = length(y_te)
pre_model[iii] = m1[2]/n
te = data.frame(x1=x1_te,x2=x2_te)
m2 = table((round(fit))==(round(y_te)))
pre_lm[iii] = m2[2]/n
k1=5
k2=3
l1 = out$value

l2 = -sum(log(dnorm(fit_1-y,0,sqrt(sum((fit_1-y)^2)/(length(y)-2)))))
l_model[iii] = out$value
l_lm[iii] = out_lm$value
aic_model[iii] = 2*k1+2*l1
bic_model[iii] = k1*log(n)+2*l1
aic_lm[iii] = 2*k2+2*l2
bic_lm[iii] = k2*log(n)+2*l2
print(iii)
}
mean(mse_model)
mean(pre_model)
mean(l_model)
mean(aic_model)
mean(bic_model)

mean(mse_lm)
mean(pre_lm)
mean(l_lm)
mean(aic_lm)
mean(bic_lm)



library("Deriv")
DD <- function(fname, x, order = 1) {
   if(order < 1) stop("'order' must be >= 1")
   if(order == 1) Deriv(fname,x)
   else DD(Deriv(fname, x), x, order - 1)
}

lf_norm=function(par1,par2,par3,par4,par5){
 w0=exp(par1+par2*x1+par3*x2)/
    (1+exp(par1+par2*x1+par3*x2))
 ww2=1-(w0)
  d <- (w0 * ((y>8)&(y<9)) + ww2 * ((par4/par5)*(y/par5)^(par4-1)*
exp(-(y/par5)^par4)))            
  -sum(log(d))
}

##############################################beta0
x=c("par1")#vector of unknown parameters
DD(lf_norm,x,2)

dd_beta0 = function (par1, par2, par3, par4, par5) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2)
    .e2 <- 1 + .e1
    .e3 <- y/par5
    .e4 <- 1 - .e1/.e2
    .e7 <- y > 8 & y < 9
    .e9 <- .e3^(par4 - 1)
    .e10 <- exp(-.e3^par4)
    .e11 <- .e1 * .e7
    .e12 <- .e2 * (.e11/.e2 + par4 * .e4 * .e10 * .e9/par5)
    -sum((.e4/.e12 - .e11/.e12^2) * (.e7 - par4 * .e10 * .e9/par5) * 
        .e4 * .e1)
}



##############################################beta1
x=c("par2")#vector of unknown parameters
DD(lf_norm,x,2)

dd_beta1 = function (par1, par2, par3, par4, par5) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2)
    .e2 <- 1 + .e1
    .e3 <- y/par5
    .e4 <- 1 - .e1/.e2
    .e7 <- y > 8 & y < 9
    .e9 <- .e3^(par4 - 1)
    .e10 <- exp(-.e3^par4)
    .e11 <- .e1 * .e7
    .e12 <- .e2 * (.e11/.e2 + par4 * .e4 * .e10 * .e9/par5)
    -sum(x1^2 * (.e4/.e12 - .e11/.e12^2) * (.e7 - par4 * .e10 * 
        .e9/par5) * .e4 * .e1)
}


##############################################beta2
x=c("par3")#vector of unknown parameters
DD(lf_norm,x,2)

dd_beta2 = function (par1, par2, par3, par4, par5) 
{
    .e1 <- exp(par1 + par2 * x1 + par3 * x2)
    .e2 <- 1 + .e1
    .e3 <- y/par5
    .e4 <- 1 - .e1/.e2
    .e7 <- y > 8 & y < 9
    .e9 <- .e3^(par4 - 1)
    .e10 <- exp(-.e3^par4)
    .e11 <- .e1 * .e7
    .e12 <- .e2 * (.e11/.e2 + par4 * .e4 * .e10 * .e9/par5)
    -sum(x2^2 * (.e4/.e12 - .e11/.e12^2) * (.e7 - par4 * .e10 * 
        .e9/par5) * .e4 * .e1)
}



##############################################alpha
x=c("par4")#vector of unknown parameters
DD(lf_norm,x,2)

dd_alpha = function (par1, par2, par3, par4, par5) 
{
    .e1 <- y/par5
    .e2 <- exp(par1 + par2 * x1 + par3 * x2)
    .e3 <- .e1^(par4 - 1)
    .e4 <- .e1^par4
    .e5 <- 1 + .e2
    .e8 <- log(y) - log(par5)
    .e10 <- .e1^(2 * par4 - 1)
    .e11 <- 1 - .e2/.e5
    .e12 <- exp(-.e4)
    .e14 <- .e3 + par4 * (.e3 - .e10) * .e8
    .e18 <- par4 * .e8
    .e19 <- par5 * (.e2 * (y > 8 & y < 9)/.e5 + par4 * .e11 * 
        .e12 * .e3/par5)
    -sum(((2 * .e3 + par4 * (.e3 - 2 * .e10) * .e8 - (.e14 * 
        .e4 + .e10)) * .e8/.e19 - ((1 - .e18 * .e4) * .e3 + .e18 * 
        .e3) * .e14 * .e11 * .e12/.e19^2) * .e11 * .e12)
}



##############################################gamma
x=c("par5")#vector of unknown parameters
DD(lf_norm,x,2)

dd_gamma = function (par1, par2, par3, par4, par5) 
{
    .e1 <- y/par5
    .e2 <- exp(par1 + par2 * x1 + par3 * x2)
    .e3 <- par4 - 1
    .e4 <- 1 + .e2
    .e5 <- .e1^.e3
    .e8 <- 2 * .e3
    .e10 <- par4 - 2
    .e12 <- par4 * (1 - .e2/.e4) * exp(-.e1^par4)
    .e13 <- .e1^.e8
    .e14 <- .e1^.e10
    .e18 <- .e2 * (y > 8 & y < 9)/.e4 + .e12 * .e5/par5
    .e23 <- y * (par4 * .e13 - .e3 * .e14)/par5 - .e5
    -sum(.e12 * (y * ((2 * .e14 - y * (2 * (par4 * .e1^(.e8 - 
        1)) - .e10 * .e1^(par4 - 3))/par5) * .e3 + par4 * (.e23 * 
        .e5 - .e13))/(par5^4 * .e18) - par5 * (2 * .e18 + .e12 * 
        .e23/par5) * .e23/(par5^2 * .e18)^2))
}

par = out$par
#SE(beta0)
round(1/dd_beta0(par[1],par[2],par[3],par[4],par[5]),4)

#SE(beta1)
round(1/dd_beta1(par[1],par[2],par[3],par[4],par[5]),4)

#SE(beta2)
round(1/dd_beta2(par[1],par[2],par[3],par[4],par[5]),4)

#SE(alpha)
round(1/dd_alpha(par[1],par[2],par[3],par[4],par[5]),4)

#SE(gamma)
round(1/dd_gamma(par[1],par[2],par[3],par[4],par[5]),4)




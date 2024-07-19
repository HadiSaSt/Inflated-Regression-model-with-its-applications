rm(list=ls(all=TRUE))

library(survival)

data(flchain)
data = flchain
#hist(log(flchain$futime),breaks=100)
# creating training dataset
mse_model = mse_lm = pre_model = pre_lm = c()
l_model = l_lm = aic_model = bic_model = aic_lm = bic_lm = c()
for(iii in 1:1000){
sample <- sample(c(TRUE,FALSE), nrow(data),  
                 replace=TRUE, prob=c(0.7,0.3)) 
train  <- data[sample, ] 
  
# creating testing dataset 
test  <- data[!sample, ] 

y=train$futime
y=log(y+2)
x1=train$death
x2=train$age
med=median(y)
y_te=test$futime
y_te=log(y_te+2)
x1_te=test$death
x2_te=test$age
######################################################
##################### Gamma

lf_norm=function(par){
  w0 = (exp(par[1]+par[2]*x1+par[3]*x2))/
(exp(par[1]+par[2]*x1+par[3]*x2)+1)
  etww=1-w0
  d <- (w0 * ((y>8)&(y<9)) + etww * dgamma(y,par[4],par[5]))            
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
    else if(ran==1) out[i]= rgamma(1,par[4],par[5])
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
k1=4
k2=3
l1 = out$value

l2 = -sum(log(dnorm(fit_1-y,0,sqrt(sum((fit_1-y)^2)/(length(y)-2)))))
l_model[iii] = out$value
l_lm[iii] = out_lm$value
aic_model[iii] = 2*k1+2*l1
bic_model[iii] = k1*log(n)+2*l1
aic_lm[iii] = 2*k2+2*l2
bic_lm[iii] = k2*log(n)+2*l2
#print(iii)
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
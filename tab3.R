rm(list=ls(all=TRUE))

library(extraDistr)

m=10000
n=1000

rzkip <- function(n,par) {
  p=c(par[1],par[2],(1-(par[1]+par[2])))
  outt=rep(0,n)
  for(i in 1:n)
   {
    ran=sample(c(0,1,2), size=1, prob=p)
    if(ran==0) outt[i]=0
    else if(ran==1) outt[i]=2
    else if(ran==2) outt[i]=rpois(1,exp(par[3]+par[4]*x1[i]+par[5]*x2[i]))
   }
return(outt)
}  

out_pw_b0=out_pw_b1=out_pw_b2=c()
out_p_inf_b0=out_p_inf_b1=out_p_inf_b2=c()
out_pcom_inf_b0=out_pcom_inf_b1=out_pcom_inf_b2=c()

out_pw_w0=out_pw_w2=c()
out_p_inf_w0=out_p_inf_w2=c()
out_pcom_inf_w0=out_pcom_inf_w2=c()

out_pw=c()
out_p_inf=c()
out_pcom_inf=c()

for(ii in 1:m){                                                                                                                                                                                                                                                                                                                                                                                                               
x1=rnorm(n)
x2=x2=rbinom(n,1,.5)
par=c(.05,.05,-.5,.5,.5)
y=rzkip(n,par)
table(y)
###################################################
#####################################zero-2 poisson
lf_p<- function(par) {
  ww2=1-(par[1]+par[2])
  la2=exp(par[3]+par[4]*x1+par[5]*x2)
  d0 <- (par[1] * (y == 0) + par[2] * (y == 2) +
  ww2 * dpois(y,la2))            
  -sum(log(d0))
}


 c_i=c(-1,0,0)
 u_i=rbind(c(-1,-1,0,0,0),c(1,0,0,0,0),c(0,1,0,0,0))
 init=c(.1,.1,.1,.1,.1)

 out2_p=constrOptim(init, lf_p, NULL, ui=u_i, ci=c_i)
################################################################### EM

lf_p_comp=function(par){
  ww2=1-(par[1]+par[2])
  d4 <- (par[1]/(par[1]+ww2*dpois(y,exp(par[3]+par[4]*x1+par[5]*x2)))*log(par[1]/ww2)*(y==0)+
         par[2]/(par[2]+ww2*dpois(y,exp(par[3]+par[4]*x1+par[5]*x2)))*log(par[2]/ww2)*(y==2)-
         log(1+(par[1]+par[2])/ww2)+ww2*log(dpois(y,exp(par[3]+par[4]*x1+par[5]*x2)))
)            
  -sum(d4)
}

 c_i=c(-1,0,0)
 u_i=rbind(c(-1,-1,0,0,0),c(1,0,0,0,0),c(0,1,0,0,0))
 init=c(.1,.1,.1,.1,.1)
 out2_p_comp=constrOptim(init, lf_p_comp, NULL, ui=u_i, ci=c_i)
################################################################### Poiss
lf_p_w<- function(par) {
  la2=exp(par[1]+par[2]*x1+par[3]*x2)
  d0 <- dpois(y,la2)           
  -sum(log(d0))
}


 c_i=c(-100)
 u_i=rbind(c(-1,0,0))
 init=c(.1,.1,.1)

 out2_p_w=constrOptim(init, lf_p_w, NULL, ui=u_i, ci=c_i)



pp_w=out2_p_w$par
pp=out2_p$par
pp_com=out2_p_comp$par

out_pw_b0[ii]=pp_w[1]
out_pw_b1[ii]=pp_w[2]
out_pw_b2[ii]=pp_w[3]
out_p_inf_b0[ii]=pp[3]
out_p_inf_b1[ii]=pp[4]
out_p_inf_b2[ii]=pp[5]
out_pcom_inf_b0[ii]=pp_com[3]
out_pcom_inf_b1[ii]=pp_com[4]
out_pcom_inf_b2[ii]=pp_com[5]

out_p_inf_w0[ii]=pp[1]
out_p_inf_w2[ii]=pp[2]
out_pcom_inf_w0[ii]=pp_com[1]
out_pcom_inf_w2[ii]=pp_com[2]

out_pw[ii]=out2_p_w$value
out_p_inf[ii]=out2_p$value
out_pcom_inf[ii]=out2_p_comp$value

print(ii)
}

round(mean(out_pw_b0),4)+.5
round(mean(out_pw_b1),4)-.5
round(mean(out_pw_b2),4)-.5
round(mean(out_p_inf_b0),4)+.5
round(mean(out_p_inf_b1),4)-.5
round(mean(out_p_inf_b2),4)-.5
round(mean(out_pcom_inf_b0),4)+.5
round(mean(out_pcom_inf_b1),4)-.5
round(mean(out_pcom_inf_b2),4)-.5

round(mean(out_p_inf_w0),4)-.05
round(mean(out_p_inf_w2),4)-.05
round(mean(out_pcom_inf_w0),4)-.05
round(mean(out_pcom_inf_w2),4)-.05

round(mean(out_pw),4)
round(mean(out_p_inf),4)
round(mean(out_pcom_inf),4)

round(mean((out_pw_b0+.5)^2),4)
round(mean((out_pw_b1-.5)^2),4)
round(mean((out_pw_b2-.5)^2),4)
round(mean((out_p_inf_b0+.5)^2),4)
round(mean((out_p_inf_b1-.5)^2),4)
round(mean((out_p_inf_b2-.5)^2),4)
round(mean((out_pcom_inf_b0+.5)^2),4)
round(mean((out_pcom_inf_b1-.5)^2),4)
round(mean((out_pcom_inf_b2-.5)^2),4)

round(mean((out_p_inf_w0-.05)^2),4)
round(mean((out_p_inf_w2-.05)^2),4)
round(mean((out_pcom_inf_w0-.05)^2),4)
round(mean((out_pcom_inf_w2-.05)^2),4)

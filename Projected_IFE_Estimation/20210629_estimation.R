rm(list = ls())

setwd("~/Dropbox/InteractivefixedEffect/2020_Projection_IFE/Code/Application")

library(splines)
library(lmtest)
library(sandwich)
library(xtable)

#Load data
data = read.csv("data_merged.csv")
data = subset(data,year>1990)
data = subset(data,country=="Austria"|country=="Australia"|country=="Belgium"|country=="Canada"|country=="Chile"
       |country=="Colombia"|country=="Denmark"|country=="Finland"|country=="France"|country=="Germany"
       |country=="Greece"|country=="Hungary"|country=="Iceland"|country=="Ireland"|country=="Israel"
       |country=="Italy"|country=="Japan"|country=="Korea"|country=="Luxembourg"|country=="Mexico"
       |country=="Netherlands"|country=="New Zealand"|country=="Norway"|country=="Poland"|country=="Portugal"
       |country=="Spain"|country=="Sweden"|country=="Switzerland"|country=="Turkey"|country=="United Kingdom"
       |country=="United States")

#Pooled OLS with White errors
y = data$gdp.growth
x = as.matrix(data[,c(7:13,16)])
#All restricted
#x = x[,c(1,2,5,6)]
#OECD resricted
#x = x[,c(1,3,5,6)]

fit1 = lm(y~x)
coeftest(fit1, vcov = vcovHC(fit1, type = "HC0"))

#P-IFE estimation
N = length(unique(data$country))
T = length(unique(data$year))
Q = ncol(x)

xbar = aggregate(data[,c(7:13,16)], list(data$country), mean)[,-1]
#All restricted
#xbar = aggregate(data[,c(7,8,11,12)], list(data$country), mean)[,-1]
#OECD restricted
#xbar = aggregate(data[,c(7,9,11,12)], list(data$country), mean)[,-1]

J = floor(2*sqrt(N))
designX = array(0,dim =c(N,J,Q))
for (q in 1: Q){
  designX[,,q] = bs(xbar[,q], intercept = T, degree=J-1)
}
phi = matrix(0,N,J*Q)
for (i in 1:N){
  phi[i,]   = as.vector(designX[i,,])
}
#J=2
#phi = cbind(as.matrix(xbar),as.matrix(xbar)^2)

a = solve(t(phi) %*% phi+ 0.00000001*diag(J*Q))
p = phi %*% a %*% t(phi)
m = diag(N) - p

Phi = phi[rep(1:N,times=T),]
A = solve(t(Phi) %*% Phi+ 0.00000001*diag(J*Q))

P = Phi %*% A %*% t(Phi)
M = diag(N*T) - P

fit2 =lm(M%*%y~M%*%x)
coeftest(fit2, vcov = vcovHC(fit2))

ytilde = matrix(fit2$residuals,nrow=N,ncol=T,byrow=T)

#Pure factor model (PPCA)
#ytilde = matrix(y,nrow=N,ncol=T,byrow=T)

#eigen(t(matrix(y,nrow=N,byrow=T))%*%matrix(y,nrow=N,byrow=T))$values[-T]/eigen(t(matrix(y,nrow=N,byrow=T))%*%matrix(y,nrow=N,byrow=T))$values[-1]
eigen(t(ytilde)%*%p%*%ytilde)$values[-T]/eigen(t(ytilde)%*%p%*%ytilde)$values[-1]
K=4

pdf(file = "factor_oecd.pdf", width = 7, height = 4, family = "Helvetica") # defaults to 7 x 7 inches
par(mar=c(2,2,1,1))
plot(1:29,eigen(t(ytilde)%*%p%*%ytilde)$values,cex=ifelse((1:29)>K,2,2),pch=16,type="b")
dev.off()

Fhat = eigen(t(ytilde)%*%p%*%ytilde)$vectors[,1:K]*sqrt(T)
Ghat = 1/T * p%*%ytilde%*%Fhat
lambdahat = 1/T * ytilde%*%Fhat

bhat = 1/T * solve(t(phi)%*%phi+0.0000000001*diag(J*Q)) %*% t(phi) %*% ytilde %*% Fhat

rownames(Ghat)=unique(data$country)
rownames(lambdahat)=unique(data$country)

#Plot OECD
pdf(file = "F_oecd.pdf", width = 7, height = 3, family = "Helvetica")
par(mar=c(2,2,2,2),mfrow=c(1,2))
plot(Fhat[,1]~data$year[1:T],type="l")
plot(Fhat[,2]~data$year[1:T],type="l")
dev.off()
#plot(Fhat[,3]~data$year[1:T],type="l")
#plot(Fhat[,4]~data$year[1:T],type="l")
#abline(v="1991",col="red")

#Plot all
pdf(file = "F_all.pdf", width = 4, height = 3, family = "Helvetica")
par(mar=c(2,2,1,1))
plot(Fhat[,1]~data$year[1:T],type="l",ylim=c(-2.2,3.5))
lines(Fhat[,2]~data$year[1:T],col="red")
lines(Fhat[,3]~data$year[1:T],col="blue")
dev.off()

#Average growth rates per year
pdf(file = "growthseries.pdf", width = 7, height = 4, family = "Helvetica") # defaults to 7 x 7 inches
par(mar=c(2,2,1,1))
plot(by(data$gdp.growth,data$year,mean)~data$year[1:T],type="l",xlab="",ylab="",ylim=c(-0.18,0.21))
lines(by(data$gdp.growth,data$year,quantile,probs=c(0.05))~data$year[1:T],col="black",lty="dashed")
lines(by(data$gdp.growth,data$year,quantile,probs=c(0.95))~data$year[1:T],col="black",lty="dashed")
dev.off()

#Descriptive statistics
table = matrix(0,9,4)
table[1,] = 100*c(mean(data$gdp.growth),median(data$gdp.growth),min(data$gdp.growth),max(data$gdp.growth))
table[2,] = c(mean(data$young),median(data$young),min(data$young),max(data$young))
table[3,] = c(mean(data$fert),median(data$fert),min(data$fert),max(data$fert))
table[4,] = c(mean(data$life),median(data$life),min(data$life),max(data$life))
table[5,] = 100*c(mean(data$pop.growth),median(data$pop.growth),min(data$pop.growth),max(data$pop.growth))
table[6,] = c(mean(data$pl_i),median(data$pl_i),min(data$pl_i),max(data$pl_i))
table[7,] = c(mean(data$csh_c),median(data$csh_c),min(data$csh_c),max(data$csh_c))
table[8,] = c(mean(data$csh_g),median(data$csh_g),min(data$csh_g),max(data$csh_g))
table[9,] = c(mean(data$csh_i),median(data$csh_i),min(data$csh_i),max(data$csh_i))
xtable(table,digits=2)

#Regression output
xtable(rbind(coeftest(fit2,vcov=vcovHC(fit2))[-1,1],coeftest(fit2,vcov=vcovHC(fit2))[-1,3]),digits=4)

#Ghat, Lambdahat
max(abs(Ghat))
max(abs(lambdahat-Ghat))
sqrt(sum((Ghat)^2))
sqrt(sum((lambdahat-Ghat)^2))

######
#betahat = solve(t(x)%*%kronecker(M,matrix(1,T,T))%*%x)%*%t(x)%*%kronecker(M,matrix(1,T,T))%*%y
#fit2 = lm(kronecker(M,matrix(1,T,T))%*%y~kronecker(M,matrix(1,T,T))%*%x)
#coeftest(fit2, vcov = vcovHC(fit2, type = "HC0"))


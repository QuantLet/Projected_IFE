rm(list = ls())

setwd("~/Dropbox/InteractivefixedEffect/2020_Projection_IFE/Code/Simulation")
library(splines)
library(mgcv)
library(gsynth)
library(tsDyn)
library(xtable)

##simulation function

simulation = function(S, N, T, Q, K){
  
  #1.1 result matrix
  result2 = result3 = vector(length=4)
  
  #1.2 fix regression coefficients
  beta = c(2,1,-1)

  #1.3 simulate regressors
  set.seed(1)
  xbar = matrix(rnorm(N*Q,1,0.5),N,Q)
  set.seed(2)
  pi = array(rnorm(N*T*Q,sd=0.5),dim=c(N,T,Q))
  x = array(0,dim=c(N,T,Q))
  for (t in 1:T){
    x[,t,] = pi[,t,]+xbar
  }

  #1.4 generate functions g_kq
  set.seed(3)
  a = matrix(runif(K*Q,-1,1),K,Q)
  set.seed(4)
  b = matrix(runif(K*Q,-1,1),K,Q)

  g = matrix(0,N,K)
  g[,1] = a[1,1]*xbar[,1]^2+b[1,2]*xbar[,2]
  g[,2] = a[2,2]*xbar[,2]^2+b[2,3]*xbar[,3]
  g[,3] = a[3,3]*xbar[,3]^2+b[3,1]*xbar[,1]
  
  #1.5 get basis functions, projection matrix, residual maker matrix
  #J = floor(sqrt(N))
  J = 3
  designX = array(0,dim =c(N,J,Q))
  for (q in 1: Q){
    designX[,,q] = bs(xbar[,q], intercept = T, degree=J-1)
  }
  phi = matrix(0,N,J*Q)
  for (i in 1:N){
    phi[i,]   = as.vector(designX[i,,])
  }
  A = solve(t(phi) %*% phi + 0.00001*diag(J*Q))
  
  phi = cbind(xbar,xbar^2)
  A = solve(t(phi) %*% phi)
  
  P = phi %*% A %*% t(phi)
  M = diag(N) - P

  for (s in 1:S){
    print(paste("s = ",s))
  #1.6 simulate factors
  #F = VAR.sim(diag(0.5,K,K),T,include="none",varcov=diag(0.75,K))
  F = VAR.sim(diag(0,K,K),T,include="none",varcov=diag(1,K))

    #2.1 simulate error term u_it, idiosyncratic effects gamma_i
    u = matrix(rnorm(N*T,0,1),N,T)
    gamma = matrix(0,N,K)
    #u = matrix(rt(N*T,df=10),N,T)
    #gamma = matrix(rt(N*K,df=10)/20,N,K)
  
    #2.2 generate y
    y = matrix(0,N,T)
    for(i in 1:N){
      y[i,] = x[i,,]%*%beta + as.matrix(F,T,K)%*% (g[i,]+gamma[i]) + u[i,]
    }

    #2.3 Projection-based IFE estimator
    sumaa   = matrix(0, Q, Q)
    sumbb   = matrix(0, Q, 1)
    for(t in 1:T){
      aa =  (t(x[,t,])%*%M%*%x[,t,])
      bb =  t(x[,t,])%*%M%*%y[,t]
      sumaa = aa+ sumaa
      sumbb = bb+ sumbb
    } 
    betahat = solve(sumaa)%*%sumbb
   
    #2.4 Estimation of interactive fixed effect components
    ytilde = matrix(0,N,T)
    for (t in 1:T){
      ytilde[,t] = y[,t]-x[,t,]%*%betahat  
    }
    Fhat = eigen(t(ytilde)%*%P%*%ytilde)$vectors[,1:K]*sqrt(T)
    Ghat = 1/T * P%*%ytilde%*%Fhat
    lambdahat = 1/T * ytilde%*%Fhat
    
    Fhat_PCA = eigen(t(ytilde)%*%ytilde)$vectors[,1:K]*sqrt(T)
    lambdahat_PCA = 1/T * ytilde%*%Fhat_PCA
    
    F0 = sqrt(T)*eigen(F%*%t(g)%*%g%*%t(F))$vector[,1:K]
    G0 = 1/T*g%*%t(F)%*%F0
    
    Fmax = max(abs(F0-Fhat))
    FFrob = sqrt(sum((F0-Fhat)^2))/sqrt(T)
    Gmax = max(abs(G0-Ghat))
    GFrob = sqrt(sum((G0-Ghat)^2))/sqrt(N)
    
    result2[1] = result2[1] + Fmax
    result2[2] = result2[2] + FFrob
    result2[3] = result2[3] + Gmax 
    result2[4] = result2[4] + GFrob
    
    Fmax_PCA = max(abs(F0-Fhat_PCA))
    FFrob_PCA = sqrt(sum((F0-Fhat_PCA)^2))/sqrt(T)
    lambdamax_PCA = max(abs(G0+gamma-lambdahat_PCA))
    lambdaFrob_PCA = sqrt(sum((G0+gamma-lambdahat_PCA)^2))/sqrt(N)
    
    result3[1] = result3[1] + Fmax_PCA
    result3[2] = result3[2] + FFrob_PCA
    result3[3] = result3[3] + lambdamax_PCA
    result3[4] = result3[4] + lambdaFrob_PCA
  }
  result2 = result2/S
  result3 = result3/S
  
  result = rbind(result2,result3) 
  return(result)
}


##1
#T=10 fixed
results = array(0,dim=c(11,2,4))

results[1,,] = simulation(1000,25,10,3,3)
results[2,,] = simulation(1000,50,10,3,3)
results[3,,] = simulation(1000,100,10,3,3)
results[4,,] = simulation(1000,150,10,3,3)
results[5,,] = simulation(1000,200,10,3,3)
results[6,,] = simulation(1000,250,10,3,3)
results[7,,] = simulation(1000,300,10,3,3)
results[8,,] = simulation(1000,350,10,3,3)
results[9,,] = simulation(1000,400,10,3,3)
results[10,,] = simulation(1000,450,10,3,3)
results[11,,] = simulation(1000,500,10,3,3)

#F max error
pdf(file="FMax_T10.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results[,1,1]~c(25,seq(50,500,50)),type="l",ylim=c(0,4.2),col="red",lwd=3,xlab="",ylab="")
lines(results[,2,1]~c(25,seq(50,500,50)),col="blue",lty="dashed",lwd=3)
dev.off()

#F Frobenius error
pdf(file="FFrob_T10.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results[,1,2]~c(25,seq(50,500,50)),type="l",ylim=c(0,2),col="red",lwd=3,xlab="",ylab="")
lines(results[,2,2]~c(25,seq(50,500,50)),col="blue",lty="dashed",lwd=3)
dev.off()

#G max error
pdf(file="GMax_T10.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results[,1,3]~c(25,seq(50,500,50)),type="l",ylim=c(0,3.5),col="red",lwd=3,xlab="",ylab="")
lines(results[,2,3]~c(25,seq(50,500,50)),col="blue",lty="dashed",lwd=3)
dev.off()

#G Frobenius error
pdf(file="GFrob_T10.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results[,1,4]~c(25,seq(50,500,50)),type="l",ylim=c(0,1.5),col="red",lwd=3,xlab="",ylab="")
lines(results[,2,4]~c(25,seq(50,500,50)),col="blue",lty="dashed",lwd=3)
dev.off()



##2
#T=50 fixed
results2 = array(0,dim=c(11,2,4))

results2[1,,] = simulation(1000,25,50,3,3)
results2[2,,] = simulation(1000,50,50,3,3)
results2[3,,] = simulation(1000,100,50,3,3)
results2[4,,] = simulation(1000,150,50,3,3)
results2[5,,] = simulation(1000,200,50,3,3)
results2[6,,] = simulation(1000,250,50,3,3)
results2[7,,] = simulation(1000,300,50,3,3)
results2[8,,] = simulation(1000,350,50,3,3)
results2[9,,] = simulation(1000,400,50,3,3)
results2[10,,] = simulation(1000,450,50,3,3)
results2[11,,] = simulation(1000,500,50,3,3)

#F max error
pdf(file="FMax_T50.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results2[,1,1]~c(25,seq(50,500,50)),type="l",ylim=c(0,4.2),col="red",lwd=3,xlab="",ylab="")
lines(results2[,2,1]~c(25,seq(50,500,50)),col="blue",lty="dashed",lwd=3)
dev.off()

#F Frobenius error
pdf(file="FFrob_T50.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results2[,1,2]~c(25,seq(50,500,50)),type="l",ylim=c(0,2),col="red",lwd=3,xlab="",ylab="")
lines(results2[,2,2]~c(25,seq(50,500,50)),col="blue",lty="dashed",lwd=3)
dev.off()

#G max error
pdf(file="GMax_T50.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results2[,1,3]~c(25,seq(50,500,50)),type="l",ylim=c(0,3.5),col="red",lwd=3,xlab="",ylab="")
lines(results2[,2,3]~c(25,seq(50,500,50)),col="blue",lty="dashed",lwd=3)
dev.off()

#G Frobenius error
pdf(file="GFrob_T50.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results2[,1,4]~c(25,seq(50,500,50)),type="l",ylim=c(0,1.5),col="red",lwd=3,xlab="",ylab="")
lines(results2[,2,4]~c(25,seq(50,500,50)),col="blue",lty="dashed",lwd=3)
dev.off()



##3
#N=200 fixed
results3 = array(0,dim=c(10,2,4))

results3[1,,] = simulation(1000,200,10,3,3)
results3[2,,] = simulation(1000,200,20,3,3)
results3[3,,] = simulation(1000,200,30,3,3)
results3[4,,] = simulation(1000,200,40,3,3)
results3[5,,] = simulation(1000,200,50,3,3)
results3[6,,] = simulation(1000,200,60,3,3)
results3[7,,] = simulation(1000,200,70,3,3)
results3[8,,] = simulation(1000,200,80,3,3)
results3[9,,] = simulation(1000,200,90,3,3)
results3[10,,] = simulation(1000,200,100,3,3)

#F max error
pdf(file="FMax_N200.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results3[,1,1]~seq(10,100,10),type="l",ylim=c(0,5),col="red",lwd=3,xlab="",ylab="")
lines(results3[,2,1]~seq(10,100,10),col="blue",lty="dashed",lwd=3)
dev.off()

#F Frobenius error
pdf(file="FFrob_N200.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results3[,1,2]~seq(10,100,10),type="l",ylim=c(0,2),col="red",lwd=3,xlab="",ylab="")
lines(results3[,2,2]~seq(10,100,10),col="blue",lty="dashed",lwd=3)
dev.off()

#G max error
pdf(file="GMax_N200.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results3[,1,3]~seq(10,100,10),type="l",ylim=c(0,3.5),col="red",lwd=3,xlab="",ylab="")
lines(results3[,2,3]~seq(10,100,10),col="blue",lty="dashed",lwd=3)
dev.off()

#G Frobenius error
pdf(file="GFrob_N200.pdf",width=4,height=3)
par(mar=c(2,2,1,1))
plot(results3[,1,4]~seq(10,100,10),type="l",ylim=c(0,1.5),col="red",lwd=3,xlab="",ylab="")
lines(results3[,2,4]~seq(10,100,10),col="blue",lty="dashed",lwd=3)
dev.off()
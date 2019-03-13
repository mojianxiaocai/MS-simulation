setwd("H:/chengjie/TNAR simulation/MS-NAR simulation")
install.packages("markovchain")
source('MS-NAR_helper.R')
ptm<-proc.time()

#############Generate the Responses
set.seed(1992)
Nsize=100 #node number 
Time=50 #sample number 
beta=c(0.3,0.2,0.3,0.7) #parameter
p11=0.7   #Transition probability
p22=0.4   
Beta = beta[2:4]
sig=1
W = getDyadW(N = Nsize, N1 = 10, delta = -0.8, normalize = TRUE)
G1 = beta[2] * W + beta[3] * diag(1, nrow(W))
G2 = beta[2] * W + beta[4] * diag(1, nrow(W))
G=list(G1,G2)
gamma0 = c(-0.5, 0.3, 0.8, 0, 0)  ### true parameter for gamma
theta=c(beta,p11,p22,sig,sig,gamma0)
ZSigma = 0.5^abs(outer(1:5, 1:5, "-"))  ### Sigma_z (covariance for Z)
Z=mvrnorm(n=Nsize,mu=rep(0,nrow(ZSigma)), Sigma=ZSigma) ## Z~N(0,ZSigma)




#theta.start <- c(0.3,0.2,0.3,0.7,0.7,0.4,1,1,gamma0)
theta.start <- c(0.5,0.5,0.5,0.5,0.5,0.5,2,2,0.5,0.5,0.5,0.5,0.5)



Nrep=10
num <- matrix(rep(NA,Nrep*length(theta)),nrow=length(theta))
num2<-matrix(rep(NA,Nrep*length(theta)),nrow=length(theta))
for(k in 1:Nrep){
  Ymat=simu.Ymat(W, beta0 = beta[1], Beta = beta[2:4], Time = Time, G = G, Z = Z, sig = 1,Nsize = Nsize)
  max.lik.optim <- try(optim(theta.start, lik ,Ymat=Ymat, W=W,method= "BFGS", hessian=FALSE))
  num[,k]<-max.lik.optim$par
  num2[,k]<-(theta-max.lik.optim$par)^2
}

estimation<-(apply(num,1,mean))
error<-apply(num2,1,sum)/Nrep
CI=matrix(rep(NA,5*length(theta)),nrow=length(theta))
for( k in 1:length(theta)){
  CI[k,]<-quantile(num[k,],probs = c(0.05,0.25,0.5,0.75,0.95))
}
colnames(CI)<-c("0.05","0.25","0.5","0.75","0.95")
(result=rbind(estimation,error,t(CI)))
colnames(result)<-c("beta0","beta1","beta2_1","beta2_2","p11","p22","delta1","delta2","gamma1","gamma2","gamma3","gamma4","gamma5")
proc.time()-ptm

write.csv(result,"MS-Simulation_result_T50_N100.csv")


setwd("H:/chengjie/TNAR simulation/MS-NAR simulation/N500T100")
source('MS-NAR_helper.R')
input_fileName<-sapply('num_v', function(x){
  paste(x,1:10,".csv",sep='')
})
num=read.csv(input_fileName[1])[,-1]
for(i in 2:10){
  num=cbind(num,read.csv(input_fileName[i])[,-1])
}
num=as.matrix(num)
theta=c(0.3,0.2,0.3,0.7,0.7,0.4,1.0,1.0,-0.5,0.3,0.8,0.0,0.0)
num2<-matrix(rep(NA,ncol(num)*length(theta)),nrow=length(theta))
num3<-matrix(rep(NA,ncol(num)*length(theta)),nrow=length(theta))

for(k in 1:ncol(num)){
  num2[,k]<-(theta-num[,k])^2
  num3[,k]<-abs(theta-num[,k]) 
}

estimation<-(apply(num,1,mean))
error<-apply(num2,1,sum)/ncol(num)
(aae<-apply(num3,2,mean))
CI=matrix(rep(NA,5*length(theta)),nrow=length(theta))
for( k in 1:length(theta)){
  (CI[k,]=quantile(num[k,],probs = c(0.05,0.25,0.5,0.75,0.95)))
}
colnames(CI)<-c("0.05","0.25","0.5","0.75","0.95")
(result=rbind(estimation,error,t(CI)))
colnames(result)<-c("beta0","beta1","beta2_1","beta2_2","p11","p22","delta1","delta2","gamma1","gamma2","gamma3","gamma4","gamma5")
write.csv(result,"MS-Simulation_result_T100_N500.csv")
write.csv(aae,"T100_N500_aae.csv")

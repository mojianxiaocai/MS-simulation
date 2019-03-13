num<-num[,1:300]
num2<-num2[,1:300]

num<-num[,-which(abs(num[6,])>10)]
num2<-num2[,-which(num2[6,]>10)]
estimation<-(apply(num[,1:79],1,mean))
error<-apply(num2[,1:79],1,sum)/Nrep
CI=matrix(rep(NA,5*length(theta)),nrow=length(theta))
for( k in 1:length(theta)){
  CI[k,]<-quantile(num[k,1:79],probs = c(0.05,0.25,0.5,0.75,0.95))
}
colnames(CI)<-c("0.05","0.25","0.5","0.75","0.95")
(result=rbind(estimation,error,t(CI)))
colnames(result)<-c("beta0","beta1","beta2_1","beta2_2","p11","p22","delta1","delta2","gamma1","gamma2","gamma3","gamma4","gamma5")

write.csv(result,"MS-Simulation_result_T20_N300_v2.csv")
36387.51
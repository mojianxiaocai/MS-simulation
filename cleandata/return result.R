estimation<-(apply(num,1,mean))
error<-apply(num2,1,sum)/Nrep
(aae<-apply(num3,2,mean))
CI=matrix(rep(NA,5*length(theta)),nrow=length(theta))
for( k in 1:length(theta)){
  CI[k,]<-quantile(num[k,],probs = c(0.05,0.25,0.5,0.75,0.95))
}
colnames(CI)<-c("0.05","0.25","0.5","0.75","0.95")
(result=rbind(estimation,error,t(CI)))
colnames(result)<-c("beta0","beta1","beta2_1","beta2_2","p11","p22","delta1","delta2","gamma1","gamma2","gamma3","gamma4","gamma5")
write.csv(num,"num_v1.csv")


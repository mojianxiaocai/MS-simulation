for(k in 1:ncol(num)){
  if(num[3,k]>num[4,k]){
       sub=num[3,k]
       num[3,k]=num[4,k]
       num[4,k]=sub
     }
     if(num[5,k]<num[6,k]){
       sub=num[5,k]
       num[5,k]=num[6,k]
       num[6,k]=sub
     }
  num2[,k]<-(theta-num[,k])^2
  num3[,k]<-abs(theta-num[,k]) 
}
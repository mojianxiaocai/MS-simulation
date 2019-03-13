num
hist(num)

library(tidyverse)
num<-num[,-which(abs(num[6,])>10)]
num2<-num2[,-which(num2[5,]>10)]
num2<-num2[,-which(is.na(num2[1,]))]
tnum<-t(num2)
#mean(num2)
(aae<-apply(tnum,1,mean))

ggplot()
tnum<-cbind(tnum,aae)
colnames(tnum)<-c("beta0","beta1","beta2_1","beta2_2","p11","p22","delta1","delta2","gamma1","gamma2","gamma3","gamma4","gamma5","AAE")
data<-as.data.frame(tnum)

ggplot(data,aes(aae))+
  geom_density()
ggplot(data,aes(aae))+
  geom_histogram(binwidth = 0.001)

write.csv(aae,"T20_N100_aae.csv")

s=read.csv("T20_N100_aae.csv")
s[,1]=rep("T20N100",nrow(s))
a=read.csv("T20_N300_aae.csv")
a[,1]=rep("T20N300",nrow(a))
b=read.csv("T20_N400_aae.csv")
b[,1]=rep("T20N400",nrow(b))
c=read.csv("T50_N100_aae.csv")
c[,1]=rep("T50N100",nrow(c))
d=read.csv("T50_N300_aae.csv")
d[,1]=rep("T50N300",nrow(d))
e=read.csv("T50_N400_aae.csv")
e[,1]=rep("T50N400",nrow(e))
f=read.csv("T100_N100_aae.csv")
f[,1]=rep("T100N100",nrow(f))
g=read.csv("T100_N300_aae.csv")
g[,1]=rep("T100N300",nrow(g))
aae2=rbind(s,a,b,c,d,e,f,g)
colnames(aae2)<-c("size","aae")
aae=list(T20N300=a,T20N400=b,T50N100=c,T50N300=d,T50N400=e,T100N100=f,T100N300=g)
data=as.data.frame(aae)

ggplot(aae2,aes(size,aae,file=size))+
  geom_boxplot(alpha=0.5)

#####delete na column
delete.na<-function(num){
  num<-num[,-which(is.na(apply(num,2,sum)))]
  return(num)
}



num<-delete.na(num)
num2<-delete.na(num2)
num3<-delete.na(num3)
lik<-function(theta,Ymat,W){
  Ymat1=W%*%Ymat
  Time=ncol(Ymat)-1
  N=nrow(Ymat)
  X = cbind(rep(1, nrow(Ymat) * Time), as.vector(Ymat1[, -ncol(Ymat)]), 
            as.vector(Ymat[, -ncol(Ymat)]))
  Yvec = as.vector(Ymat[, -1])
  temp1=X%*%c(theta[-c(4,5:8)])
  temp2=X%*%c(theta[-c(3,5:8)])
  alpha1 <- theta[1]#beta0
  alpha2 <- theta[2]#beta1
  alpha3 <- theta[3]#beta2_1
  alpha4 <- theta[4]#beta2_2
  alpha5 <- theta[7] #delta_1
  alpha6 <- theta[8] #delta_2
  p11 <-theta[5]
  p22 <-theta[6]
  error1=(Yvec-temp1)
  error2=(Yvec-temp2)
  sigma1=diag(rep(alpha5^2,N),N)
  sigma2=diag(rep(alpha6^2,N),N)
  
  dist.1<-0
  for(i in 1:Time){
    dist.1<-rbind(dist.1,dmvnorm(error1[((i-1)*N+1):(i*N)],log = TRUE))
  }
  dist.2 <- 0
  for(i in 1:Time){
    dist.2<-rbind(dist.2,dmvnorm(error2[((i-1)*N+1):(i*N)],log = TRUE))
  }  
  dist <- cbind(dist.1[-1], dist.2[-1])
  o.v <- c(1,1)
  P <- matrix(c(p11, 1-p11, 1- p22, p22), nrow=2, ncol=2)
  nstates=2
  xi.a <- rep(0,2*Time)
  xi.a <- matrix(xi.a, nrow=Time,ncol=2)
  xi.b <- rep(0,2*Time)
  xi.b <- matrix(xi.b, nrow=Time,ncol=2)
  model.lik <- rep(0, Time)
  A=rbind(diag(nstates) - P, rep(1, nstates))
  xi.a[1,] <- (solve(t(A) %*% A) %*% t(A))[, nstates + 1]
  #xi.a[1,]<-(c(p11,p22)*dist[1,])/as.numeric(o.v%*%(c(p11,p22)*dist[1,]))
  for(i in 1:(Time-1)){
    xi.b[i+1,]<-P%*%xi.a[i,]
    xi.a[i+1,]<-(xi.b[i+1,]*dist[i+1,])/as.numeric(o.v%*%(xi.b[i+1,]*dist[i+1,]))
    model.lik[i+1]<-o.v%*%(xi.b[i+1,]*dist[i+1,])
  }
  logl <- sum(log(model.lik[2:length(model.lik)]))
  return(-logl)
}

theta.start <- c(0.3,0.2,0.3,0.7,0.7,0.4,1,1)
theta <- rep(0.5,length(theta.start))

max.lik.optim <- optim(theta.start, lik ,Ymat=Ymat, W=W,method= "BFGS", hessian=T)
max.lik.optim$par

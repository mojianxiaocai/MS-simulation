# lnoil x
# lnng y
lik <- function(theta,Ymat,W){
  Ymat1 = W %*% Ymat
  Time = ncol(Ymat) - 1
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
  dist.1 <- 0
  dist.1 <- (1/(alpha5*sqrt(2*pi)))*exp((-(Yvec-temp1)^2)/(2*alpha5^2))
  dist.2 <- 0
  dist.2 <- (1/(alpha6*sqrt(2*pi)))*exp((-(Yvec-temp2)^2)/(2*alpha6^2))
  dist <- cbind(dist.1, dist.2)
  o.v <- c(1,1)
  P <- matrix(c(p11, 1-p11, 1- p22, p22), nrow=2, ncol=2)
  xi.a <- rep(0,2*nrow(Ymat)*Time)
  xi.a <- matrix(xi.a, nrow=nrow(Ymat)*Time,ncol=2)
  xi.b <- rep(0,2*nrow(Ymat)*Time)
  xi.b <- matrix(xi.a, nrow=nrow(Ymat)*Time,ncol=2)
  model.lik <- rep(0, nrow(Ymat)*Time)
  xi.a[1,] <- (c(p11,p22)*dist[1,])/as.numeric(o.v%*%(c(p11,p22)*dist[1,]))
  ## Here is the Hamilton filter
  for (i in 1:(nrow(Ymat)*Time-1)){
    xi.b[i+1,] <- P%*%xi.a[i,]
    xi.a[i+1,] <- (xi.b[i+1,]*dist[i+1,])/as.numeric(o.v%*%(xi.b[i+1,]*dist[i+1,]))
    model.lik[i+1] <- o.v%*%(xi.b[i+1,]*dist[i+1,])
  }
  logl <- sum(log(model.lik[2:length(model.lik)]))
  return(-logl)
}

theta.start <- c(0.3,0.2,0.3,0.7,0.7,0.4,1,1)
theta.start2 <- rep(0.9,length(theta.start))

max.lik.optim <- optim(theta.start, lik ,Ymat=Ymat, W=W,method= "BFGS", hessian=T)
max.lik.optim$par
##max.lik.optim$par ## the estimation of p11 and p22 isnot good
##[1] 0.3902751 0.1550396 0.3145242 0.7027228 0.9839558 0.9838522 1.0157611 0.9645398
OI <- solve(max.lik.optim$hessian)

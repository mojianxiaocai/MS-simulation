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
  
  nlags <- 1
  nstates <- 2 ^ (nlags+1)
  lagstate <- 1 + outer( 1:nstates, 1:(nlags+1), FUN=function(i,j) { trunc((i-1) / 2 ^ (nlags + 1-j) ) %% 2 } )
  
  
  transit <- outer( X=1:nstates, Y=1:nstates, FUN=function(i,j) {
    (( 2 *lagstate[i,1] + lagstate[j,1] - 1) - 1) * (((i-1) %% (2^nlags)) == trunc((j-1)/2)) + 1
  } ) 
  
  tp <- c( 0, p11, 1-p22, 1-p11, p22)
  P <- array(tp[transit], c(nstates, nstates))
  A <- rbind( diag(nstates) - P, rep(1, nstates) )  # bottom of page 684
  ergodic.pi <- (solve( t(A) %*% A ) %*% t(A)) [,nstates + 1] # [22.2.26]
  
  xi.t.t <- ergodic.pi %o% rep(1,nlags)
  xi.t.t_1 <- xi.t.t
  log.likelihood <- 0
  for ( tt in 1:Time*nrow(Ymat) )
  {
  residuals <-c((Yvec-temp1)[tt], (Yvec-temp2)[tt]) #[22.4.24]
  eta.t <- dnorm(residuals, mean = 0, sd = sigma)       # [22.4.2 ]
  fp <- eta.t * xi.t.t_1[,tt]         # numerator [22.4.5]
  fpt <- sum(fp)                          # [22.4.8]
  xi.t.t <- cbind( xi.t.t, fp / fpt )                 # [22.4.5]
  log.likelihood <- log.likelihood + log(fpt)   # [22.4.7]
  xi.t.t_1 <- cbind( xi.t.t_1, P %*% xi.t.t[,tt] )                 # [22.4.6]
  }
  return(-log.likelihood)
}

theta.start <- c(0.3,0.2,0.3,0.7,0.7,0.4,1,1)
theta.start2 <- rep(0.9,length(theta.start))

max.lik.optim <- optim(theta.start, lik ,Ymat=Ymat, W=W,method= "BFGS", hessian=T,gr=NULL)
max.lik.optim$par
##max.lik.optim$par ## the estimation of p11 and p22 isnot good
##[1] 0.3902751 0.1550396 0.3145242 0.7027228 0.9839558 0.9838522 1.0157611 0.9645398
OI <- solve(max.lik.optim$hessian)
?optim

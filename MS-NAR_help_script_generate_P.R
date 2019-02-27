generate.P<-function(nlags=1){
  nstates <- 2 ^ (nlags+1)
  lagstate <- 1 + outer( 1:nstates, 1:(nlags+1), FUN=function(i,j) { trunc((i-1) / 2 ^ (nlags + 1-j) ) %% 2 } )
  transit <- outer( X=1:nstates, Y=1:nstates, FUN=function(i,j) {
    (( 2 *lagstate[i,1] + lagstate[j,1] - 1) - 1) * (((i-1) %% (2^nlags)) == trunc((j-1)/2)) + 1
  } ) 
  list(transit=transit,nstates=nstates)
}

nstates=generate.P(0)$nstates
lik.fun<-function(THETA,YT,W){
  p11star <- THETA[  "p11star" ]
  p22star <- THETA[  "p22star" ]
  tp <- c( 0, p11star, 1-p22star, 1-p11star, p22star )
  transit=generate.P(0)$transit
  nstates=generate.P(0)$nstates
  P <- array(tp[transit], c(nstates, nstates))
  A=rbind(diag(nstates) - P, rep(1, nstates))
  ergodic.pi <- (solve(t(A) %*% A,tol = 1e-19) %*% t(A))[, nstates + 1]
  
  xi.t.t<-rep(0,nstates*length(YT[,-ncol(YT)]))
  xi.t.t<-matrix(xi.t.t,nrow=length(YT[,-ncol(YT)]),ncol=nstates)
  xi.t.t_1<-rep(0,nstates*length(YT[,-ncol(YT)]))
  xi.t.t_1<-matrix(xi.t.t_1,nrow=length(YT[,-ncol(YT)]),ncol=nstates)
  model.lik<-rep(0,length(YT[,-ncol(YT)]))
  
  xi.t.t[1,]<-
  xi.t.t_1[1,]<-xi.t.t[1,]
  
  alpha5 <- THETA["sigma1"] #delta_1
  alpha6 <- THETA["sigma2"] #delta_2
  Ymat1 = W %*% YT
  Yvec = as.vector(YT[, -1])
  Time = ncol(YT) - 1
  X = cbind(rep(1, nrow(YT) * Time), as.vector(Ymat1[, -ncol(YT)]), 
            as.vector(YT[, -ncol(YT)]))
  temp1=X%*%c(THETA[c(3:5)])
  temp2=X%*%c(THETA[c(3,4,6)])
  dist.1 <- (1/(alpha5*sqrt(2*pi)))*exp((-(Yvec-temp1)^2)/(2*alpha5^2))
  dist.2 <- (1/(alpha6*sqrt(2*pi)))*exp((-(Yvec-temp1)^2)/(2*alpha6^2))
  dist <- cbind(dist.1, dist.2,dist.1, dist.2)
  
  for(i in 1:(length(YT[,-ncol(YT)])-1)){
    xi.t.t_1[i+1,]<-P%*%xi.t.t[i,]
    fp<-dist[i,]*xi.t.t_1[i,]
    fpt<-sum(fp)
    xi.t.t[i+1,]<-fp/fpt
    model.lik[i+1]<-fpt
  }
  logl<-sum(log(model.lik[2:length(model.lik)]))
  return(-logl)
}


THETA <- c( p11star=.85, p22star=.70, beta=c(0.4,0.4,0.4,0.9),
            sigma=c(1,1) )
max.lik.optim<- optim(par=THETA, hessian=TRUE, fn=lik.fun, gr=NULL, YT=Ymat,W=W,method="BFGS")
max.lik.optim$par
